#ifndef SIMULATION_H
#define SIMULATION_H

#include "geometric.hpp"
#include "glm.hpp"

#include <chrono>
#include <thread>
#include <iostream>

using namespace glm;

struct particle_t
{
	vec2 position = vec2(0., 0.);
	vec2 velocity = vec2(0., 0.);
	float mass = 1.;
};

struct domain_t
{
	vec2 size = vec2(1., 1.);
	bool radial_gravity = false;
	vec2 gravity = vec2(0., -9.81);
	float bounciness = .9;
	std::chrono::time_point<std::chrono::high_resolution_clock> last_update = std::chrono::high_resolution_clock::now();
};

class Simulation
{
	private:
		std::thread t;
		unsigned int iteration = 0;

	public:
		int n_particles;
		float p_radius = 1.;
		particle_t *particles;
		particle_t *new_particles;
		domain_t domain;
		float speed = 1.;
		//std::chrono::duration<long int, std::ratio<1, 1000000000> > last_delta_t = std::chrono::milliseconds(0);
		float last_delta_t = 0.;
		unsigned int n_threads = 4;
		unsigned int sim_hertz = 300;
		bool should_run = true;

		float particles_bounciness = .9;

		void spawn_particles_as_rect()
		{
			float side_length = ceil(sqrt(n_particles));
			float dx = domain.size.x/side_length;
			float dy = domain.size.y/side_length;

			std::cout << "N_particles: " << this->n_particles << std::endl;

			for(int i = 0; i < n_particles; i++)
			{
				particle_t* p = &particles[i];

				*p = particle_t();

				p->position = vec2(
					mod((float)i, side_length) * dx,
					(int)(i / side_length) * dy
				);


				p->velocity = vec2(
						p->position.x - domain.size.x / 2.,
						p->position.y - domain.size.y / 2.
				);

				p->mass = 1.;

				new_particles[i] = particles[i];
			}
		}

		vec2 total_forces_on_particle(particle_t* p)
		{
			if(domain.radial_gravity)
			{
				vec2 to_center = vec2(domain.size.x / 2., domain.size.y / 2.) - p->position;
				float distance = glm::length(to_center);
				float gravity_norm = glm::length(domain.gravity);

				if(distance == 0.) distance = 0.01;

				to_center.x *= gravity_norm / distance;
				to_center.y *= gravity_norm / distance;

				return to_center;
			}

			return domain.gravity;
		}

		void collision_check(particle_t* A, int i)
		{
			for(int p_i = 0; p_i < n_particles; p_i++)
			{
				if(i == p_i) continue;

				particle_t* B = &particles[p_i];

				// Early continue with cheap test
				float delta_x = A->position.x - B->position.x;
				float delta_y = A->position.y - B->position.y;
				if(delta_x > 2.*p_radius || delta_x < -2.*p_radius || delta_y > 2.*p_radius || delta_y < -2.*p_radius) continue;

				// Vector going from B's center to A's center,
				// normal to the surface of both
				vec2 normal = A->position - B->position;

				float dist = glm::length(normal);
				if(dist == 0.) dist = 0.01;

				// Normalize
				normal /= dist;
				
				// Collision happens
				if(dist < 2.*p_radius)
				{
					vec2 relative_vel = B->velocity - A->velocity;

					// Bounce direction
					vec2 a_newdir = A->velocity + normal*dot(relative_vel, normal);

					// Dampen the bounce
					A->velocity = a_newdir*particles_bounciness;

					// A is inside B's radius, move it out
					vec2 a_newpos = A->position - normal * (dist - 2.f*p_radius);
					A->position = a_newpos;
				}
			}
		}

		void domain_boundaries(particle_t* p)
		{
			// Floor & Ceiling
			if(p->position.y < 1.)
			{
				p->velocity.y = -p->velocity.y * domain.bounciness;
				p->position.y = 1.;
			}
			if(p->position.y > domain.size.y - 1)
			{
				p->velocity.y = -p->velocity.y * domain.bounciness;
				p->position.y = domain.size.y - 1.;
			}

			// Walls
			if(p->position.x < 1.)
			{
				p->velocity.x = -p->velocity.x * domain.bounciness;
				p->position.x = 1.;
			}
			if(p->position.x > domain.size.x - 1.)
			{
				p->velocity.x = -p->velocity.x * domain.bounciness;
				p->position.x = domain.size.x - 1.;
			}
		}


	//public:
		Simulation(int n_particles, vec2 domain_size)
		{
			// Allocate needed space for particles
			this->particles 	= (particle_t*)malloc(sizeof(particle_t)*n_particles);
			this->new_particles = (particle_t*)malloc(sizeof(particle_t)*n_particles);
			this->n_particles = n_particles;

			this->domain = domain_t{.size = domain_size};

			spawn_particles_as_rect();
		}

		void run()
		{
			if(t.joinable())
			{
				// Stop infinite loop
				if(should_run) should_run = false;
				t.join();
			}

			// Seems like ok behavior to set it to true since the intent
			// behind calling run() should be to, well, run the sim
			should_run = true;

			t = std::thread([&]()
			{
				while(should_run)
				{
					auto sleep_time = 1.e9 * std::chrono::nanoseconds(1) / ((float)sim_hertz);

					if(std::chrono::high_resolution_clock::now() - sleep_time >= domain.last_update)
						tick();
				}
			});
		}

		void tick()
		{
			auto now = std::chrono::high_resolution_clock::now();

			float delta_t = (now - domain.last_update).count() / 1.e9; // count gives nanoseconds, *1e9 to get seconds
			last_delta_t = delta_t;
			delta_t *= speed; // Scale delta_t depending on sim speed

			domain.last_update = now;

			// Saved number of used threads to leave none behind when cleaning up
			unsigned int used_threads = n_threads;
			std::thread ts[used_threads];

			// Ceil to get every particles in case it's not round
			int p_per_thread = ceil((float)n_particles/(float)used_threads);

			// For every thread
			// Gotta pass t_i by value, else it'll continue changing while the thread is running
			for(int t_i = 0; t_i < used_threads; t_i++) ts[t_i] = std::thread([&, t_i]()
			{
				// Iterate on a given list of particles
				for(int p_i = t_i*p_per_thread; p_i < (t_i+1)*p_per_thread && p_i < n_particles; p_i++)
				{
					/*
					 * If the order stays the same, some particles
					 * will have collision checks before later ones,
					 * which will move them always in the same order
					 * and can lead to those being stuck.
					 * Every two steps, invert the iteration order
					 */
					int i = p_i;
					if(mod((float)iteration, 2.f) == 0.) i = n_particles - p_i - 1;

					particle_t* n_p = &new_particles[i];
					particle_t* p = &particles[i];

					n_p->position = p->position;
					n_p->velocity = p->velocity;

					// Assuming gravity is the only force acting on particles:
					// The following formula gives the analytic solution to the
					// particle's position
					// a(t) = g
					// v(t) = V0 + g*t
					// p(t) = P0 + V0*t + g/2*t²
					//
					// To take into account force changes and collision, the sim cannot
					// use the analytic solution, so using the formula above, it gives:
					// v(t) = v(t-1) + a(t)*t
					// p(t) = p(t-1) + v(t)*t
					//
					// a(t) is known on time but cannot be predicted (without running the sim in advance)
					// a(t) in the first line and v(t) in the second are treated as if they
					// are constant in the interval (t, t+1). The higher the tickrate, the 
					// more accurate the sim will be as this interval will be reduced

					collision_check(n_p, i);
					domain_boundaries(n_p);

					n_p->velocity += total_forces_on_particle(p) * delta_t;		//  v(t-1) + a(t)*t
					n_p->position += n_p->velocity * delta_t; 					//  p(t-1) + v(t)*t
				}
			});

			// End threads
			for(int t_i= 0; t_i < used_threads; t_i++)
				ts[t_i].join();

			// Swap both particle buffers
			particle_t* last_frame_data = particles;
			particles = new_particles;
			new_particles = last_frame_data;

			iteration++;
		}

		~Simulation()
		{
			if(t.joinable())
			{
				should_run = false;
				t.join();
			}

			free(particles);
			free(new_particles);
		}
};

#endif
