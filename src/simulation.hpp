#ifndef SIMULATION_H
#define SIMULATION_H

#include "glm.hpp"

#include <chrono>
#include <thread>
#include <iostream>

using namespace glm;

struct particle_t
{
	vec2 position = vec2(0., 0.);
	vec2 velocity = vec2(0., 0.);
	vec2 acceleration = vec2(0., 0.);
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
		particle_t *new_particles;
		unsigned int nbx = 0;
		unsigned int nby = 0;

	public:
		particle_t *particles;
		domain_t domain;

		float last_delta_t = 0.;
		int n_particles;

		float p_radius = 1.;
		float speed = 1.;

		unsigned int n_threads = 4;
		unsigned int sim_hertz = 300;
		bool should_run = true;

		float particles_bounciness = .9;

		unsigned int n_bounding_boxes_x = 40;
		unsigned int n_bounding_boxes_y = 40;
		particle_t* * * p_per_b; // Flattened 2D array of arrays of pointers to particles
		int* n_p_per_b; // Flattened 2D array

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

		void total_forces_on_particle(particle_t* p)
		{
			if(domain.radial_gravity)
			{
				vec2 to_center = domain.size/2.f - p->position;

				float distance = glm::length(to_center);
				float gravity_norm = glm::length(domain.gravity);

				if(distance == 0.) distance = 0.01;
				to_center *= gravity_norm / distance;

				p->acceleration += to_center;
				return;
			}

			p->acceleration += domain.gravity;
		}

		void collision_check(particle_t* A, particle_t* old_A)
		{
			// Position in bounding boxes
			int A_x = (int)floor(A->position.x/domain.size.x * nbx);
			int A_y = (int)floor(A->position.y/domain.size.y * nby);
			int A_xy = A_x*nby+ A_y;

			int B_x = A_x - 1;
			if(B_x < 0) B_x = 0;
			for(;B_x <= A_x+1 && B_x < nbx; B_x++)
			{

				int B_y = A_y - 1;
				if(B_y < 0) B_y = 0;
				for(;B_y <= A_y+1 && B_y < nby; B_y++)
				{
					int B_xy = B_x*nby+ B_y;
					for(int p_i = 0; p_i < n_p_per_b[B_xy]; p_i++)
					{
						particle_t* B = p_per_b[B_xy][p_i];

						// Early continue with cheap test
						float delta_x = A->position.x - B->position.x;
						float delta_y = A->position.y - B->position.y;
						if(delta_x > 2.*p_radius || delta_x < -2.*p_radius || delta_y > 2.*p_radius || delta_y < -2.*p_radius) continue;

						// Continue if the test is with itself
						if(B == old_A) continue;

						
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
							vec2 a_newpos = A->position - normal * (dist - 2.f * p_radius);
							A->position = a_newpos;
						}
					}
				}
			}
		}

		void domain_boundaries(particle_t* p)
		{
			// Floor & Ceiling
			if(p->position.y < p_radius)
			{
				// Equal and opposite reaction blablabla
				if(p->acceleration.y < 0.) p->acceleration.y = 0.;

				p->velocity.y = -p->velocity.y * domain.bounciness;
				p->position.y = p_radius;
			}
			if(p->position.y > domain.size.y - p_radius)
			{
				if(p->acceleration.y > 0.) p->acceleration.y = 0.;

				p->velocity.y = -p->velocity.y * domain.bounciness;
				p->position.y = domain.size.y - p_radius;
			}

			// Walls
			if(p->position.x < p_radius)
			{
				if(p->acceleration.x < 0.) p->acceleration.x = 0.;

				p->velocity.x = -p->velocity.x * domain.bounciness;
				p->position.x = p_radius;
			}
			if(p->position.x > domain.size.x - p_radius)
			{
				if(p->acceleration.x > 0.) p->acceleration.x = 0.;

				p->velocity.x = -p->velocity.x * domain.bounciness;
				p->position.x = domain.size.x - p_radius;
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
				should_run = false;
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

		void end()
		{
			should_run = false;
			if(t.joinable()) t.join();
		}

		void build_bounding_boxes()
		{
			// Free particles lists before creating new one
			for(int i = 0; i < nbx*nby; i++) free(p_per_b[i]);

			bool has_changed = nbx != n_bounding_boxes_x || nby != n_bounding_boxes_y;

			nbx = n_bounding_boxes_x;
			nby = n_bounding_boxes_y;

			// Realloc space, bounding boxes sizes have changed
			if(has_changed)
			{
				n_p_per_b = (int*)realloc(n_p_per_b, sizeof(int)*nbx*nby);
				p_per_b = (particle_t***)realloc(p_per_b, sizeof(particle_t**)*nbx*nby);
			}

			// Incremental list, counting current index of particle when inserting
			// (see last loop inserting pointers)
			unsigned int n_p_per_b_again[nbx*nby];

			// Assigning all 0s to 2D array containing the number
			// of particles per bounding box
			for(int i = 0; i < nbx*nby; i++)
			{
				n_p_per_b[i] = 0;
				n_p_per_b_again[i] = 0;
			}

			// Counting number of particles per bounding box
			for(int p_i = 0; p_i < n_particles; p_i++)
			{
				particle_t* p = &particles[p_i];

				int x = (int)floor(p->position.x/domain.size.x * nbx);
				int y = (int)floor(p->position.y/domain.size.y * nby);

				// Some particles might be leaking out of the domain
				if(x < 0 || y < 0 || x > nbx-1 || y > nby-1) continue;

				n_p_per_b[x*nby+ y]++;
			}

			// Alloc space needed in list of particles per bounding box
			for(int i = 0; i < nbx*nby; i++)
				p_per_b[i] = (particle_t**)malloc(sizeof(particle_t*)*n_p_per_b[i]);

			// Insert particles pointers in bounding boxes
			for(int p_i = 0; p_i < n_particles; p_i++)
			{
				particle_t* p = &particles[p_i];

				int p_x = (int)floor(p->position.x/domain.size.x * nbx);
				int p_y = (int)floor(p->position.y/domain.size.y * nby);
				int p_xy = p_x*nby + p_y;

				if(p_x < 0 || p_y < 0 || p_x > nbx-1 || p_y > nby-1) continue;
				
				// Count current particle index in bounding box
				p_per_b[p_xy][n_p_per_b_again[p_xy]++] = p;
			}
		}

		void tick()
		{
			auto now = std::chrono::high_resolution_clock::now();

			float delta_t = (now - domain.last_update).count() / 1.e9; // count gives nanoseconds, *1e9 to get seconds
			last_delta_t = delta_t;
			delta_t *= speed; // Scale delta_t depending on sim speed

			domain.last_update = now;

			int max_b_x = (int)floor(domain.size.x/(p_radius*2.));
            int max_b_y = (int)floor(domain.size.y/(p_radius*2.));
			if(n_bounding_boxes_x > max_b_x) n_bounding_boxes_x = max_b_x;
			if(n_bounding_boxes_y > max_b_y) n_bounding_boxes_y = max_b_y;
			build_bounding_boxes();

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

					n_p->mass = p->mass;
					n_p->position = p->position;
					n_p->velocity = p->velocity;
					n_p->acceleration = vec2(0., 0.); // Reset acceleration each frame

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

					collision_check(n_p, p);
					total_forces_on_particle(n_p);
					domain_boundaries(n_p);

					n_p->velocity += n_p->acceleration * delta_t;		//  v(t-1) + a(t)*t
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
			should_run = false;

			// Stop main thread
			if(t.joinable())
				t.join();

			// Free bounding boxes with particles
			for(int i = 0; i < nbx*nby; i++)
			{
					free(p_per_b[i]);
			}
			free(p_per_b);

			// Free number of particles per box
			free(n_p_per_b);

			free(particles);
			free(new_particles);
		}
};

#endif
