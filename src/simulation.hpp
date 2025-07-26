#ifndef SIMULATION_H
#define SIMULATION_H

#include "geometric.hpp"
#include "glm.hpp"

#include <chrono>
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
	std::chrono::time_point<std::chrono::system_clock> last_update = std::chrono::system_clock::now();
};

class Simulation
{
	//private:
	public:
		int n_particles;
		float p_radius = 1.;
		particle_t *particles;
		particle_t *new_particles;
		domain_t domain;
		float speed = 1.;
		float last_delta_t = 0.;

		float particles_bounciness = .9;

		void spawn_particles_as_rect()
		{
			for(int i = 0; i < n_particles; i++)
			{
				particles[i] = particle_t();
			}

			float side_length = sqrt(n_particles);

			// aspect_ratio*side_length to deform the square of coordinates
			int n_x = floor(sqrt(n_particles)* (domain.size.x / domain.size.y));

			int n_y = floor(n_particles/n_x);

			// Number of particles in last row in case not a perfect square
			int n_x_last_row = mod((float)n_particles, (float)n_x);
			
			float rect_size = 0.9;
			float delta_x = domain.size.x / n_x * rect_size;
			float delta_y = domain.size.y / (n_y+1.) * rect_size;

			// Generate a complete rectangle of particles
			for(int y = 0; y < n_y; y++)
			{
				for(int x = 0; x < n_x; x++)
				{
					int i = y*n_x + x;

					particles[i].position = vec2(x*delta_x, y*delta_y) + vec2(domain.size.x * (1.-rect_size), domain.size.y * (1.-rect_size));
					particles[i].velocity = vec2(particles[i].position.x - domain.size.x / 2., particles[i].position.y - domain.size.y/2.);
				}
			}

			// If the number of particles doens't have a perfect square root,
			// add an incomplete last row to rectangle
			for(int x = 0; x < n_x_last_row; x++)
			{
				int i = n_y*n_x + x;

				particles[i].position = vec2(x*delta_x, (n_y-1)*delta_y) + vec2(domain.size.x * (1.-rect_size), domain.size.y * (1.-rect_size));
				particles[i].velocity = vec2(particles[i].position.x - domain.size.x / 2., particles[i].position.y - domain.size.y/2.);
			}
		}

		vec2 total_forces_on_particle(particle_t* p)
		{
			if(domain.radial_gravity)
			{
				vec2 to_center = vec2(domain.size.x / 2., domain.size.y / 2.) - p->position;
				float distance = glm::length(to_center);
				float gravity_norm = glm::length(domain.gravity);

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

				// Normalize
				normal /= dist;
				
				// Collision happens
				if(dist < 2.*p_radius)
				{
					//std::cout << "Collision between " << i << " and " << p_i << " at dist " << dist << std::endl;

					// A is inside B's radius, move it out
					vec2 a_newpos = A->position + vec2(-normal.x * (dist - 2.*p_radius), -normal.y * (dist - 2.*p_radius));

					// dir = vel + 2*normal*dot(-vel, normal)
					// Symmetrize (is that a word ?) the velocity vector
					// from the normal vector (-vel because they face
					// opposite from each other)
					float dot_a_vel_normal = dot(-A->velocity, normal);
					vec2 a_newdir = A->velocity + vec2(normal.x * 2. * dot_a_vel_normal, normal.y * 2. * dot_a_vel_normal);

					A->velocity = vec2(
							a_newdir.x * particles_bounciness,
							a_newdir.y * particles_bounciness
					);
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
			this->particles = (particle_t*)malloc(sizeof(particle_t)*n_particles);
			this->new_particles = (particle_t*)malloc(sizeof(particle_t)*n_particles);
			this->n_particles = n_particles;

			this->domain = domain_t{.size = domain_size};

			spawn_particles_as_rect();
		}

		void tick()
		{
			// One weird bug, delta_t increases dramatically when the mouse goes over the top window bar
			float delta_t = (std::chrono::system_clock::now() - domain.last_update).count() / 1e9f;
			last_delta_t = delta_t;
			delta_t *= speed;

			for(int p_i = 0; p_i < n_particles; p_i++)
			{
				particle_t* n_p = &new_particles[p_i];
				particle_t* p = &particles[p_i];

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

				collision_check(n_p, p_i);
				domain_boundaries(n_p);

				n_p->velocity += total_forces_on_particle(p) * delta_t;		//  v(t-1) + a(t)*t
				n_p->position += n_p->velocity*delta_t; 					//  p(t-1) + v(t)*t
			}


			// Swap both particle buffers
			particle_t* last_frame_data = particles;
			particles = new_particles;
			new_particles = last_frame_data;

			domain.last_update = std::chrono::system_clock::now();
		}

		~Simulation()
		{
			free(particles);
		}
};

#endif
