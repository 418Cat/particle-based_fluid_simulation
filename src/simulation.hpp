#ifndef SIMULATION_H
#define SIMULATION_H

#include "glm.hpp"

#include <chrono>
#include <iostream>

using namespace glm;

struct particle_t
{
	vec2 position = vec2(0., 0.);
	vec2 velocity = vec2(0., 0.);
};

struct domain_t
{
	vec2 size = vec2(1., 1.);
	vec2 gravity = vec2(0., 9.81);
	std::chrono::time_point<std::chrono::system_clock> last_update = std::chrono::system_clock::now();
};

class Simulation
{
	private:
		int n_particles;
		particle_t *particles;
		domain_t domain;

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
			
			float delta_x = domain.size.x / n_x;
			float delta_y = domain.size.y / (n_y+1.);

			// Draw a complete rectangle of particles
			for(int y = 0; y < n_y; y++)
			{
				for(int x = 0; x < n_x; x++)
				{
					int i = y*n_x + x;

					particles[i].position = vec2(x*delta_x, y*delta_y);
					particles[i].velocity = vec2(0.);
				}
			}

			// If the number of particles doens't have a perfect square root,
			// add an incomplete last row to rectangle
			for(int x = 0; x < n_x_last_row; x++)
			{
				int i = n_y*n_x + x;

				particles[i].position = vec2(x*delta_x, (n_y-1)*delta_y);
				particles[i].velocity = vec2(0.);
			}
		}


	public:
		Simulation(int n_particles, vec2 domain_size)
		{
			// Allocate needed space for particles
			this->particles = (particle_t*)malloc(sizeof(particle_t)*n_particles);
			this->n_particles = n_particles;

			this->domain = domain_t{.size = domain_size};

			spawn_particles_as_rect();
		}

		void tick()
		{
			float delta_t = (std::chrono::system_clock::now() - domain.last_update).count() / 1e9f;

			for(int p_i = 0; p_i < n_particles; p_i++)
			{
				particle_t* p = &particles[p_i];

				// Assuming gravity is the only force acting on particles:
				// a(t) = g
				// v(t) = V0 + g*t
				// p(t) = P0 + V0*t + g/2*t²
				
				p->position += 
					p->velocity*delta_t + 								//  V0*t + 
					vec2(domain.gravity.x/2., domain.gravity.y/2.) *	//  g/2 *
					(delta_t*delta_t);									//  t²

				p->velocity += domain.gravity * delta_t;				//  g*t
			}

			domain.last_update = std::chrono::system_clock::now();
		}

		~Simulation()
		{
			free(particles);
		}
};

#endif
