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
};

struct domain_t
{
	vec2 size = vec2(1., 1.);
	vec2 gravity = vec2(0., -9.81);
	std::chrono::time_point<std::chrono::system_clock> last_update = std::chrono::system_clock::now();
};

class Simulation
{
	//private:
	public:
		int n_particles;
		particle_t *particles;
		particle_t *new_particles;
		domain_t domain;
		float speed = 1.;

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

			// Draw a complete rectangle of particles
			for(int y = 0; y < n_y; y++)
			{
				for(int x = 0; x < n_x; x++)
				{
					int i = y*n_x + x;

					particles[i].position = vec2(x*delta_x, y*delta_y) + vec2(domain.size.x * (1.-rect_size), domain.size.y * (1.-rect_size));
					particles[i].velocity = particles[i].position - vec2(domain.size.x / 2., 0.);
				}
			}

			// If the number of particles doens't have a perfect square root,
			// add an incomplete last row to rectangle
			for(int x = 0; x < n_x_last_row; x++)
			{
				int i = n_y*n_x + x;

				particles[i].position = vec2(x*delta_x, (n_y-1)*delta_y) + vec2(domain.size.x * (1.-rect_size), domain.size.y * (1.-rect_size));
				particles[i].velocity = particles[i].position - vec2(domain.size.x / 2., 0.);
			}
		}

		vec2 total_forces_on_particle(int p_i)
		{
			return domain.gravity;
		}

		void collision_check(particle_t* p, int i)
		{
			for(int p_i = 0; p_i < n_particles; p_i++)
			{
				if(i == p_i) continue;

				particle_t* p_current = &particles[p_i];

				// Vector going from p_current's center to p's center,
				// normal to the surface of both
				vec2 normal = p->position - p_current->position;
				vec2 normalized_normal = glm::normalize(normal);

				float pyth = glm::length(normal) / 2.;

				// Collision happens
				if(pyth < 1. && pyth > -1.)
				{
					p->position += vec2(normalized_normal.x*2. - normal.x, normalized_normal.y*2. - normal.y);

					float norm_vel_dot = glm::dot(normalized_normal, p->velocity);

					vec2 bounce =
						vec2(p->velocity.x, p->velocity.y) +
						vec2(normalized_normal.x*-2.*norm_vel_dot, normalized_normal.y*-2.*norm_vel_dot);

					p->velocity = bounce;
					p->velocity.x *= .9;
					p->velocity.y *= .9;
				}
			}
		}

		void domain_boundaries(particle_t* p)
		{
			// Floor & Ceiling
			if(p->position.y < 1.)
			{
				p->velocity.y = -p->velocity.y * 0.9;
				p->position.y = 1.;
			}
			if(p->position.y > domain.size.y - 1)
			{
				p->velocity.y = -p->velocity.y * 0.9;
				p->position.y = domain.size.y - 1.;
			}

			// Walls
			if(p->position.x < 1.)
			{
				p->velocity.x = -p->velocity.x * .9;
				p->position.x = 1.;
			}
			if(p->position.x > domain.size.x - 1.)
			{
				p->velocity.x = -p->velocity.x * .9;
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

			float delta_t = (std::chrono::system_clock::now() - domain.last_update).count() / 1e9f;
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

				n_p->velocity += total_forces_on_particle(p_i) * delta_t;	//  v(t-1) + a(t)*t
				n_p->position += n_p->velocity*delta_t; 					//  p(t-1) + v(t)*t

				collision_check(n_p, p_i);
				domain_boundaries(n_p);

				p->position = n_p->position;
				p->velocity = n_p->velocity;
			}

			domain.last_update = std::chrono::system_clock::now();
		}

		~Simulation()
		{
			free(particles);
		}
};

#endif
