#ifndef SIMULATION_H
#define SIMULATION_H

#include "geometric.hpp"
#include "glm.hpp"

#include <chrono>
#include <iostream>

using namespace glm;

#define p_radius 1.
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
	vec2 gravity = vec2(0., 0./*-9.81*/);
	float bounciness = 1.;
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

		float particles_bounciness = 1.;

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

		void collision_check(particle_t* p, int i)
		{
			for(int p_i = 0; p_i < n_particles; p_i++)
			{
				if(i == p_i) continue;

				particle_t* p_current = &particles[p_i];

				// Early continue with cheap test
				float delta_x = p->position.x - p_current->position.x;
				float delta_y = p->position.y - p_current->position.y;
				if(delta_x > 2.*p_radius || delta_x < -2.*p_radius || delta_y > 2.*p_radius || delta_y < -2.*p_radius) continue;

				// Vector going from p_current's center to p's center,
				// normal to the surface of both
				vec2 normal = p->position - p_current->position;

				float dist = glm::length(normal);

				// Normalize
				normal /= dist;
				
				// Collision happens
				if(dist <= 2.*p_radius && dist >= -2.*p_radius)
				{
					float v_a = glm::length(p->velocity);
					float v_b = glm::length(p_current->velocity);

					float m_a = p->mass;
					float m_b = p_current->mass;

					// p = m*v
					float momentum_a = v_a * m_a;
					float momentum_b = v_b * m_b;

					float total_momentum = momentum_a + momentum_b;

					// k = 1/2 * m * v²
					float kinetic_a = 1./2. * m_a * v_a*v_a;
					float kinetic_b = 1./2. * m_b * v_b*v_b;

					float total_kinetic = kinetic_a + kinetic_b;

					// `normal` var points from p_current's center to p's center
					vec2 normal_a = -normal;
					//vec2 normal_b = -normal;

					// glm doesn't provide scalar <-> vector multiplication
					// Dot product should be cheap but i just store it to use twice
					float vel_a_norm_dot = glm::dot(p->velocity, normal_a);
					vec2 dir_a = p->velocity + vec2(normal_a.x * 2.*vel_a_norm_dot, normal_a.y * 2.*vel_a_norm_dot);
					dir_a = glm::normalize(dir_a);

					float a = -2. * m_b / (m_a+m_b);
					float dot_a = glm::dot(p->velocity - p_current->velocity, p->position - p_current->position);
					float dist_squared = dist*dist;
					vec2 diff = vec2(normal.x*dist, normal.y*dist);

					float fact = a*dot_a/dist_squared;
					vec2 vel_a = vec2(diff.x*fact, diff.y*fact);


					//float vel_b_norm_dot = glm::dot(p_current->velocity, normal_b);
					//vec2 dir_b = p_current->velocity + vec2(normal_b.x * 2.*vel_b_norm_dot, normal_b.y * 2.*vel_b_norm_dot);
					//dir_b = glm::normalize(dir_b);
					
					float vel_a_final = v_a * (m_a-m_b) / (m_a+m_b)  +  v_b * 2.*m_b / (m_a+m_b);
					float vel_b_final = v_b * (m_a-m_b) / (m_a+m_b)  +  v_a * 2.*m_a / (m_a+m_b);

					//p->velocity = vec2(dir_a.x * vel_a_final, dir_a.y * vel_a_final);
					p->velocity = vel_a;

					printf("\n\nBounce between %d and %d. Current velocity: (%.1f, %.1f): %.1f", i, p_i, p->velocity.x, p->velocity.y, glm::length(v_a));
					// Particle is inside another, move it away
					p->position += vec2(normal.x*(2.*p_radius-dist), normal.y*(2.*p_radius-dist));

					//printf("           New velocity: (%.1f, %.1f): %.1f\n", p->velocity.x, p->velocity.y, v_a_final);
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
			//delta_t = 0.01f;
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

				n_p->velocity += total_forces_on_particle(p) * delta_t;	//  v(t-1) + a(t)*t
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
