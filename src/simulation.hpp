#ifndef SIMULATION_H
#define SIMULATION_H

#include "glm.hpp"

#include <chrono>
#include <thread>

#include "memory.h"

using glm::vec2, glm::mod;

struct particle_t
{
	vec2 position = vec2(1., 1.);
	vec2 velocity = vec2(0., 0.);
	vec2 acceleration = vec2(0., 0.);
	float mass = 1.;
};

struct domain_t
{
	vec2 size = vec2(500., 500.);
	bool radial_gravity = false;
	vec2 gravity = vec2(0., -9.81);
	float bounciness = .9;
};

struct sim_state_t
{
	unsigned int n_threads;

	// Number of ticks
	unsigned int iteration = 0;
	
	// One day I'll understand how std::chrono works
	// with the ratios and type annotations
	//std::chrono::duration<>();
	float sim_seconds = 0.; // Number of seconds since sim start

	float delta_t = 0.;
	std::chrono::time_point<std::chrono::high_resolution_clock> last_update = std::chrono::high_resolution_clock::now();

	domain_t domain;

	particle_t* particles;
	particle_t* buff_particles;
	unsigned int n_particles;
	float p_radius = 1.;
	float p_bounciness = 0.9;
};

struct bounding_boxes_t
{
	unsigned int nbx = 0;
	unsigned int nby = 0;

	// Flattened 2D array of arrays of pointers to particles
	particle_t* * * p_per_b = NULL; 
	
	// Flattened 2D array containing the
	// number of particles per bounding box
	int* n_p_per_b = NULL; 
};

struct settings_t
{
	bool run = false;
	unsigned int n_threads = 4;
	unsigned int hertz = 300;

	float speed = 1.;

	unsigned int n_bounding_boxes_x = 1;
	unsigned int n_bounding_boxes_y = 1;

	vec2 domain_size = vec2(500., 500.);
	float domain_bounciness = 0.9;
	vec2 domain_gravity = vec2(0., -9.81);
	bool domain_gravity_radial = false;

	float particle_radius = 1.;
	float particles_bounciness = 0.9;
};

class Simulation
{
	private:
		std::thread t;
		bounding_boxes_t bboxes;
		sim_state_t state;

	public:
		settings_t settings;

		// Acts as public access to particles data
		particle_t* particles;

		void update_settings()
		{
			state.domain.size 			= settings.domain_size;
			state.domain.bounciness 	= settings.domain_bounciness;
			state.domain.gravity 		= settings.domain_gravity;
			state.domain.radial_gravity = settings.domain_gravity_radial;

			state.p_radius 		= settings.particle_radius;
			state.p_bounciness 	= settings.particles_bounciness;

			// Clamp the number of bounding boxes to avoid
			// bounding boxes being of size < 2*p_radius which
			// would cause problems for collision detection
			int max_b_x = (int)floor(state.domain.size.x/(state.p_radius*2.));
            int max_b_y = (int)floor(state.domain.size.y/(state.p_radius*2.));
			if(settings.n_bounding_boxes_x > max_b_x)
				settings.n_bounding_boxes_x = max_b_x;
			if(settings.n_bounding_boxes_y > max_b_y)
				settings.n_bounding_boxes_y = max_b_y;

			// Clamp to min 1
			if(settings.n_bounding_boxes_x < 1) settings.n_bounding_boxes_x = 1;
			if(settings.n_bounding_boxes_y < 1) settings.n_bounding_boxes_y = 1;

			// Clamp n_threads to max supported threads by cpu
			unsigned int max_threads = std::thread::hardware_concurrency();
			if(settings.n_threads > max_threads)
				settings.n_threads = max_threads;

			state.n_threads = settings.n_threads;
		}

		const unsigned int& n_particles()
			{return state.n_particles;}

		const float& sim_tick_time()
			{return state.delta_t;}

		void spawn_particles_as_rect()
		{
			float side_length = ceil(sqrt(state.n_particles));
			float dx = state.domain.size.x/side_length;
			float dy = state.domain.size.y/side_length;

			for(int i = 0; i < state.n_particles; i++)
			{
				particle_t* p = &state.particles[i];

				*p = particle_t();

				p->position = vec2(
					mod((float)i, side_length) * dx + state.p_radius,
					(int)(i / side_length) * dy + state.p_radius
				);


				p->velocity = vec2(
						p->position.x - state.domain.size.x / 2.,
						p->position.y - state.domain.size.y / 2.
				);

				p->mass = 1.;

				state.buff_particles[i] = *p;
				particles[i] = *p;
			}
		}

		void total_forces_on_particle(particle_t* p)
		{
			if(state.domain.radial_gravity)
			{
				vec2 to_center = state.domain.size/2.f - p->position;

				float distance = glm::length(to_center);
				float gravity_norm = glm::length(state.domain.gravity);

				if(distance == 0.) distance = 0.01;
				to_center *= gravity_norm / distance;

				p->acceleration += to_center;
				return;
			}

			p->acceleration += state.domain.gravity;
		}

		void collision_check(particle_t* A, particle_t* old_A)
		{
			// Position in bounding boxes
			int A_x = (int)(A->position.x/state.domain.size.x * bboxes.nbx);
			int A_y = (int)(A->position.y/state.domain.size.y * bboxes.nby);
			int A_xy = A_x*bboxes.nby+ A_y;

			int B_x = A_x - 1;
			if(B_x < 0) B_x = 0;
			for(;B_x <= A_x+1 && B_x < bboxes.nbx; B_x++)
			{

				int B_y = A_y - 1;
				if(B_y < 0) B_y = 0;
				for(;B_y <= A_y+1 && B_y < bboxes.nby; B_y++)
				{
					int B_xy = B_x*bboxes.nby+ B_y;

					for(int p_i = 0; p_i < bboxes.n_p_per_b[B_xy]; p_i++)
					{
						particle_t* B = bboxes.p_per_b[B_xy][p_i];

						// Early continue with cheap test
						float delta_x = A->position.x - B->position.x;
						float delta_y = A->position.y - B->position.y;
						if(	delta_x >  2.*state.p_radius ||
							delta_x < -2.*state.p_radius ||
							delta_y >  2.*state.p_radius ||
							delta_y < -2.*state.p_radius) continue;

						// Continue if the test is with itself
						if(B == old_A) continue;

						// Vector going from B's center to A's center,
						// normal to the surface of both
						vec2 normal = A->position - B->position;

						// Squared distance for the test
						float dist = normal.x*normal.x + normal.y*normal.y;
						if(dist == 0.) dist = 0.01;
						
						// Collision happens, test with (2*p_radius)² since
						// dist is the distance squared
						if(dist < 4.*state.p_radius*state.p_radius)
						{
							// The actual distance is needed,
							// no choice about using square root
							dist = sqrt(dist);

							// Normalize
							normal /= dist;

							// Mirror the relative velocity (B_vel - A_vel) from the normal vector and dampen the bounce
							A->velocity += normal * (dot(B->velocity - A->velocity, normal)*state.p_bounciness);

							// A is inside B's radius, move it out
							A->position -= normal * (dist - 2.f * state.p_radius);
						}
					}
				}
			}
		}

		void domain_boundaries(particle_t* p)
		{
			// Floor & Ceiling
			if(p->position.y < state.p_radius)
			{
				//if(p->acceleration.y < 0.) p->acceleration.y = -p->acceleration.y;

				p->velocity.y = -p->velocity.y * state.domain.bounciness;
				p->position.y = state.p_radius;
			}
			if(p->position.y > state.domain.size.y - state.p_radius)
			{
				//if(p->acceleration.y > 0.) p->acceleration.y = -p->acceleration.y;

				p->velocity.y = -p->velocity.y * state.domain.bounciness;
				p->position.y = state.domain.size.y - state.p_radius;
			}

			// Walls
			if(p->position.x < state.p_radius)
			{
				//if(p->acceleration.x < 0.) p->acceleration.x = -p->acceleration.x;

				p->velocity.x = -p->velocity.x * state.domain.bounciness;
				p->position.x = state.p_radius;
			}
			if(p->position.x > state.domain.size.x - state.p_radius)
			{
				//if(p->acceleration.x > 0.) p->acceleration.x = -p->acceleration.x;

				p->velocity.x = -p->velocity.x * state.domain.bounciness;
				p->position.x = state.domain.size.x - state.p_radius;
			}
		}


		Simulation(int n_particles)
		{
			// Allocate needed space for particles
			particles 			= (particle_t*)malloc(sizeof(particle_t)*n_particles);
			state.particles 	= (particle_t*)malloc(sizeof(particle_t)*n_particles);
			state.buff_particles= (particle_t*)malloc(sizeof(particle_t)*n_particles);
			state.n_particles 	= n_particles;

			spawn_particles_as_rect();
		}

		void reset(int n_particles)
		{
			std::clog << "Resetting sim with " << n_particles << " particles" << std::endl;

			bool was_running = settings.run;
			end();

			particles 			= (particle_t*)realloc(particles, 				sizeof(particle_t)*n_particles);
			state.particles 	= (particle_t*)realloc(state.particles, 		sizeof(particle_t)*n_particles);
			state.buff_particles= (particle_t*)realloc(state.buff_particles, 	sizeof(particle_t)*n_particles);

			state.n_particles 	= n_particles;

			spawn_particles_as_rect();

			if(was_running) run();
		}

		void run()
		{
			if(t.joinable())
			{
				// Stop infinite loop
				settings.run = false;
				t.join();
			}

			// Seems like ok behavior to set it to true since the intent
			// behind calling run() should be to, well, run the sim
			settings.run = true;

			// Reset timing to avoid end() <-> run() time difference issues
			state.last_update = std::chrono::high_resolution_clock::now();

			t = std::thread([&]()
			{
				while(settings.run)
				{
					auto sleep_time = 1.e9 * std::chrono::nanoseconds(1) / ((float)settings.hertz);

					// Sleep until it's time to tick() again, based on settings.hertz
					std::this_thread::sleep_until(state.last_update + sleep_time);
					tick();
				}
			});
		}

		void end()
		{
			settings.run = false;
			if(t.joinable()) t.join();
		}

		void build_bounding_boxes()
		{
			// Free particles lists before creating new one.
			// Since malloc wasn't called for lists of length
			// 0, don't call free for those
			for(int i = 0; i < bboxes.nbx*bboxes.nby; i++)
				if(bboxes.n_p_per_b[i] > 0) free(bboxes.p_per_b[i]);

			bool has_changed =
				bboxes.nbx != settings.n_bounding_boxes_x ||
				bboxes.nby != settings.n_bounding_boxes_y;

			bboxes.nbx = settings.n_bounding_boxes_x;
			bboxes.nby = settings.n_bounding_boxes_y;

			// Realloc space, bounding boxes sizes have changed
			if(has_changed)
			{
				bboxes.n_p_per_b = (int*)realloc(bboxes.n_p_per_b,
						sizeof(int)*bboxes.nbx*bboxes.nby);

				bboxes.p_per_b   = (particle_t***)realloc(bboxes.p_per_b,
						sizeof(particle_t**)*bboxes.nbx*bboxes.nby);
			}

			// Incremental list, counting current index of particle when inserting
			// (see last loop inserting pointers)
			unsigned int n_p_per_b_again[bboxes.nbx*bboxes.nby];

			// Assigning all 0s to 2D array containing the number
			// of particles per bounding box
			memset(bboxes.n_p_per_b, 0, sizeof(unsigned int)*bboxes.nbx*bboxes.nby);
			memset(n_p_per_b_again , 0, sizeof(unsigned int)*bboxes.nbx*bboxes.nby);

			// Counting number of particles per bounding box
			for(int p_i = 0; p_i < state.n_particles; p_i++)
			{
				particle_t* p = &state.particles[p_i];

				int x = (int)(p->position.x/state.domain.size.x * bboxes.nbx);
				int y = (int)(p->position.y/state.domain.size.y * bboxes.nby);
				int xy = x*bboxes.nby + y;
				
				// Some particles might be leaking out of the domain
				if(x < 0 || y < 0 || x > bboxes.nbx-1 || y > bboxes.nby-1) continue;

				bboxes.n_p_per_b[xy]++;
			}

			// Alloc space needed in list of particles per bounding box.
			// Only call malloc if lists contain particles
			for(int i = 0; i < bboxes.nbx*bboxes.nby; i++)
				if(bboxes.n_p_per_b[i] > 0) bboxes.p_per_b[i] = (particle_t**)malloc(sizeof(particle_t*)*bboxes.n_p_per_b[i]);

			// Insert particles pointers in bounding boxes
			for(int p_i = 0; p_i < state.n_particles; p_i++)
			{
				particle_t* p = &state.particles[p_i];

				int p_x = (int)(p->position.x/state.domain.size.x * bboxes.nbx);
				int p_y = (int)(p->position.y/state.domain.size.y * bboxes.nby);
				int p_xy = p_x*bboxes.nby + p_y;

				if(p_x < 0 || p_y < 0 || p_x > bboxes.nbx-1 || p_y > bboxes.nby-1) continue;

				// Count current particle index in bounding box
				bboxes.p_per_b[p_xy][n_p_per_b_again[p_xy]++] = p;
			}
		}

		void tick()
		{
			auto now = std::chrono::high_resolution_clock::now();
			state.delta_t = (now - state.last_update).count() / 1.e9; 	// count gives nanoseconds, *1e9 to get seconds
			float delta_t = state.delta_t * settings.speed; 			// Scale this tick's delta_t depending on sim speed
			state.last_update = now;

			update_settings();
			build_bounding_boxes();

			std::thread ts[state.n_threads];

			// Ceil to get every particles in case it's not round
			int p_per_thread = ceil((float)state.n_particles/(float)state.n_threads);

			// For every thread
			// Gotta pass t_i by value, else it'll continue changing while the thread is running
			for(int t_i = 0; t_i < state.n_threads; t_i++) ts[t_i] = std::thread([&, t_i]()
			{
				// Iterate on a given list of particles
				for(int p_i = t_i*p_per_thread; p_i < (t_i+1)*p_per_thread && p_i < state.n_particles; p_i++)
				{
					/*
					 * If the order stays the same, some particles
					 * will have collision checks before later ones,
					 * which will move them always in the same order
					 * and can lead to those being stuck.
					 * Every two steps, invert the iteration order
					 * to try to reduce this effect
					 */
					int i = p_i;
					if(mod((float)state.iteration, 2.f) < 1.)
						i = state.n_particles - p_i - 1;

					particle_t* n_p = &state.buff_particles[i];
					particle_t* p = &state.particles[i];

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

					total_forces_on_particle(n_p);
					collision_check(n_p, p);
					domain_boundaries(n_p);

					n_p->velocity += n_p->acceleration * delta_t;		//  v(t-1) + a(t)*t
					n_p->position += n_p->velocity * delta_t; 					//  p(t-1) + v(t)*t
					
					particles[i] = *p;
				}
			});

			// End threads
			for(int t_i= 0; t_i < state.n_threads; t_i++)
				ts[t_i].join();

			// Swap both particle buffers
			particle_t* last_frame_data = state.particles;
			state.particles = state.buff_particles;
			state.buff_particles = last_frame_data;

			state.iteration++;
		}

		~Simulation()
		{
			settings.run = false;

			// Stop main thread
			if(t.joinable())
				t.join();

			// Free bounding boxes with particles
			for(int i = 0; i < bboxes.nbx*bboxes.nby; i++)
					if(bboxes.n_p_per_b[i] > 0) free(bboxes.p_per_b[i]);
			free(bboxes.p_per_b);

			// Free number of particles per bounding box
			free(bboxes.n_p_per_b);

			free(particles);
			free(state.particles);
			free(state.buff_particles);
		}
};

#endif
