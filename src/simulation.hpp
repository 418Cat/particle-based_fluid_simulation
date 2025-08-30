#ifndef SIMULATION_H
#define SIMULATION_H

#include <cstdlib>
#include <iostream>

#include "ext/scalar_constants.hpp"
#include "glm.hpp"

#include <chrono>
#include <thread>

#include "memory.h"

using glm::mod, glm::dvec3, glm::ivec3, glm::vec;

// Easy precision change
#define num double
#define vec3 dvec3

struct particle_t
{
	vec3 position 		= vec3(1.);
	vec3 velocity 		= vec3(0.);
	vec3 acceleration 	= vec3(0.);
	num mass 			= 1.;

	ivec3 bbox_xyz		= glm::vec<3, int>(0);

	num density 		= 1.;
	num pressure 		= 1.;
};

struct tree_node_t
{
	num mass = 0.;
	unsigned int n_p = 0;
	particle_t** contained_p = NULL;// Must only contain particles
									// if node is LEAF or is being
									// built

	vec3 center_of_mass = vec3(0.);

	vec3 coords = vec3(0.);
	vec3 size = vec3(0.);
	num volume = 0.;

	tree_node_t* parent = NULL;
	tree_node_t* children = NULL; // Array of size 8

	enum type_t
	{
		ROOT,	// Doesn't have parents but has children
		BRANCH, // Has both 
		LEAF,	// Doesn't have children
		EMPTY, 	// Has no particles (So no children)
	} type;
};

struct octree_state_t
{
	unsigned int max_depth = 5;
	unsigned int min_p_in_node = 4;

	// (Gotta find a better name)
	// Threshold above which particles_gravity()
	// goes down one level deeper in the tree, based
	// on the ratio between the volume of the node and
	// the distance of the node's center of mass
	// from the particle
	num vol_dit_thresh = 0.01;
	num gravity_factor = 1.e10;
	tree_node_t* root;
};

struct domain_t
{
	vec3 size = vec3(500., 500., 500.);
	bool radial_gravity = false;
	bool gravity_axis[3] = {false, true, false};
	vec3 gravity = vec3(0., -9.81, 0.);
	num bounciness = .9;
};

struct sim_state_t
{
	unsigned int n_threads;

	num delta_t = 0.;
	std::chrono::time_point<std::chrono::high_resolution_clock> last_update = std::chrono::high_resolution_clock::now();

	domain_t domain;

	particle_t* particles;
	particle_t* buff_particles;
	unsigned int n_particles;

	bool p_collisions;
	enum p_collision_type_t 
	{
		VELOCITY,
		ACCELERATION
	} p_collision_type;
	num p_radius;
	num p_bounciness;

	bool p_gravity;
	bool p_gravity_inverse;

	bool liquid = true;
	num liquid_rest_density = 1.;
	num smoothing_length_h = 2.;
	enum liquid_kernel_f
	{
		CUBIC_SPLINE,
	} p_liquid_kernel = sim_state_t::CUBIC_SPLINE;
	num stiffness = 1.;
};

struct bounding_boxes_t
{
	unsigned int nbx = 0;
	unsigned int nby = 0;
	unsigned int nbz = 0;

	// Flattened 3D array of arrays of pointers to particles
	particle_t* * * p_per_b = NULL; 
	
	// Flattened 3D array containing the
	// number of particles per bounding box
	int* n_p_per_b = NULL; 
};

struct settings_t
{
	bool run = false;
	unsigned int n_threads = 4;
	unsigned int hertz = 1500;

	float speed = 1.;

	unsigned int n_bounding_boxes_x = 20;
	unsigned int n_bounding_boxes_y = 20;
	unsigned int n_bounding_boxes_z = 20;

	vec3 domain_size = vec3(200., 200., 200.);
	num domain_bounciness = 0.85;
	vec3 domain_gravity = vec3(0., -9.81, 0.);
	bool domain_gravity_axis[3] = {false, true, false};
	bool domain_gravity_radial = false;

	bool particles_collisions = true;
	sim_state_t::p_collision_type_t collision_type = sim_state_t::VELOCITY;
	num particle_radius = 1.;
	num particles_bounciness = 0.99;
	bool particle_gravity = false;
	num particles_gravity_factor = 1.e10;
	bool particles_gravity_inverse = false;

	unsigned int octree_max_depth = 7;
	num volume_to_distance_threshold = 0.5;
	unsigned int octree_min_particles_in_node = 4;

	bool liquid_sim = true;
	num liquid_rest_density = .5;
	num smoothing_length_h = 6.;
	sim_state_t::liquid_kernel_f liquid_kernel = sim_state_t::CUBIC_SPLINE;
	num stiffness_constant_k = 0.5;
};

class Simulation
{
	private:
		std::thread t;
		std::thread *ts;
		sim_state_t state;

		bounding_boxes_t bboxes;
		octree_state_t octree_state;

	public:
		settings_t settings;

		// Acts as public access to particles data
		particle_t* particles;
		tree_node_t* octree_root;

		void update_settings()
		{
			state.domain.size 			= settings.domain_size;
			state.domain.bounciness 	= settings.domain_bounciness;
			state.domain.gravity 		= settings.domain_gravity;
			for(int i = 0; i < 3; i++) state.domain.gravity_axis[i] = settings.domain_gravity_axis[i];
			state.domain.radial_gravity = settings.domain_gravity_radial;

			state.p_collisions = settings.particles_collisions;
			state.p_collision_type = settings.collision_type;
			state.p_radius 		= settings.particle_radius;
			state.p_bounciness 	= settings.particles_bounciness;

			// Clamp the number of bounding boxes to avoid
			// bounding boxes being of size < 2*p_radius which
			// would cause problems for collision detection
			int max_b_x = (int)floor(state.domain.size.x/(state.p_radius*2.));
            int max_b_y = (int)floor(state.domain.size.y/(state.p_radius*2.));
            int max_b_z = (int)floor(state.domain.size.z/(state.p_radius*2.));
			if(settings.n_bounding_boxes_x > max_b_x)
				settings.n_bounding_boxes_x = max_b_x;
			if(settings.n_bounding_boxes_y > max_b_y)
				settings.n_bounding_boxes_y = max_b_y;
			if(settings.n_bounding_boxes_z > max_b_z)
				settings.n_bounding_boxes_z = max_b_z;

			// Clamp to min 1
			if(settings.n_bounding_boxes_x < 1) settings.n_bounding_boxes_x = 1;
			if(settings.n_bounding_boxes_y < 1) settings.n_bounding_boxes_y = 1;
			if(settings.n_bounding_boxes_z < 1) settings.n_bounding_boxes_z = 1;

			// Clamp n_threads to max supported threads by cpu
			unsigned int max_threads = std::thread::hardware_concurrency();
			if(settings.n_threads > max_threads)
				settings.n_threads = max_threads;

			state.n_threads = settings.n_threads;

			state.p_gravity = settings.particle_gravity;
			state.p_gravity_inverse = settings.particles_gravity_inverse;

			octree_state.gravity_factor = settings.particles_gravity_factor;
			octree_state.max_depth = settings.octree_max_depth;

			octree_state.vol_dit_thresh = settings.volume_to_distance_threshold;
			octree_state.min_p_in_node = settings.octree_min_particles_in_node;

			state.liquid = settings.liquid_sim;
			state.liquid_rest_density = settings.liquid_rest_density;
			if(state.liquid_rest_density <= 0.) state.liquid_rest_density = 1.;
			state.smoothing_length_h = settings.smoothing_length_h;
			if(state.smoothing_length_h	<= 0.) state.smoothing_length_h = 1.;
			state.p_liquid_kernel = settings.liquid_kernel;
			state.stiffness = settings.stiffness_constant_k;
			if(state.stiffness <= 0.) state.stiffness = 1.;
		}

		const unsigned int& n_particles()
			{return state.n_particles;}

		const num& sim_tick_time()
			{return state.delta_t;}

		void spawn_particles_as_rect(bool with_vel=true)
		{
			num side_length = ceil(powf(state.n_particles, 1./3.));
			num dx = state.domain.size.x/side_length;
			num dy = state.domain.size.y/side_length;
			num dz = state.domain.size.z/side_length;

			for(int i = 0; i < state.n_particles; i++)
			{
				particle_t* p = &state.particles[i];

				*p = particle_t();

				p->position = vec3(
					mod((num)i, side_length) * dx + state.p_radius,
					(int)mod(i/side_length, side_length) * dy + state.p_radius,
					(int)(i/(side_length*side_length)) * dz + state.p_radius
				);

				if(with_vel) p->velocity = p->position - state.domain.size / 2.;
				p->mass = 1.;

				state.buff_particles[i] = *p;
				particles[i] = *p;
			}
		}

		void domain_interactions(particle_t* n_p, particle_t* p)
		{
			if(state.domain.radial_gravity)
			{
				vec3 to_center = state.domain.size/2. - n_p->position;

				num dist = distance(p->position, state.domain.size/2.);

				num gravity_norm = sqrt(
					(state.domain.gravity_axis[0] ? glm::abs(state.domain.gravity.x) : 0.) +
					(state.domain.gravity_axis[1] ? glm::abs(state.domain.gravity.y) : 0.) +
					(state.domain.gravity_axis[2] ? glm::abs(state.domain.gravity.z) : 0.)
				);

				if(dist == 0.) dist = 0.01;
				to_center *= gravity_norm / dist;

				n_p->acceleration += to_center / n_p->mass;
			}
			else
				n_p->acceleration += vec3(
					state.domain.gravity_axis[0] ? state.domain.gravity.x : 0.,
					state.domain.gravity_axis[1] ? state.domain.gravity.y : 0.,
					state.domain.gravity_axis[2] ? state.domain.gravity.z : 0.
				);


			// Sphere shaped domain ______________________________________________
			//double radius = 50.;

			//double dist_x = p->position.x - state.domain.size.x/2.;
			//double dist_y = p->position.y - state.domain.size.y/2.;
			//double dist_z = p->position.z - state.domain.size.z/2.;

			//double dist_sqrd = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

			//if(dist_sqrd > radius*radius)
			//{
			//	double dist = sqrt(dist_sqrd);

			//	vec3 normal = p->position - state.domain.size/2.;
			//	normal /= dist;

			//	p->position = state.domain.size/2. + normal * radius;

			//	double dot_pos_norm = dot(normal, p->velocity);
			//	p->velocity -= 2.*normal*dot_pos_norm;
			//	p->velocity *= state.domain.bounciness;
			//}
			// ___________________________________________________________________

			// For every component of the position, check if out of domain
			for(int i = 0; i < 3; i++)
			{
				num* p_pos   = &((num*)&(n_p->position))[i];
				num* p_vel   = &((num*)&(n_p->velocity))[i];
				num* p_accel = &((num*)&(n_p->acceleration))[i];

				num* w_coord = &((num*)&(state.domain.size))[i];

				// Outer wall check
				if(*p_pos >  *w_coord - state.p_radius)
				{
					//*p_accel = 0.;
					*p_vel  *= -state.domain.bounciness;
					*p_pos   = *w_coord - state.p_radius;
				}

				if(*p_pos < state.p_radius)
				{
					//*p_accel = 0.;
					*p_vel  *= -state.domain.bounciness;
					*p_pos   = state.p_radius;
				}
			}
		}

		void particles_collision(particle_t* A, particle_t* old_A)
		{
			// Position in bounding boxes
			int A_x = A->bbox_xyz.x;
			int A_y = A->bbox_xyz.y;
			int A_z = A->bbox_xyz.z;
			int A_xyz = A_x*bboxes.nby*bboxes.nbz + A_y*bboxes.nbz + A_z;

			int B_x = A_x - 1;
			if(B_x < 0) B_x = 0;
			for(;B_x <= A_x+1 && B_x < bboxes.nbx; B_x++)
			{

				int B_y = A_y - 1;
				if(B_y < 0) B_y = 0;
				for(;B_y <= A_y+1 && B_y < bboxes.nby; B_y++)
				{

					int B_z = A_z - 1;
					if(B_z < 0) B_z = 0;
					for(;B_z <= A_z+1 && B_z < bboxes.nbz; B_z++)
					{

						int B_xyz = B_x*bboxes.nby*bboxes.nbz + B_y*bboxes.nbz + B_z;

						for(int p_i = 0; p_i < bboxes.n_p_per_b[B_xyz]; p_i++)
						{
							particle_t* B = bboxes.p_per_b[B_xyz][p_i];

							// Early continue with cheap test
							num delta_x = glm::abs(A->position.x - B->position.x);
							num delta_y = glm::abs(A->position.y - B->position.y);
							num delta_z = glm::abs(A->position.z - B->position.z);
							if(	delta_x >  2.*state.p_radius ||
								delta_y >  2.*state.p_radius ||
								delta_z >  2.*state.p_radius) continue;

							// Continue if the test is with itself
							if(B == old_A) continue;

							// Vector going from B's center to A's center,
							// normal to the surface of both
							vec3 normal = A->position - B->position;

							// Squared distance for the test to avoid expensive sqrt
							num dist_sqrd = normal.x*normal.x + normal.y*normal.y + normal.z*normal.z;
							if(dist_sqrd == 0.) dist_sqrd = 0.01;
							
							// Collision happens, test with (2*p_radius)²
							if(dist_sqrd < 4.*state.p_radius*state.p_radius)
							{
								// The actual distance is needed
								num dist = sqrt(dist_sqrd);

								// Normalize
								normal /= dist;

								// Mirror the relative velocity (B_vel - A_vel) from 
								// the normal vector and dampen the bounce
								vec3 delta_vel = -normal * dot(A->velocity - B->velocity, normal) * state.p_bounciness;

								if(state.p_collision_type == sim_state_t::VELOCITY)
								{
									A->velocity += delta_vel;
								}
								else if(state.p_collision_type == sim_state_t::ACCELERATION)
								{
									// Add the derivative to the acceleration
									A->acceleration += delta_vel / state.delta_t;
								}

								// A is inside B's radius, move it out
								A->position += -normal*(dist - 2. * state.p_radius);

								// If A too close from B, move in a random direction
								// TODO: Make this vector always be of length 1
								// without using expensive operations like sqrt
								// other TODO: Check if rand() is expensive
								//if(dist <= state.p_radius * 1.1)
									//A->acceleration += -normal*(dist - 2. * state.p_radius) / 2.;

								if(dist <= state.p_radius*0.2)
									A->position += vec3(
											(num)std::rand() / RAND_MAX * state.p_radius,
											(num)std::rand() / RAND_MAX * state.p_radius,
											(num)std::rand() / RAND_MAX * state.p_radius
									);
							}
						}
					}
				}
			}
		}

		vec3 particles_gravity(particle_t* A, tree_node_t* N, unsigned int depth=0)
		{
			if(N->type == tree_node_t::EMPTY) return vec3(0.);

			const num G = 6.6743e-11 * octree_state.gravity_factor *
				(state.p_gravity_inverse ? -1. : 1.);

			vec3 normal = N->center_of_mass - A->position;
			num dist = length(normal);

			// Node is far enough, use center of mass and total mass
			if(N->volume / dist < octree_state.vol_dit_thresh)
			{
				return normal * G * (A->mass*N->mass)/(dist*dist);
			}

			vec3 total_gravity = vec3(0.);

			// Node doesn't have any children, iterate particles
			if(N->type == tree_node_t::LEAF)
			{
				for(int p_i = 0; p_i < N->n_p; p_i++)
				{
					particle_t* B = N->contained_p[p_i];

					if(A == B) continue;

					normal = B->position - A->position;
					dist = length(normal);

					total_gravity += normal * G * (A->mass*B->mass)/(dist*dist);
				}

				return total_gravity;
			}

			for(int n_i = 0; n_i < 8; n_i++)
			{
				total_gravity += particles_gravity(A, &(N->children[n_i]), depth+1);
			}

			return total_gravity;
		}

		void compute_liquid_density(particle_t* I, particle_t* old_I)
		{
			num density = 0.;


			// Position in bounding boxes
			int I_x = I->bbox_xyz.x;
			int I_y = I->bbox_xyz.y;
			int I_z = I->bbox_xyz.z;
			int I_xyz = I_x*bboxes.nby*bboxes.nbz + I_y*bboxes.nbz + I_z;

			// Ratio between bboxes' size and liquid_h to avoid resizing bboxes
			int h_to_bbox_x = ceil(state.smoothing_length_h / (state.domain.size.x / bboxes.nbx));
			int h_to_bbox_y = ceil(state.smoothing_length_h / (state.domain.size.y / bboxes.nby));
			int h_to_bbox_z = ceil(state.smoothing_length_h / (state.domain.size.z / bboxes.nbz));

			int J_x = I_x - h_to_bbox_x;
			if(J_x < 0) J_x = 0;
			for(;J_x <= I_x + h_to_bbox_x && J_x < bboxes.nbx; J_x++)
			{

				int J_y = I_y - h_to_bbox_y;
				if(J_y < 0) J_y = 0;
				for(;J_y <= I_y+h_to_bbox_y && J_y < bboxes.nby; J_y++)
				{

					int J_z = I_z - h_to_bbox_z;
					if(J_z < 0) J_z = 0;
					for(;J_z <= I_z+h_to_bbox_z && J_z < bboxes.nbz; J_z++)
					{

						int J_xyz = J_x*bboxes.nby*bboxes.nbz + J_y*bboxes.nbz + J_z;

						for(int p_i = 0; p_i < bboxes.n_p_per_b[J_xyz]; p_i++)
						{
							particle_t* J = bboxes.p_per_b[J_xyz][p_i];

							// Early continue with cheap test
							num delta_x = glm::abs(I->position.x - J->position.x);
							num delta_y = glm::abs(I->position.y - J->position.y);
							num delta_z = glm::abs(I->position.z - J->position.z);
							if(	delta_x >  2.*state.smoothing_length_h ||
								delta_y >  2.*state.smoothing_length_h ||
								delta_z >  2.*state.smoothing_length_h) continue;

							// Continue if the test is with itself
							//if(J == old_I) continue;

							// Vector going from J's center to I's center,
							// normal to the surface of both
							vec3 normal = I->position - J->position;

							// Squared distance for the test to avoid expensive sqrt
							num dist_sqrd = normal.x*normal.x + normal.y*normal.y + normal.z*normal.z;
							if(dist_sqrd == 0.) dist_sqrd = 0.01;
							
							// J is in radius h of I
							if(dist_sqrd < state.smoothing_length_h*state.smoothing_length_h)
							{
								// The actual distance is needed
								num dist = sqrt(dist_sqrd);

								// Normalize
								normal /= dist;

								// Wij = 1/(h^d) * f(q)
								// with q=dist/h
								// and d the number of dimensions
								num q = dist / state.smoothing_length_h;

								num f_q = 48./22. * (
									q < 1. ? 	2./3. - q*q+ 1./2. * q*q*q	:
									(q < 2. ? 	1./6. * pow(2. - q, 3.) 	:  0.
									));

								num w_ij = 1. / (pow(state.smoothing_length_h, 3)) * f_q;

								density += w_ij * J->mass;
							}
						}
					}
				}
			}

			// p_i = k * ( (P_i / P_0)^7  -  1)
			// with p_i the pressure exerted by current particle
			// P_i the density at location of current particle
			// P_0 the rest density of fluid
			// k a stiffness constant
			num pressure = state.stiffness * (
					pow(density / state.liquid_rest_density, 7.) - 1.
			);

			I->density = density;
			I->pressure = pressure;
		}

		void compute_liquid_pressure_viscosity(particle_t* I, particle_t* old_I)
		{
			
			vec3 sum_pressure = vec3(0.);
			vec3 sum_viscosity = vec3(0.);

			// Position in bounding boxes
			int I_x = I->bbox_xyz.x;
			int I_y = I->bbox_xyz.y;
			int I_z = I->bbox_xyz.z;
			int I_xyz = I_x*bboxes.nby*bboxes.nbz + I_y*bboxes.nbz + I_z;

			// Ratio between bboxes' size and liquid_h to avoid resizing bboxes
			int h_to_bbox_x = ceil(state.smoothing_length_h / (state.domain.size.x / bboxes.nbx));
			int h_to_bbox_y = ceil(state.smoothing_length_h / (state.domain.size.y / bboxes.nby));
			int h_to_bbox_z = ceil(state.smoothing_length_h / (state.domain.size.z / bboxes.nbz));

			int J_x = I_x - h_to_bbox_x;
			if(J_x < 0) J_x = 0;
			for(;J_x <= I_x + h_to_bbox_x && J_x < bboxes.nbx; J_x++)
			{

				int J_y = I_y - h_to_bbox_y;
				if(J_y < 0) J_y = 0;
				for(;J_y <= I_y+h_to_bbox_y && J_y < bboxes.nby; J_y++)
				{

					int J_z = I_z - h_to_bbox_z;
					if(J_z < 0) J_z = 0;
					for(;J_z <= I_z+h_to_bbox_z && J_z < bboxes.nbz; J_z++)
					{

						int J_xyz = J_x*bboxes.nby*bboxes.nbz + J_y*bboxes.nbz + J_z;

						for(int p_i = 0; p_i < bboxes.n_p_per_b[J_xyz]; p_i++)
						{
							particle_t* J = bboxes.p_per_b[J_xyz][p_i];

							// Early continue with cheap test
							num delta_x = glm::abs(I->position.x - J->position.x);
							num delta_y = glm::abs(I->position.y - J->position.y);
							num delta_z = glm::abs(I->position.z - J->position.z);
							if(	delta_x >  2.*state.smoothing_length_h ||
								delta_y >  2.*state.smoothing_length_h ||
								delta_z >  2.*state.smoothing_length_h) continue;

							// Continue if the test is with itself
							if(J == old_I) continue;

							// Vector going from J's center to I's center,
							// normal to the surface of both
							vec3 normal = I->position - J->position;

							// Squared distance for the test to avoid expensive sqrt
							num dist_sqrd = normal.x*normal.x + normal.y*normal.y + normal.z*normal.z;
							if(dist_sqrd == 0.) dist_sqrd = 0.01;
							
							// J is in radius h of I
							if(dist_sqrd < state.smoothing_length_h*state.smoothing_length_h)
							{
								// The actual distance is needed
								num dist = sqrt(dist_sqrd);

								// Normalize
								normal /= dist;

								vec3 total_force = vec3(0.);

								auto w_ij = [&](vec3 xyz)
								{
									// Wij = 1/(h^d) * f(q)
									// with q=dist/h
									// and d the number of dimensions
									num q = distance(I->position, xyz) / state.smoothing_length_h;

									num f_q = 48./22. * (
										q < 1. ? 	2./3. - q*q+ 1./2. * q*q*q	:
										(q < 2. ? 	1./6. * pow(2. - q, 3.) 	:  0.
										));

									return 1. / (pow(state.smoothing_length_h, 3)) * f_q;
								};

								const num delta = 0.05*state.smoothing_length_h;
								vec3 w_ij_gradient = w_ij(J->position) - vec3(
									w_ij(J->position + vec3(delta, 0., 0.) * (num)state.smoothing_length_h),
									w_ij(J->position + vec3(0., delta, 0.) * (num)state.smoothing_length_h),
									w_ij(J->position + vec3(0., 0., delta) * (num)state.smoothing_length_h)
								);

								//printf("Wij gradient: %f  ,  %f  ,  %f        \n", w_ij_gradient.x, w_ij_gradient.y, w_ij_gradient.z);

								sum_pressure +=
									J->mass *
									(I->pressure / (I->density*I->density) + J->pressure / (J->density*J->density))
									* w_ij_gradient;
								
								sum_viscosity += J->mass / J->density * (I->velocity - J->velocity)
									* (I->position - J->position) * w_ij_gradient /
									( (I->position - J->position)*(I->position - J->position) + 0.01 * (state.smoothing_length_h * state.smoothing_length_h));
							}
						}
					}
				}
			}

			//printf("Sum_pressure: %.1f , %.1f , %.1f           ", sum_pressure.x, sum_pressure.y, sum_pressure.z);
			//printf("Sum_viscosity: %.1f , %.1f , %.1f           ", sum_viscosity.x, sum_viscosity.y, sum_viscosity.z);
			vec3 pressure_force = -(I->pressure / I->density) * (num)I->density * sum_pressure;
			vec3 viscosity_force = I->mass * 1.e-6 * 2. * sum_viscosity;

			//printf("Pressure force: %.1f , %.1f , %.1f         \n", pressure_force.x, pressure_force.y, pressure_force.z);
			//printf("Viscosity force: %.1f , %.1f , %.1f         \n", viscosity_force.x, viscosity_force.y, viscosity_force.z);
			
			I->acceleration += (pressure_force + viscosity_force) / I->mass;
		}

		void particles_interactions(particle_t* A, particle_t* old_A)
		{
			if(state.p_collisions)
				particles_collision(A, old_A);

			if(state.p_gravity)
				A->acceleration += particles_gravity(old_A, octree_state.root);
		}

		Simulation(int n_particles)
		{
			// Allocate needed space for particles
			particles 			= (particle_t*)malloc(sizeof(particle_t)*n_particles);
			state.particles 	= (particle_t*)malloc(sizeof(particle_t)*n_particles);
			state.buff_particles= (particle_t*)malloc(sizeof(particle_t)*n_particles);
			state.n_particles 	= n_particles;

			octree_state.root = new tree_node_t;
			octree_root = new tree_node_t;

			ts = new std::thread[std::thread::hardware_concurrency()];

			update_settings();
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
					auto sleep_time = 1.e9 * std::chrono::nanoseconds(1) / ((num)settings.hertz);

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

		void build_octree(tree_node_t* node, unsigned int current_depth=0)
		{
			unsigned int allocated_p = 0;

			if(node->type == tree_node_t::ROOT)
			{
				octree_state.root->contained_p = (particle_t**)malloc(sizeof(particle_t*)*state.n_particles);

				// Fill octree root with pointers to particles
				for(int p_i = 0; p_i < state.n_particles; p_i++) octree_state.root->contained_p[p_i] = &state.particles[p_i];

				node->size = state.domain.size;
				node->volume = state.domain.size.x * state.domain.size.y * state.domain.size.z;
				node->n_p = state.n_particles;
			}

			node->children = new tree_node_t[8];

			vec3 size = node->size / 2.;

			for(unsigned int x = 0; x < 2; x++)
			{
				for(unsigned int y = 0; y < 2; y++)
				{
					for(unsigned int z = 0; z < 2; z++)
					{
						node->children[x*4+y*2+z] = tree_node_t
						{
							.mass = 0.,
							.n_p = 0,
							.contained_p = NULL,
							.coords = node->coords + vec3(x,y,z)*size,
							.size = size,
							.volume = size.x*size.y*size.z,
							.parent = node,
							.children = NULL,
						};

						tree_node_t* curr_node = &node->children[x*4+y*2+z];

						// Alloc for worst case scenario, all 
						// particles from parent not yet
						// allocated are in the current node
						curr_node->contained_p = (particle_t**)malloc(
								(node->n_p - allocated_p) * sizeof(particle_t)
						);

						// Still need to iterate on the parent's list
						for(int p_i = 0; p_i < node->n_p; p_i++)
						{
							particle_t* p = node->contained_p[p_i];

							// If particle is already allocated
							// in other child, continue
							if(p == NULL) continue;

							bool is_in_node =
								p->position.x >= curr_node->coords.x &&
								p->position.y >= curr_node->coords.y &&
								p->position.x < curr_node->coords.x + curr_node->size.x &&
								p->position.y < curr_node->coords.y + curr_node->size.y;

							if(is_in_node)
							{
								// Add to child
								curr_node->contained_p[curr_node->n_p] = p;
								curr_node->mass += p->mass;

								// Weighted center of mass (literally)
								curr_node->center_of_mass +=
									p->position*p->mass / curr_node->mass;

								// Remove from parent
								node->contained_p[p_i] = NULL;

								curr_node->n_p++;
							}
						}
						
						if(curr_node->n_p == 0)
						{
							free(curr_node->contained_p);
							curr_node->type = tree_node_t::EMPTY;
							continue;
						}

						if(curr_node->n_p <= octree_state.min_p_in_node || current_depth >= octree_state.max_depth)
						{
							curr_node->type = tree_node_t::LEAF;

							// Resize to fit the actual number of particles in node
							curr_node->contained_p = (particle_t**)realloc(
									curr_node->contained_p,
									curr_node->n_p * sizeof(particle_t)
							);
						}
						else
						{
							curr_node->type = tree_node_t::BRANCH;
							build_octree(curr_node, current_depth+1);
						}


						allocated_p += curr_node->n_p;
					} // z
				} // y
			} // x

			// Only leaf nodes point to particles
			free(node->contained_p);
		}

		void destroy_octree(tree_node_t* node)
		{
			for(int n_i = 0; n_i < 8; n_i++)
			{
				tree_node_t* curr_node = &node->children[n_i];

				switch(curr_node->type)
				{
				// Has kids
				case tree_node_t::ROOT:
				case tree_node_t::BRANCH:
					destroy_octree(curr_node);
					delete curr_node->children;
					break;

				// Has particles
				case tree_node_t::LEAF:
					free(curr_node->contained_p);
					break;

				// Has none of those
				case tree_node_t::EMPTY:
					break;
				}
			}
		}

		void build_bounding_boxes()
		{
			// Free particles lists before creating new one.
			// Since malloc wasn't called for lists of length
			// 0, don't call free for those
			for(int i = 0; i < bboxes.nbx*bboxes.nby*bboxes.nbz; i++)
				if(bboxes.n_p_per_b[i] > 0) free(bboxes.p_per_b[i]);

			bool has_changed =
				bboxes.nbx != settings.n_bounding_boxes_x ||
				bboxes.nby != settings.n_bounding_boxes_y ||
				bboxes.nbz != settings.n_bounding_boxes_z;

			// TODO: Fix race condition where user changes n_bounding_boxes
			// with imgui between update_settings() railguards and now
			// causing settings.n_bounding_boxes to possibly be < 1
			bboxes.nbx = settings.n_bounding_boxes_x;
			bboxes.nby = settings.n_bounding_boxes_y;
			bboxes.nbz = settings.n_bounding_boxes_z;

			// Realloc space, bounding boxes sizes have changed
			if(has_changed)
			{
				bboxes.n_p_per_b = (int*)realloc(bboxes.n_p_per_b,
						sizeof(int)*bboxes.nbx*bboxes.nby*bboxes.nbz);

				bboxes.p_per_b   = (particle_t***)realloc(bboxes.p_per_b,
						sizeof(particle_t**)*bboxes.nbx*bboxes.nby*bboxes.nbz);
			}

			// Incremental list, counting current index of particle when inserting
			// (see last loop inserting pointers)
			unsigned int n_p_per_b_again[bboxes.nbx*bboxes.nby*bboxes.nbz];

			// Assigning all 0s to 3D array containing the number
			// of particles per bounding box
			memset(bboxes.n_p_per_b, 0, sizeof(unsigned int)*bboxes.nbx*bboxes.nby*bboxes.nbz);
			memset(n_p_per_b_again , 0, sizeof(unsigned int)*bboxes.nbx*bboxes.nby*bboxes.nbz);

			// Counting number of particles per bounding box
			for(int p_i = 0; p_i < state.n_particles; p_i++)
			{
				particle_t* p = &state.particles[p_i];

				int x = (int)floor(p->position.x/state.domain.size.x * bboxes.nbx);
				int y = (int)floor(p->position.y/state.domain.size.y * bboxes.nby);
				int z = (int)floor(p->position.z/state.domain.size.z * bboxes.nbz);
				int xyz = x*bboxes.nby*bboxes.nbz + y*bboxes.nbz + z;

				bool is_out =   x < 0 || y < 0 || z < 0 ||
								x > bboxes.nbx-1 || y > bboxes.nby-1 || z > bboxes.nbz-1;
				
				// Some particles might be leaking out of the domain
				if(is_out)
				{
					p->bbox_xyz = ivec3(-1);
					continue;
				}
				p->bbox_xyz = glm::vec<3, int>(x, y, z);

				bboxes.n_p_per_b[xyz]++;
			}

			// Alloc space needed in list of particles per bounding box.
			// Only call malloc if lists contain particles
			for(int i = 0; i < bboxes.nbx*bboxes.nby*bboxes.nbz; i++)
				if(bboxes.n_p_per_b[i] > 0) bboxes.p_per_b[i] = (particle_t**)malloc(sizeof(particle_t*)*bboxes.n_p_per_b[i]);

			// Insert particles pointers in bounding boxes
			for(int p_i = 0; p_i < state.n_particles; p_i++)
			{
				particle_t* p = &state.particles[p_i];

				if(p->bbox_xyz == ivec3(-1)) continue;

				int xyz = 	p->bbox_xyz.x*bboxes.nby*bboxes.nbz +
							p->bbox_xyz.y*bboxes.nbz 			+
							p->bbox_xyz.z;

				// Count current particle index in bounding box
				bboxes.p_per_b[xyz][n_p_per_b_again[xyz]++] = p;
			}
		}

		void tick()
		{
			auto now = std::chrono::high_resolution_clock::now();
			state.delta_t = (now - state.last_update).count() / 1.e9; 	// count gives nanoseconds, *1e9 to get seconds
			num delta_t = state.delta_t * settings.speed; 			// Scale this tick's delta_t depending on sim speed
			state.last_update = now;

			update_settings();
			if(state.p_collisions || state.liquid) build_bounding_boxes();

			octree_state.root->type = tree_node_t::ROOT;
			if(state.p_gravity) build_octree(octree_state.root);

			// Ceil to get every particles in case it's not round
			int p_per_thread = ceil((num)state.n_particles/(num)state.n_threads);

			// For every thread
			// Gotta pass t_i by value, else it'll continue changing while the thread is running
			for(int t_i = 0; t_i < state.n_threads; t_i++) ts[t_i] = std::thread([&, t_i]()
			{
				// Iterate on a given list of particles
				for(int p_i = t_i*p_per_thread; p_i < (t_i+1)*p_per_thread && p_i < state.n_particles; p_i++)
				{
					particle_t* n_p = &state.buff_particles[p_i];
					particle_t* p = &state.particles[p_i];

					*n_p = *p;
					n_p->acceleration = vec3(0., 0., 0.); // Reset acceleration each frame

					particles_interactions(n_p, p);
					domain_interactions(n_p, p);

					// If liquid sim, compute density and
					// delay integration to second loop
					if(state.liquid) compute_liquid_density(n_p, p);

					// Else integrate
					else
					{
						n_p->velocity += n_p->acceleration * (num)delta_t;		//  v(t) = v(t-1) + a(t)*t
						n_p->position += n_p->velocity * (num)delta_t; 			//  p(t) = p(t-1) + v(t)*t
						particles[p_i] = *p;
					}
				}

				// If liquid sim, second loop to compute
				// pressure and viscosity
				if(state.liquid) for(int p_i = t_i*p_per_thread; p_i < (t_i+1)*p_per_thread && p_i < state.n_particles; p_i++)
				{
					particle_t* n_p = &state.buff_particles[p_i];
					particle_t* p = &state.particles[p_i];

					if(n_p->density > 0.)
						compute_liquid_pressure_viscosity(n_p, p);

					n_p->velocity += n_p->acceleration * (num)delta_t;		//  v(t) = v(t-1) + a(t)*t
					n_p->position += n_p->velocity * (num)delta_t; 			//  p(t) = p(t-1) + v(t)*t

					particles[p_i] = *p;
				}

			});

			// End threads
			for(int t_i= 0; t_i < state.n_threads; t_i++)
				ts[t_i].join();

			tree_node_t* last_frame_root = octree_root;
			octree_root = octree_state.root;
			octree_state.root = octree_root;
			if(state.p_gravity) destroy_octree(octree_state.root);

			// Swap both particle buffers
			particle_t* last_frame_data = state.particles;
			state.particles = state.buff_particles;
			state.buff_particles = last_frame_data;
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

			delete octree_state.root;
			//delete octree_root;

			free(particles);
			free(state.particles);
			free(state.buff_particles);
		}
};

#endif
