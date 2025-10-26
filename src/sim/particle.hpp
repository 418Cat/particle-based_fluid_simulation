#ifndef PARTICLE_H
#define PARTICLE_H

#include "util/maths.hpp"

#include "sim/collision_info.hpp"

class Particle
{
	public:
		vec3 position 		= vec3(1.);
		vec3 velocity 		= vec3(0.);
		vec3 acceleration 	= vec3(0.);
		num mass 			= 1.;
		num radius			= 1.;

		vec3 bbox_xyz		= glm::vec<3, int>(0);

		num density 		= 1.;
		num pressure 		= 1.;

		bool collides(Particle &p, collision_info_t *collision_info);
};

#endif
