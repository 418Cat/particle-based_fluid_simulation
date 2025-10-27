#ifndef DEFAULT_DOMAIN_H
#define DEFAULT_DOMAIN_H

#include <stdexcept>

#include "util/maths.hpp"
#include "util/geometry/shapes/rectangle.hpp"

#include "sim/particle.hpp"
#include "sim/domains/domain.hpp"

class RectangleDomain : public Domain
{
	public:
		bool radial_gravity = false;
		bool gravity_axis[3] = {false, true, false};
		vec3 gravity = vec3(0., -9.81, 0.);
		num bounciness = .9;
		Rectangle rect;

	public:

		RectangleDomain(vec3 size) : rect(size)
		{}

		// lkdjhqdsf
		// qdlkjqfd
		// qdlfkjqhsdf
		//

		~RectangleDomain() {}

		void interactions(Particle* new_p, const Particle* old_p) const override
		{
			if(this->radial_gravity)
			{
				vec3 to_center = this->rect.size/2. - old_p->position;

				num dist = distance(new_p->position, this->rect.size/2.);

				num gravity_norm = sqrt(
					(this->gravity_axis[0] ? glm::abs(this->gravity.x) : 0.) +
					(this->gravity_axis[1] ? glm::abs(this->gravity.y) : 0.) +
					(this->gravity_axis[2] ? glm::abs(this->gravity.z) : 0.)
				);

				if(dist == 0.) dist = 0.01;
				to_center *= gravity_norm / dist;

				new_p->acceleration += to_center / new_p->mass;
			}
			else
				new_p->acceleration += vec3(
					this->gravity_axis[0] ? this->gravity.x : 0.,
					this->gravity_axis[1] ? this->gravity.y : 0.,
					this->gravity_axis[2] ? this->gravity.z : 0.
				);


			// Sphere shaped domain ______________________________________________
			//double radius = 50.;

			//double dist_x = p->position.x - this->size.x/2.;
			//double dist_y = p->position.y - this->size.y/2.;
			//double dist_z = p->position.z - this->size.z/2.;

			//double dist_sqrd = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

			//if(dist_sqrd > radius*radius)
			//{
			//	double dist = sqrt(dist_sqrd);

			//	vec3 normal = p->position - this->size/2.;
			//	normal /= dist;

			//	p->position = this->size/2. + normal * radius;

			//	double dot_pos_norm = dot(normal, p->velocity);
			//	p->velocity -= 2.*normal*dot_pos_norm;
			//	p->velocity *= this->bounciness;
			//}
			// ___________________________________________________________________

			// For every component of the position, check if out of domain
			for(int i = 0; i < 3; i++)
			{
				num* p_pos    = &((num*)&(new_p->position))[i];
				num* p_vel    = &((num*)&(new_p->velocity))[i];
				num* p_accel  = &((num*)&(new_p->acceleration))[i];
				num* p_radius = &new_p->radius;

				num* w_coord  = &((num*)&(this->rect.size))[i];

				// Outer wall check
				if(*p_pos > *w_coord - *p_radius)
				{
					//*p_accel = 0.;
					*p_vel  *= -this->bounciness;
					*p_pos   = *w_coord - *p_radius;
				}

				if(*p_pos < *p_radius)
				{
					//*p_accel = 0.;
					*p_vel  *= -this->bounciness;
					*p_pos   = *p_radius;
				}
			}
		}


		const Rectangle& get_shape() const override
		{
			return this->rect;
		}

		num get_volume() const override
		{ return this->rect.volume(); }

		num get_bounciness() const override
		{ return this->bounciness; }



		void set_shape(const Shape& shape) override
		{
			if(shape.shape_type() != shape_t::RECTANGLE)
				throw new std::invalid_argument(
					"Err: Size of RectangleDomain:Domain instance must be set with Rectangle:Shape instance"
			);

			this->rect = (Rectangle&)(shape);
		}

};

#endif
