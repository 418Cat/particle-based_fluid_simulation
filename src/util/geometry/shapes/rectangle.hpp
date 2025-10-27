#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "util/maths.hpp"
#include "util/geometry/shapes/shape.hpp"

class Rectangle : public Shape
{
	public:
		vec3 size;

		~Rectangle() override {}

		Rectangle(vec3 size)
		{
			this->size = size;
		}

		num volume() const override
		{
			return
				this->size.x *
				this->size.y *
				this->size.z;
		}

		shape_t shape_type() const override
		{
			return shape_t::RECTANGLE;
		}

		vec3 fitting_rectangle() const override
		{
			return this->size;
		}
};

#endif
