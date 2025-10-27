#ifndef SHAPE_H
#define SHAPE_H

#include "util/maths.hpp"

/**
 * Possible 3d shapes having non-0 volume
 */
enum shape_t
{
	SPHERE,
	RECTANGLE,
};

class Shape
{
	public:
		virtual ~Shape() = default;

		/**
		* Get the total volume contained in the shape
		* @return The shape's volume
		*/
		virtual num volume() const = 0;

		/**
		* Get the shape's type
		* @return the shape's type
		*/
		virtual shape_t shape_type() const = 0;

		/**
		* Get the size of the smallest rectangle containing the domain
		* @return The size of the rectangle
		*/
		virtual vec3 fitting_rectangle() const = 0;
};

#endif
