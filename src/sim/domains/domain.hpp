#ifndef DOMAIN_H
#define DOMAIN_H

#include "sim/particle.hpp"
#include "util/maths.hpp"
#include "util/geometry/shapes/shape.hpp"

class Domain
{
	public:

		/**
		* Applies the forces the domain is exercing on a particle
		* @param new_p Pointer to the new version of a given particle on which the forces will be applied
		* @param old_p Read only pointer to the t-1 version of a given particle
		*/
		virtual void interactions(Particle* new_p, const Particle* old_p) const = 0;

		/*-----------------*\
		*					*
		*		GETTERS		*
		*					*
		\*-----------------*/

		/**
		* Get a const reference to the domain's shape
		* @return Reference to the domain's shape
		*/
		virtual const Shape& get_shape() const = 0;

		/**
		* Compute the total volume occupied by the domain
		* @return The total volume
		*/
		virtual num get_volume() const = 0;

		/**
		* Get the domain's bounciness
		* @return the domain's bounciness
		*/
		virtual num get_bounciness() const = 0;


		/*-----------------*\
		*					*
		*		SETTERS		*
		*					*
		\*-----------------*/

		/**
		* Set the domain's size using a shape
		* @param shape The new domain's dimensions
		*/
		virtual void set_shape(const Shape& shape) = 0;

		// Virtual destructor
		virtual ~Domain() = default;
};

#endif
