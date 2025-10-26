#ifndef BOUNDING_BOXES_H
#define BOUNDING_BOXES_H

#include "sim/particle.hpp"

class BoundingBoxes
{
	public:
		unsigned int nbx = 0;
		unsigned int nby = 0;
		unsigned int nbz = 0;

		// Flattened 3D array of arrays of pointers to particles
		Particle* * * p_per_b = NULL; 
		
		// Flattened 3D array containing the
		// number of particles per bounding box
		int* n_p_per_b = NULL; 
};

#endif
