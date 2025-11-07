#ifndef TREE_NODE_H
#define TREE_NODE_H

#include "util/maths.hpp"
#include "sim/particle.hpp"

class TreeNode
{
	public:
		num mass = 0.;
		unsigned int n_p = 0;
		Particle** contained_p = NULL;// Must only contain particles
										// if node is LEAF or is being
										// built

		vec3 center_of_mass = vec3(0.);

		vec3 coords = vec3(0.);
		vec3 size = vec3(0.);
		num volume = 0.;

		TreeNode* parent = NULL;
		TreeNode* children = NULL; // Array of size 8

		enum node_type
		{
			ROOT,	// Doesn't have parents but has children
			BRANCH, // Has both 
			LEAF,	// Doesn't have children
			EMPTY, 	// Has no particles (So no children)
		} type;

		static TreeNode* NEW_ROOT(vec3 size)
		{ return new TreeNode(ROOT, size); }

		static TreeNode* NEW_EMPTY()
		{ return new TreeNode(EMPTY); }

		TreeNode() {}

		~TreeNode()
		{
			
		}
		
	private:
		TreeNode(node_type type, vec3 size)
		{
			this->type = type;
			this->size = size;
		}

		TreeNode(node_type t) {this->type = t;}
};

#endif
