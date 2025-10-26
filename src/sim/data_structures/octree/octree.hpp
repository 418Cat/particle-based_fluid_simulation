#ifndef OCTREE_H
#define OCTREE_H

#include "util/maths.hpp"
#include "sim/particle.hpp"
#include "sim/data_structures/octree/tree_node.hpp"

class Octree
{
	public:
		unsigned int max_depth = 5;
		unsigned int min_p_in_node = 4;

		Particle** particles;
		const unsigned int& n_particles;

		const vec3& size;

		TreeNode* root;

		Octree(const vec3& size, Particle** particles, const unsigned int& n_particles) : size(size), n_particles(n_particles)
		{
			this->particles = particles;
			this->root = TreeNode::NEW_ROOT(size);

		}

		void build(TreeNode* node=NULL, unsigned int current_depth=0)
		{
			if(node == NULL)
				node = this->root;

			unsigned int allocated_p = 0;

			if(node->type == TreeNode::ROOT)
			{
				node->contained_p = (Particle**)malloc(sizeof(Particle*)*n_particles);

				// Fill octree root with pointers to particles
				for(int p_i = 0; p_i < n_particles; p_i++) node->contained_p[p_i] = particles[p_i];

				node->size = size;
				node->volume = size.x * size.y * size.z;
				node->n_p = n_particles;
			}

			node->children = new TreeNode[8];

			vec3 size = node->size / 2.;

			for(unsigned int x = 0; x < 2; x++)
			{
				for(unsigned int y = 0; y < 2; y++)
				{
					for(unsigned int z = 0; z < 2; z++)
					{
						TreeNode child = TreeNode();

						child.mass = 0.;
						child.n_p = 0;
						child.contained_p = NULL;
						child.coords = node->coords + vec3(x,y,z)*size;
						child.size = size;
						child.volume = size.x*size.y*size.z;
						child.parent = node;
						child.children = NULL;

						node->children[x*4 + y*2 + z] = child;

						// Alloc for worst case scenario, all 
						// particles from parent not yet
						// allocated are in the current node
						child.contained_p = (Particle**)malloc(
								(node->n_p - allocated_p) * sizeof(Particle)
						);

						// Still need to iterate on the parent's list
						for(int p_i = 0; p_i < node->n_p; p_i++)
						{
							Particle* p = node->contained_p[p_i];

							// If particle is already allocated
							// in other child, continue
							if(p == NULL) continue;

							bool is_in_node =
								p->position.x >= child.coords.x &&
								p->position.y >= child.coords.y &&
								p->position.x <  child.coords.x + child.size.x &&
								p->position.y <  child.coords.y + child.size.y;

							if(is_in_node)
							{
								// Add to child
								child.contained_p[child.n_p] = p;
								child.mass += p->mass;

								// Weighted center of mass (literally)
								child.center_of_mass +=
									p->position*p->mass / child.mass;

								// Remove from parent
								node->contained_p[p_i] = NULL;

								child.n_p++;
							}
						}
						
						if(child.n_p == 0)
						{
							free(child.contained_p);
							child.type = TreeNode::EMPTY;
							continue;
						}

						if(child.n_p <= this->min_p_in_node || current_depth >= this->max_depth)
						{
							child.type = TreeNode::LEAF;

							// Resize to fit the actual number of particles in node
							child.contained_p = (Particle**)realloc(
									child.contained_p,
									child.n_p * sizeof(Particle)
							);
						}
						else
						{
							child.type = TreeNode::BRANCH;
							build(&child, current_depth+1);
						}


						allocated_p += child.n_p;
					} // z
				} // y
			} // x

			// Only leaf nodes point to particles
			free(node->contained_p);
		}

		static int total_nodes(TreeNode* node)
		{
			int total = 0;

			for(int n_i = 0; n_i < 8; n_i++)
			{
				TreeNode* curr_node = &node->children[n_i];
				
				switch(curr_node->type)
				{
					case TreeNode::ROOT:
					case TreeNode::BRANCH:
						total += total_nodes(curr_node);
						total++;
						break;

					case TreeNode::LEAF:
						total++;
						break;
					
					case TreeNode::EMPTY:
						break;
				}
			}

			return total;
		}

		static void destroy_octree(TreeNode* node)
		{
			for(int n_i = 0; n_i < 8; n_i++)
			{
				TreeNode* curr_node = &node->children[n_i];

				switch(curr_node->type)
				{
				// Has kids
				case TreeNode::ROOT:
				case TreeNode::BRANCH:
					destroy_octree(curr_node);
					delete curr_node->children;
					break;

				// Has particles
				case TreeNode::LEAF:
					free(curr_node->contained_p);
					break;

				// Has none of those
				case TreeNode::EMPTY:
					break;
				}
			}
		}

		~Octree()
		{
			// TODO: Fix memory leak
			destroy_octree(this->root);
		}
};

#endif
