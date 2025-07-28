#ifndef RENDER_H
#define RENDER_H


#include <iostream>

#include "glad/glad.h"
#include <GLFW/glfw3.h>

#include "simulation.hpp"
#include "shader.h"

class Render
{
	private:
		unsigned int vao_1;
		unsigned int vao_2;
		unsigned int vbo;
		unsigned int ebo;

		Shader* shaders;

		GLFWwindow* window;

		float tri[8] = 
		{
			-1.f, -1.f,
			 1.f, -1.f,
			 1.f,  1.f,
			-1.f,  1.f,
		};

		unsigned int indices[6] =
		{
			0, 1, 3,
			1, 2, 3
		};

		unsigned int positions_vbo;

		Simulation* simulation;

	public:
		float zoom = 1.;
		int win_x, win_y;
		float arrow_max_vel = 20.;

		Render(GLFWwindow* win, Simulation* sim)
		{
			this->simulation = sim;

			if(!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
				throw("Failed to load Glad GL Loader");
			else
				std::cout << "Loaded Glad correctly" << std::endl;

			this->window = win;

			glGenVertexArrays(1, &vao_1);
			glGenBuffers(1, &vbo);
			glGenBuffers(1, &ebo);
			glGenBuffers(1, &positions_vbo);

			glBindVertexArray(vao_1);

			glBindBuffer(GL_ARRAY_BUFFER, vbo);
			glBufferData(GL_ARRAY_BUFFER, sizeof(tri), tri, GL_STATIC_DRAW);

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

			// First attribute is the vertex's position for the plane
			glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
			glEnableVertexAttribArray(0);

			// Declare position buffer
			glBindBuffer(GL_ARRAY_BUFFER, positions_vbo);
			glBufferData(GL_ARRAY_BUFFER, simulation->n_particles*sizeof(particle_t), NULL, GL_DYNAMIC_DRAW);
			
			// 2nd attribute is position
			glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(particle_t), (void*)0);
			glEnableVertexAttribArray(1);

			// 3rd attribute is particle velocity
			glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(particle_t), (void*)sizeof(vec2));
			glEnableVertexAttribArray(2);

			std::cout << "GL Buffers generated" << std::endl;

			const char* vs_path = "shaders/vs.glsl";
			const char* fs_path = "shaders/fs.glsl";
			shaders = new Shader(vs_path, fs_path);
		}

		void frame()
		{
			glBindVertexArray(vao_1);
			glBindBuffer(GL_ARRAY_BUFFER, positions_vbo);
			glBufferData(GL_ARRAY_BUFFER, simulation->n_particles*sizeof(particle_t), simulation->particles, GL_DYNAMIC_DRAW);

			glVertexAttribDivisor(0, 0); // First attribute (Quad mesh) doesn't change
			glVertexAttribDivisor(1, 1); // Second attribute (Particle pos) changes every 1 instance
			glVertexAttribDivisor(2, 1); // Third one (velocity) too

			shaders->use();

			glfwGetWindowSize(this->window, &win_x, &win_y);
			shaders->setVec2(
					"window_size",
					(float)this->win_x, (float)this->win_y
			);

			shaders->setFloat("particle_radius",
					this->simulation->p_radius
			);

			shaders->setFloat("zoom",
					this->zoom
			);

			shaders->setFloat("arrow_max_vel",
					this->arrow_max_vel
			);

			glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0, simulation->n_particles);
		}

		~Render()
		{
			//free(shaders);
		}
};

#endif
