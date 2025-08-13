#ifndef RENDER_H
#define RENDER_H

#include <iostream>

#include <glad/glad.h>

#include <GLFW/glfw3.h>

#include <gtc/type_ptr.hpp>

#include "ext/matrix_clip_space.hpp"
#include "simulation.hpp"
#include "shader.h"
#include "camera.hpp"

class Render
{
	private:
		unsigned int particles_vao;
		unsigned int particles_ebo;
		unsigned int particles_vbo;
		unsigned int particles_positions_vbo;

		unsigned int domain_vao;
		unsigned int domain_ebo;
		unsigned int domain_vbo;

		Shader* particles_shaders = NULL;
		Shader* domain_shaders = NULL;

		GLFWwindow* window = NULL;

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

		Simulation* simulation;
		Camera* camera;

	public:
		int win_x = 1;
		int win_y = 1;
		
		bool show_vel = false;
		float arrow_max_vel = 20.;

		bool show_accel = false;
		float arrow_max_accel = 20.;

		bool show_borders = true;
		int border_size = 1;
		bool show_boxes = true;
		int box_line_size = 1;

		Render(GLFWwindow* win, Simulation* sim, Camera* camera)
		{
			std::cout << "\n=============== Starting OpenGL init" << std::endl;
			this->simulation = sim;
			this->camera = camera;

			if(!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
				throw("Failed to load Glad GL Loader");
			else
				std::cout << "Loaded Glad correctly" << std::endl;

			domain_shaders = new Shader(
					"shaders/domain/vs.glsl",
					"shaders/domain/fs.glsl"
			);
			std::cout << "Domain shaders set" << std::endl;

			particles_shaders = new Shader(
					"shaders/particles/vs.glsl",
					"shaders/particles/fs.glsl"
			);
			std::cout << "Particles shaders set" << std::endl;

			particles_buffers();
			std::cout << "Particles GL Buffers generated" << std::endl;

			domain_buffers();
			std::cout << "Domain GL Buffers generated" << std::endl;

			glEnable(GL_DEPTH_TEST);

			this->window = win;

			std::cout << "=============== OpenGL init finished\n" << std::endl;
		}

		void frame()
		{
			glfwGetWindowSize(this->window, &win_x, &win_y);

			draw_domain();
			draw_particles();
		}

		void domain_buffers()
		{
			glGenVertexArrays(1, &domain_vao);
			glGenBuffers(1, &domain_vbo);
			glGenBuffers(1, &domain_ebo);

			glBindVertexArray(domain_vao);

			glBindBuffer(GL_ARRAY_BUFFER, domain_vbo);
			glBufferData(GL_ARRAY_BUFFER, sizeof(tri), tri, GL_STATIC_DRAW);

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, domain_ebo);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

			// First attribute is the vertex's position for the plane
			glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
			glEnableVertexAttribArray(0);
		}

		void draw_domain()
		{
			glBindVertexArray(domain_vao);

			domain_shaders->use();

			domain_shaders->setVec2("domain_size",
					simulation->settings.domain_size.x,
					simulation->settings.domain_size.y
			);

			domain_shaders->setVec2("window_size",
					win_x,
					win_y
			);

			domain_shaders->setVec2("n_bounding_boxes",
					simulation->settings.n_bounding_boxes_x,
					simulation->settings.n_bounding_boxes_y
			);

			domain_shaders->setBool("show_borders", this->show_borders);
			domain_shaders->setFloat("border_size", this->border_size);

			domain_shaders->setBool("show_boxes", this->show_boxes);
			domain_shaders->setFloat("box_line_size", this->box_line_size);

			glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
		}

		void particles_buffers()
		{
			glGenVertexArrays(1, &particles_vao);
			glGenBuffers(1, &particles_vbo);
			glGenBuffers(1, &particles_ebo);
			glGenBuffers(1, &particles_positions_vbo);

			glBindVertexArray(particles_vao);

			glBindBuffer(GL_ARRAY_BUFFER, particles_vbo);
			glBufferData(GL_ARRAY_BUFFER, sizeof(tri), tri, GL_STATIC_DRAW);

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, particles_ebo);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

			// First attribute is the vertex's position for the plane
			glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
			glEnableVertexAttribArray(0);

			// Declare position buffer
			glBindBuffer(GL_ARRAY_BUFFER, particles_positions_vbo);
			glBufferData(GL_ARRAY_BUFFER, simulation->n_particles()*sizeof(particle_t), NULL, GL_DYNAMIC_DRAW);
			
			// 2nd attribute is position
			glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(particle_t), (void*)0);
			glEnableVertexAttribArray(1);

			// 3rd attribute is particle velocity
			glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(particle_t), (void*)sizeof(vec3));
			glEnableVertexAttribArray(2);
			
			// 4th attribute is particle acceleration
			glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(particle_t), (void*)(sizeof(vec3)*2));
			glEnableVertexAttribArray(3);
		}

		void draw_particles()
		{
			glBindVertexArray(particles_vao);
			glBindBuffer(GL_ARRAY_BUFFER, particles_positions_vbo);
			glBufferData(GL_ARRAY_BUFFER, simulation->n_particles()*sizeof(particle_t), simulation->particles, GL_DYNAMIC_DRAW);

			glVertexAttribDivisor(0, 0); // First attribute (Quad mesh) doesn't change
			glVertexAttribDivisor(1, 1); // Second attribute (Particle pos) changes every 1 instance
			glVertexAttribDivisor(2, 1); // Third one (velocity) too
			glVertexAttribDivisor(3, 1); // Acceleration, yes also

			particles_shaders->use();

			particles_shaders->setFloat("particle_radius",
					this->simulation->settings.particle_radius
			);

			particles_shaders->setBool("show_vel", this->show_vel);
			particles_shaders->setFloat("arrow_max_vel", this->arrow_max_vel);

			particles_shaders->setBool("show_accel", this->show_accel);
			particles_shaders->setFloat("arrow_max_accel", this->arrow_max_accel);
			particles_shaders->setMat4("view_mat",
					(float*)glm::value_ptr(this->camera->view_mat())
			);

			glm::mat4 project_mat = glm::perspectiveFov(
					camera->fov*glm::pi<float>()/180.f, (float)win_x, (float)win_y, 0.1f, 1000.f
			);

			particles_shaders->setMat4("project_mat",
					(float*)glm::value_ptr(project_mat)
			);

			particles_shaders->setVec3("camera_pos",
					camera->pos.x, camera->pos.y, camera->pos.z
			);

			glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0, simulation->n_particles());
		}

		~Render()
		{
			delete particles_shaders;
			delete domain_shaders;
		}
};

#endif
