#ifndef RENDER_H
#define RENDER_H

#include <iostream>

#include <glad/glad.h>

#include <GLFW/glfw3.h>

#include <gtc/type_ptr.hpp>

#include "ext/matrix_clip_space.hpp"
#include "sim/simulation.hpp"

#include "shader.h"

#include "math.h"

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

		unsigned int boxes_vao;
		unsigned int boxes_ebo;
		unsigned int boxes_vbo;

		Shader* particles_shaders = NULL;
		Shader* domain_shaders = NULL;

		GLFWwindow* window = NULL;

		glm::mat4 project_mat = glm::perspectiveFov(70., 10., 10., 0.1, 100.);

		float quad_vert[8] = 
		{
			-1.f, -1.f,
			 1.f, -1.f,
			 1.f,  1.f,
			-1.f,  1.f,
		};

		unsigned int quad_ind[6] =
		{
			0, 1, 3,
			1, 2, 3
		};


		float cube_vert[24]
		{
			-1, -1,  1, //0
			 1, -1,  1, //1
			-1,  1,  1, //2
			 1,  1,  1, //3
			-1, -1, -1, //4
			 1, -1, -1, //5
			-1,  1, -1, //6
			 1,  1, -1  //7
		};

		unsigned int cube_ind[36]
		{
			//Top
			2, 6, 7,
			2, 3, 7,

			//Bottom
			0, 4, 5,
			0, 1, 5,

			//Left
			0, 2, 6,
			0, 4, 6,

			//Right
			1, 3, 7,
			1, 5, 7,

			//Front
			0, 2, 3,
			0, 1, 3,

			//Back
			4, 6, 7,
			4, 5, 7
		};

		Simulation* simulation;
		Camera* camera;

	public:
		int win_x = 1;
		int win_y = 1;
		
		bool show_vel = false;
		float color_max_vel = 20.;

		bool show_accel = false;
		float color_max_accel = 20.;

		bool show_borders = true;
		float border_size = 1.;

		bool show_boxes = false;
		bool show_density = true;

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

			project_mat = glm::perspectiveFov(camera->fov*glm::pi<float>()/180.f, (float)win_x, (float)win_y, 0.1f, 1000.f);

			draw_domain();
			draw_particles();
			// boxes_buffers(); // Causes issues when called each frame, gotta rewrite this class
		}

		void boxes_buffers()
		{
			glGenVertexArrays(1, &boxes_vao);
			glGenBuffers(1, &boxes_vbo);
			glGenBuffers(1, &boxes_ebo);

			glBindVertexArray(boxes_vao);


			glBindBuffer(GL_ARRAY_BUFFER, boxes_vao);
			glBufferData(GL_ARRAY_BUFFER, 24*sizeof(float), cube_vert, GL_STATIC_DRAW);

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, boxes_ebo);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, 16*sizeof(unsigned int), cube_ind, GL_STATIC_DRAW);

			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
			glEnableVertexAttribArray(0);

			
			std::cout << "Octree has " << Octree::total_nodes(simulation->octree_root) << " nodes" << std::endl;
		}

		void domain_buffers()
		{
			glGenVertexArrays(1, &domain_vao);
			glGenBuffers(1, &domain_vbo);
			glGenBuffers(1, &domain_ebo);

			glBindVertexArray(domain_vao);

			glBindBuffer(GL_ARRAY_BUFFER, domain_vbo);
			glBufferData(GL_ARRAY_BUFFER, sizeof(cube_vert), cube_vert, GL_STATIC_DRAW);

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, domain_ebo);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cube_ind), cube_ind, GL_STATIC_DRAW);

			// First attribute is the vertex's position for the cube
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
			glEnableVertexAttribArray(0);
		}

		void draw_domain()
		{
			glBindVertexArray(domain_vao);

			domain_shaders->use();

			domain_shaders->setVec3("domain_size",
					simulation->settings.domain_size.x,
					simulation->settings.domain_size.y,
					simulation->settings.domain_size.z
			);

			domain_shaders->setVec2("window_size",
					win_x,
					win_y
			);

			domain_shaders->setBool("show_borders", this->show_borders);
			domain_shaders->setFloat("border_size", this->border_size/100.);
			domain_shaders->setMat4("view_mat", (float*)glm::value_ptr(camera->view_mat()));
			domain_shaders->setMat4("project_mat", (float*)glm::value_ptr(project_mat));

			glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
		}

		void particles_buffers()
		{
			glGenVertexArrays(1, &particles_vao);
			glGenBuffers(1, &particles_vbo);
			glGenBuffers(1, &particles_ebo);
			glGenBuffers(1, &particles_positions_vbo);

			glBindVertexArray(particles_vao);

			glBindBuffer(GL_ARRAY_BUFFER, particles_vbo);
			glBufferData(GL_ARRAY_BUFFER, sizeof(quad_vert), quad_vert, GL_STATIC_DRAW);

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, particles_ebo);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quad_ind), quad_ind, GL_STATIC_DRAW);

			// First attribute is the vertex's position for the plane
			glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
			glEnableVertexAttribArray(0);

			// Declare position buffer
			glBindBuffer(GL_ARRAY_BUFFER, particles_positions_vbo);
			glBufferData(GL_ARRAY_BUFFER, simulation->n_particles()*sizeof(Particle), NULL, GL_DYNAMIC_DRAW);
			
			// 2nd attribute is position
			glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, sizeof(Particle), (void*)0);
			glEnableVertexAttribArray(1);

			// 3rd attribute is particle velocity
			glVertexAttribPointer(2, 3, GL_DOUBLE, GL_FALSE, sizeof(Particle), (void*)sizeof(vec3));
			glEnableVertexAttribArray(2);
			
			// 4th attribute is particle acceleration
			glVertexAttribPointer(3, 3, GL_DOUBLE, GL_FALSE, sizeof(Particle), (void*)(sizeof(vec3)*2));
			glEnableVertexAttribArray(3);

			// 5th is density
			glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(Particle), (void*)(sizeof(vec3)*3 + sizeof(num)));
			glEnableVertexAttribArray(4);

			// Bounding box xyz
			glVertexAttribPointer(5, 3, GL_INT, GL_FALSE, sizeof(Particle), (void*)(sizeof(vec3)*3 + sizeof(num) + sizeof(float)));
			glEnableVertexAttribArray(5);
		}

		void draw_particles()
		{
			glBindVertexArray(particles_vao);
			glBindBuffer(GL_ARRAY_BUFFER, particles_positions_vbo);
			glBufferData(GL_ARRAY_BUFFER, simulation->n_particles()*sizeof(Particle), simulation->particles, GL_DYNAMIC_DRAW);

			glVertexAttribDivisor(0, 0); // First attribute (Quad mesh) doesn't change
			glVertexAttribDivisor(1, 1); // Second attribute (Particle pos) changes every 1 instance
			glVertexAttribDivisor(2, 1); // Third one (velocity) too
			glVertexAttribDivisor(3, 1); // Acceleration, yes also
			glVertexAttribDivisor(4, 1); // Density
			glVertexAttribDivisor(5, 1); // Bounding box xyz

			particles_shaders->use();

			particles_shaders->setFloat("particle_radius",
					this->simulation->settings.particle_radius
			);

			particles_shaders->setBool("show_vel", this->show_vel);
			particles_shaders->setFloat("color_max_vel", this->color_max_vel);

			particles_shaders->setBool("show_accel", this->show_accel);
			particles_shaders->setFloat("color_max_accel", this->color_max_accel);
			particles_shaders->setMat4("view_mat",
					(float*)glm::value_ptr(this->camera->view_mat())
			);


			particles_shaders->setMat4("project_mat",
					(float*)glm::value_ptr(project_mat)
			);

			particles_shaders->setVec3("camera_pos",
					camera->pos.x, camera->pos.y, camera->pos.z
			);

			particles_shaders->setBool("show_bboxes",
					this->show_boxes
			);

			particles_shaders->setBool("show_density",
					this->show_density
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
