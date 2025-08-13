#ifndef CAMERA_H
#define CAMERA_H

#include "geometric.hpp"
#include "ext/matrix_transform.hpp"

#include <imgui.h>
#include <GLFW/glfw3.h>

using glm::vec2,
	  glm::vec3,
	  glm::mat2,
	  glm::mat4;

class Camera
{
	private:
		double mouse_x = 0;
		double mouse_y = 0;
		bool is_dragged = false;

	public:
		float camera_speed = 20.;
		vec3 view_dir;
		vec3 pos;
		float fov = 70.;

		Camera(vec3 pos = vec3(-1.), vec3 view_dir = vec3(0., 0., 1.))
		{
			this->pos 		= pos;
			this->view_dir 	= glm::normalize(view_dir);
		}

		void turn(float x_rad, float y_rad)
		{
			mat2 rot_mat_y_xz = mat2
			(
				cos(y_rad), sin(y_rad),
				-sin(y_rad), cos(y_rad)
			);

			mat2 rot_mat_xz = mat2
			(
				cos(x_rad), sin(x_rad),
				-sin(x_rad), cos(x_rad)
			);

			vec2 y_xz = vec2(view_dir.y, view_dir.z);
			vec2 xz = vec2(view_dir.x, view_dir.z);

			y_xz = rot_mat_y_xz * y_xz;
			
			xz = rot_mat_xz * xz;

			view_dir.x = xz.x;
			view_dir.y = y_xz.x;
			view_dir.z = xz.y;

			if(view_dir.y >  0.9) view_dir.y = 0.9;
			if(view_dir.y < -0.9) view_dir.y = -0.9;

			view_dir = glm::normalize(view_dir);
		}

		mat4 view_mat()
		{
			return glm::lookAt(
					pos,
					pos + view_dir,
					vec3(0., 1., 0.)
			);
		}


		void treat_inputs(GLFWwindow* win, float delta_t)
		{
			float speed = camera_speed * delta_t *
				glfwGetKey(win, GLFW_KEY_LEFT_SHIFT) ? 2. : 1.;


			if(glfwGetKey(win, GLFW_KEY_W))
				pos += view_dir * speed;

			if(glfwGetKey(win, GLFW_KEY_S))
				pos -= view_dir * speed;

			if(glfwGetKey(win, GLFW_KEY_A))
				pos -= glm::cross(view_dir, vec3(0., 1., 0.)) * speed;

			if(glfwGetKey(win, GLFW_KEY_D))
				pos += glm::cross(view_dir, vec3(0., 1., 0.)) * speed;

			if(glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_LEFT) && !ImGui::GetIO().WantCaptureMouse)
			{
				double old_x = mouse_x;
				double old_y = mouse_y;
				glfwGetCursorPos(win, &mouse_x, &mouse_y);

				if(is_dragged)
					turn((mouse_x-old_x)/100., (mouse_y-old_y)/100.);

				is_dragged = true;
			}
			else is_dragged = false;
		}
};

#endif
