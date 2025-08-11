#ifndef CAMERA_H
#define CAMERA_H

#include "ext/matrix_clip_space.hpp"
#include "geometric.hpp"
#include "glm.hpp"
#include "ext/matrix_transform.hpp"

using glm::vec2,
	  glm::vec3,
	  glm::mat2,
	  glm::mat4;

class Camera
{
	public:
		float camera_speed = 20.;
		vec3 view_dir;
		vec3 pos;

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
			
			xz = rot_mat_xz * xz * y_xz.y;

			view_dir.x = xz.x;
			view_dir.y = y_xz.x;
			view_dir.z = xz.y;

			if(	view_dir.y >  0.9 ||
				view_dir.y < -0.9)
			{
				view_dir.y = 0.95;
				x_rad = 0.;
			}

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
};

#endif
