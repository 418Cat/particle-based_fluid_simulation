#version 460 core

layout (location=0) in vec2 vs_pos;
layout (location=1) in vec3 particle_pos;
layout (location=2) in vec3 particle_velocity;
layout (location=3) in vec3 particle_acceleration;

uniform vec2 window_size;
uniform float particle_radius;
uniform float zoom;
uniform mat4 view_mat;
uniform mat4 project_mat;

out vec2 fs_uv;
out vec3 velocity;
out vec3 acceleration;

void main()
{
	gl_Position =
		project_mat * view_mat *
		vec4(vec3(vs_pos*particle_radius, 0.) + particle_pos, 1./zoom);

	fs_uv = vs_pos;
	velocity = particle_velocity;
	acceleration = particle_acceleration;
}
