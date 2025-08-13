#version 460 core

layout (location=0) in vec2 vs_pos;
layout (location=1) in vec3 particle_pos;
layout (location=2) in vec3 particle_velocity;
layout (location=3) in vec3 particle_acceleration;

uniform float particle_radius;
uniform mat4 view_mat;
uniform mat4 project_mat;
uniform vec3 camera_pos;

out vec2 fs_uv;
out vec3 velocity;
out vec3 acceleration;

void main()
{
	vec3 z = normalize(camera_pos - particle_pos);
	vec3 y = normalize(cross(z, vec3(0., 1., 0.)));
	vec3 x = normalize(cross(y, z));

	vec3 pos = vs_pos.x*x + vs_pos.y*y;
	pos *= particle_radius;
	pos += particle_pos;

	gl_Position =
		project_mat * view_mat *
		vec4(pos, 1.);

	fs_uv = vs_pos;
	velocity = particle_velocity;
	acceleration = particle_acceleration;
}
