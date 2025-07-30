#version 460 core

layout (location=0) in vec2 vs_pos;
layout (location=1) in vec2 particle_pos;
layout (location=2) in vec2 particle_velocity;
layout (location=3) in vec2 particle_acceleration;

uniform vec2 window_size;
uniform float particle_radius;
uniform float zoom;

out vec2 fs_uv;
out vec2 velocity;
out vec2 acceleration;

void main()
{
	gl_Position = vec4(
		zoom*
		(
			vs_pos*particle_radius	+
			particle_pos
		) / window_size *
		2. - 1.

		, 0., 1.
	);

	fs_uv = vs_pos;
	velocity = particle_velocity;
	acceleration = particle_acceleration;
}
