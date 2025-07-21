#version 460 core

layout (location=0) in vec2 vs_pos;
layout (location=1) in vec2 particle_pos;

uniform vec2 domain_size;

out vec2 fs_uv;

void main()
{
	gl_Position = vec4(
			(vs_pos + particle_pos)/domain_size*2. - 1.
			, 0., 1.);

	fs_uv = vs_pos;
}
