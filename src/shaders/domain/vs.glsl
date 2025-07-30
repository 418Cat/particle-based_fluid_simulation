#version 460 core

layout (location=0) in vec2 vs_pos;

uniform vec2 domain_size;
uniform vec2 window_size;
uniform float zoom;

out vec2 vs_uv;

void main()
{
	gl_Position = vec4(
			zoom *
			(
				vs_pos / (window_size/domain_size)
			) * 2.02 - 1.
	, 0., 1.);
	vs_uv = vs_pos;
}
