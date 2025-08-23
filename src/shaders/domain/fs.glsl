#version 460 core

uniform vec2 window_size;

uniform vec3 domain_size;

uniform bool show_borders;
uniform float border_size;

in vec3 vs_uv;

out vec4 col;

void main()
{
	if( show_borders && (
		abs(vs_uv.x) > 1.-border_size && abs(vs_uv.y) > 1.-border_size ||
		abs(vs_uv.x) > 1.-border_size && abs(vs_uv.z) > 1.-border_size ||
		abs(vs_uv.y) > 1.-border_size && abs(vs_uv.z) > 1.-border_size
	)) col = vec4(1., 0., 0., 1.);
	else discard;
}
