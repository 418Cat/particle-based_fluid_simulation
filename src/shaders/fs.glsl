#version 460 core

in vec2 fs_uv;

uniform vec2 domain_size;

out vec4 col;

void main()
{
	float aspect_ratio = domain_size.x / domain_size.y;

	float pyth = fs_uv.x*fs_uv.x + fs_uv.y*fs_uv.y*aspect_ratio;

	if(pyth > 1.) discard;

	if(pyth > .7) col = vec4(1., 0., 0., 1.);
	else col = vec4(1.);
}
