#version 460 core

uniform vec2 domain_size;
uniform vec2 n_bounding_boxes;
uniform float zoom;

in vec2 vs_uv;

out vec4 col;

void main()
{
	col = vec4(0., 0., 0., 1.);

	float aspect_ratio = domain_size.x / domain_size.y;

	float border_size = 0.01;

	if(		abs(vs_uv.x) > 1. - border_size ||
			abs(vs_uv.y) > 1. - border_size
	) col = vec4(0., 1., 0., 1.);

	float boxes_size = 0.001;
	if(	mod(vs_uv.x*n_bounding_boxes.x, 1.) < boxes_size*n_bounding_boxes.x ||
		mod(vs_uv.y*n_bounding_boxes.y, 1.) < boxes_size*n_bounding_boxes.y
	) col = vec4(1., 0., 0., 1.);
}
