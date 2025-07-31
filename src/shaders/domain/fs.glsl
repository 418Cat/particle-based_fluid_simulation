#version 460 core

uniform vec2 domain_size;
uniform vec2 n_bounding_boxes;
uniform float zoom;

uniform bool show_borders;
uniform float border_size;

uniform bool show_boxes;
uniform float box_line_size;

in vec2 vs_uv;

out vec4 col;

void main()
{
	col = vec4(0., 0., 0., 1.);
	//float aspect_ratio = domain_size.x / domain_size.y;
	vec2 d_uv = vs_uv;

	// Domain borders
	if((abs(vs_uv.x) > 1. - border_size ||
		abs(vs_uv.y) > 1. - border_size)
			&& show_borders
	) col = vec4(0., 1., 0., 1.);

	// Boxes lines
	if((mod(d_uv.x*n_bounding_boxes.x, 1.) < box_line_size*n_bounding_boxes.x ||
		mod(d_uv.y*n_bounding_boxes.y, 1.) < box_line_size*n_bounding_boxes.y)
			&& show_boxes
	) col = vec4(1., 0., 0., 1.);
}
