#version 460 core

uniform vec2 window_size;

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
	col = vec4(0.);
	//float aspect_ratio = domain_size.x / domain_size.y;
	vec2 d_uv = vs_uv;
	
	vec2 pix = vs_uv * (domain_size+border_size); 

	// Domain borders
	if((abs(pix.x) > domain_size.x||
		abs(pix.y) > domain_size.y)
			&& show_borders
	) col = vec4(0., 1., 0., 1.);

	// Boxes lines
	if((mod(pix.x, domain_size.x/n_bounding_boxes.x) < box_line_size/zoom ||
	    mod(pix.y, domain_size.y/n_bounding_boxes.y) < box_line_size/zoom)
			&& show_boxes
	) col = vec4(1., 0., 0., 1.);

	// Make transparent background in case clear_screen is false
	if(col == vec4(0.)) discard;
}
