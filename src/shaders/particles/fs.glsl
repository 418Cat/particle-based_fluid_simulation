#version 460 core

in vec2 fs_uv;
in vec3 velocity;
in vec3 acceleration;
in float density;
in vec3 bboxes_xyz;

uniform bool show_vel;
uniform float arrow_max_vel;

uniform bool show_accel;
uniform float arrow_max_accel;

uniform bool show_density;

uniform bool show_bboxes;
uniform vec3 bboxes;

out vec4 col;

// Inigo Quilez
// https://iquilezles.org/articles/distfunctions2d/
float sdSegment( in vec2 p, in vec2 a, in vec2 b )
{
    vec2 pa = p-a, ba = b-a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
}

float arrow_dist(vec2 uv, vec2 value, float max_val)
{
	float head_width = 0.4;
	vec2 head_left = (1.-head_width)*value + head_width*vec2(-value.y, value.x);
	vec2 head_right = (1.-head_width)*value - head_width*vec2(-value.y, value.x);

	return
		min(sdSegment(uv, vec2(0.), value/max_val),
			min(
				sdSegment(uv, value/max_val, head_left/max_val),
				sdSegment(uv, value/max_val, head_right/max_val)
			)
	);
}

void main()
{
	float radius = length(fs_uv);

	if(radius > 1.) discard;

	vec3 inside_col = vec3(1.);
	vec3 outside_col = vec3(0.3);

	if(show_vel)
	{
		vec3 vel_ratio = abs(velocity) / arrow_max_vel;

		inside_col.r = vel_ratio.x;
		inside_col.g = vel_ratio.y;
		inside_col.b = vel_ratio.z;
	}
	if(show_accel)
	{
		vec3 acc_ratio = abs(acceleration) / arrow_max_accel;

		outside_col.r = 1. - 1./acc_ratio.x;
		outside_col.g = 1. - 1./acc_ratio.y;
		outside_col.b = 1. - 1./acc_ratio.z;
	}

	col = vec4(
			inside_col  * (1.-radius) +
			outside_col * radius,
			1.);

	if(show_density)
		col.r += log(density+1.);
	
	if(show_bboxes)
		col.rgb = pow(bboxes_xyz / 30., vec3(1. / 2.4));
}
