#version 460 core

in vec2 fs_uv;
in vec3 velocity;
in vec3 acceleration;
in float density;
in vec3 bboxes_xyz;

uniform bool show_vel;
uniform float color_max_vel;

uniform bool show_accel;
uniform float color_max_accel;

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
	float outside_inside_ratio = 0.5;

	vec3 vel_col = vec3(1.);
	vec3 accel_col = vec3(1.);

	if(show_vel)
	{
		vel_col.r = log(1. + abs(velocity.x)) / log(1.+color_max_vel);
		vel_col.g = log(1. + abs(velocity.y)) / log(1.+color_max_vel);
		vel_col.b = log(1. + abs(velocity.z)) / log(1.+color_max_vel);
	}

	if(show_accel)
	{
		accel_col.r = log(1. + abs(acceleration.x)) / log(1.+color_max_accel);
		accel_col.g = log(1. + abs(acceleration.y)) / log(1.+color_max_accel);
		accel_col.b = log(1. + abs(acceleration.z)) / log(1.+color_max_accel);
	}

	vec3 rgb_col =
		vel_col*inside_col*(1.-radius) / (1. - outside_inside_ratio) + 
		accel_col*outside_col*radius / outside_inside_ratio;

	if(show_density)
		rgb_col.r += density+1.;

	if(show_bboxes)
		rgb_col.rgb = log(bboxes_xyz + 1) / log(100.); // TODO: Pass max bboxes

	col = vec4(rgb_col,1.);

}
