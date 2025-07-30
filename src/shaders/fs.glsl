#version 460 core

in vec2 fs_uv;
in vec2 velocity;
in vec2 acceleration;

uniform vec2 domain_size;
uniform vec2 window_size;
uniform float particle_radius;

uniform bool show_vel;
uniform float arrow_max_vel;

uniform bool show_accel;
uniform float arrow_max_accel;

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

	if(radius > .95) col = vec4(1., .3, 0., 1.);
	else col = vec4(1.);

	if(show_vel)
	{
		float vel = length(velocity);
		float ratio = vel/arrow_max_vel;

		vec2 uv_dir = fs_uv/radius;

		// Color the arrow with gradient based on velocity
		if(arrow_dist(fs_uv, velocity, arrow_max_vel) < ratio * 0.1)
			col = vec4(
					vec3(0., 0., 1.)*(1.-ratio) +
					vec3(1., 0., 0.)*ratio, 
			1.);
	}
	if(show_accel)
	{
		float accel = length(acceleration);
		float ratio = accel/arrow_max_accel;

		vec2 uv_dir = fs_uv/radius;

		// Color the arrow with gradient based on velocity
		if(arrow_dist(fs_uv, acceleration, arrow_max_accel) < accel/arrow_max_accel * 0.1)
			col = vec4(0., clamp(ratio*0.8, 0., 1.), clamp(1.-ratio*0.8, 0., 1.), 1.);
	}
}
