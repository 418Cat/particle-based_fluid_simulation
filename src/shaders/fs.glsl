#version 460 core

in vec2 fs_uv;
in vec2 velocity;

uniform vec2 domain_size;
uniform vec2 window_size;
uniform float particle_radius;
uniform float arrow_max_vel;

out vec4 col;

// Inigo Quilez
// https://iquilezles.org/articles/distfunctions2d/
float sdSegment( in vec2 p, in vec2 a, in vec2 b )
{
    vec2 pa = p-a, ba = b-a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
}

float velocity_arrow_dist(vec2 uv, vec2 velocity)
{
	float head_width = 0.4;
	vec2 head_left = (1.-head_width)*velocity + head_width*vec2(-velocity.y, velocity.x);
	vec2 head_right = (1.-head_width)*velocity - head_width*vec2(-velocity.y, velocity.x);

	return
		min(sdSegment(uv, vec2(0.), velocity/arrow_max_vel),
			min(
				sdSegment(uv, velocity/arrow_max_vel, head_left/arrow_max_vel),
				sdSegment(uv, velocity/arrow_max_vel, head_right/arrow_max_vel)
			)
	);
}

void main()
{
	float radius = length(fs_uv);

	if(radius > particle_radius) discard;

	if(radius > .95) col = vec4(1., .3, 0., 1.);
	else col = vec4(1.);

	float vel = length(velocity);

	vec2 vel_dir = velocity/vel;
	vec2 uv_dir = fs_uv/radius;

	// Color the arrow with gradient based on velocity
	if(velocity_arrow_dist(fs_uv, velocity) < vel/arrow_max_vel * 0.1)
		col = vec4(
				vec3(0., 0., 1.)*(1.-vel/arrow_max_vel) +
				vec3(1., 0., 0.)*vel/arrow_max_vel, 
		1.);
}
