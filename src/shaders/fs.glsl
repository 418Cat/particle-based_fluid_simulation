#version 460 core

in vec2 fs_uv;
in vec2 velocity;

uniform vec2 domain_size;

out vec4 col;

// Inigo Quilez
// https://iquilezles.org/articles/distfunctions2d/
float sdSegment( in vec2 p, in vec2 a, in vec2 b )
{
    vec2 pa = p-a, ba = b-a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
}

float velocity_arrow_dist(vec2 uv, vec2 velocity, float max_vel)
{
	float head_width = 0.4;
	vec2 head_left = (1.-head_width)*velocity + head_width*vec2(-velocity.y, velocity.x);
	vec2 head_right = (1.-head_width)*velocity - head_width*vec2(-velocity.y, velocity.x);

	return
		min(sdSegment(uv, vec2(0.), velocity/max_vel),
			min(
				sdSegment(uv, velocity/max_vel, head_left/max_vel),
				sdSegment(uv, velocity/max_vel, head_right/max_vel)
			)
	);
}

void main()
{
	float aspect_ratio = domain_size.x / domain_size.y;

	float pyth = fs_uv.x*fs_uv.x + fs_uv.y*fs_uv.y*aspect_ratio;

	if(pyth > 1.) discard;

	if(pyth > .9) col = vec4(1., 0., 0., 1.);
	else col = vec4(1.);

	float vel = length(velocity);
	float uv_dist = length(fs_uv);

	vec2 vel_dir = velocity/vel;
	vec2 uv_dir = fs_uv/uv_dist;

	// What velocity a pix with uv_dist of 1 will represent
	float max_vel = 10.;

	if(velocity_arrow_dist(fs_uv, velocity, max_vel) < vel/max_vel * 0.1) col = vec4(0., 1., 0., 1.);
	//if(velocity_arrow_dist(fs_uv, velocity, max_vel) < 0.05) col = vec4(0., 1., 0., 1.);
}
