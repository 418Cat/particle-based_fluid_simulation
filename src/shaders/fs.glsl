#version 460 core

in vec2 fs_uv;
in vec2 velocity;

uniform vec2 domain_size;
uniform vec2 window_size;
uniform float particle_radius;

out vec4 col;

// Inigo Quilez
// https://iquilezles.org/articles/distfunctions2d/
float sdSegment( in vec2 p, in vec2 a, in vec2 b )
{
    vec2 pa = p-a, ba = b-a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
}

vec3 hsv_to_rgb(float h, float s, float v)
{
    float c = s*v;
    float x = c * (1. - abs(mod(h / 60., 2.) - 1.));
    h = mod(h, 360.);
    
    vec3 rgb = vec3(c, x, 0.);
    
    if(0. <= h && h <= 60.   ) return rgb.rgb;
    if(60. <= h && h <= 120. ) return rgb.grb;
    if(120. <= h && h <= 180.) return rgb.brg;
    if(180. <= h && h <= 240.) return rgb.bgr;
    if(240. <= h && h <= 300.) return rgb.gbr;
    if(300. <= h && h <= 360.) return rgb.rbg;
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
	float radius = length(fs_uv);

	if(radius > particle_radius) discard;

	if(radius > .9) col = vec4(1., 0., 0., 1.);
	else col = vec4(1.);

	float vel = length(velocity);

	vec2 vel_dir = velocity/vel;
	vec2 uv_dir = fs_uv/radius;

	// What velocity a pix with uv_dist of 1 will represent
	float max_vel = 10.;

	// Color the arrow with gradient based on velocity
	if(velocity_arrow_dist(fs_uv, velocity, max_vel) < vel/max_vel * 0.1)
		col = vec4(
				vec3(0., 0., 1.)*(1.-vel/max_vel) +
				vec3(1., 0., 0.)*vel/max_vel, 
		1.);
}
