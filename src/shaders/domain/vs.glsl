#version 460 core

layout (location=0) in vec3 vs_pos;

uniform vec3 domain_size;
uniform mat4 view_mat;
uniform mat4 project_mat;

out vec3 vs_uv;

void main()
{
	gl_Position =
		project_mat * view_mat *
		vec4((vs_pos/2.+.5)*domain_size, 1.);

	vs_uv = vs_pos;
}
