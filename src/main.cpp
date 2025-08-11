#include "geometric.hpp"
#include "imgui.h"
#include "render.hpp"
#include "simulation.hpp"
#include "sim_ui.hpp"
#include "camera.hpp"

double win_x = 0;
double win_y = 0;

void win_inputs_for_camera(Camera* cam, GLFWwindow* win, float delta_t)
{
	float speed = cam->camera_speed * delta_t *
		glfwGetKey(win, GLFW_KEY_LEFT_SHIFT) ? 2. : 1.;


	if(glfwGetKey(win, GLFW_KEY_W))
		cam->pos += cam->view_dir * speed;

	if(glfwGetKey(win, GLFW_KEY_S))
		cam->pos -= cam->view_dir * speed;

	if(glfwGetKey(win, GLFW_KEY_A))
		cam->pos -= glm::normalize(glm::cross(
					cam->view_dir, vec3(0., 1., 0.) *
					speed
		));

	if(glfwGetKey(win, GLFW_KEY_D))
		cam->pos += glm::normalize(glm::cross(
					cam->view_dir, vec3(0., 1., 0.) *
					speed
		));

	if(glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_LEFT) && !ImGui::GetIO().WantCaptureMouse)
	{
		double old_x = win_x;
		double old_y = win_y;
		glfwGetCursorPos(win, &win_x, &win_y);

		//cam->turn((win_x-old_x)/10., 0.);
		cam->turn(0., (win_y-old_y)/500.);
		//cam->turn(0., 0.1);

	}

	ImGui::DragFloat3("View dir", (float*)&cam->view_dir, 0.01, -1., 1.);
}

int main()
{
	SimUI ui = SimUI();

	Simulation sim = Simulation(300);
	Camera cam = Camera();
	Render render = Render(ui.window, &sim, &cam);

	ui.sim = &sim;
	ui.render = &render;

	// Settings =============================
	sim.settings.hertz = 1000;
	sim.settings.n_threads = 8;
	sim.settings.n_bounding_boxes_x = 1;
	sim.settings.n_bounding_boxes_y = 1;
	sim.settings.particle_radius = 1.;
	sim.settings.domain_gravity_radial = false;
	sim.settings.domain_gravity = vec3(0., -5., 0.);
	sim.settings.domain_size = vec3(100., 100., 100.);
	sim.settings.domain_bounciness = 0.4;
	sim.settings.particles_bounciness = 0.6;

	render.show_boxes = false;
	render.show_borders = false;
	render.show_vel = true;
	// ======================================
	
	const unsigned int frames_benchmark = 1000;
	unsigned int n_frame = 0;
	const bool benchmark = false;

	sim.run();

	while(ui.render_start() && (n_frame < frames_benchmark || !benchmark))
	{
		win_inputs_for_camera(&cam, ui.window, 0.013);
		render.frame();
		ui.render_stop();
		n_frame++;
	}

	return ui.end();
}
