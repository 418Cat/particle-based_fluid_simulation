#include "render.hpp"
#include "simulation.hpp"
#include "sim_ui.hpp"
#include "camera.hpp"

int main()
{
	SimUI ui = SimUI();

	Simulation sim = Simulation(100);
	Camera cam = Camera();
	Render render = Render(ui.window, &sim, &cam);

	ui.sim = &sim;
	ui.render = &render;
	ui.camera = &cam;

	// Settings =============================
	sim.settings.hertz = 1000;
	sim.settings.n_threads = 8;
	sim.settings.n_bounding_boxes_x = 20;
	sim.settings.n_bounding_boxes_y = 20;
	sim.settings.particle_radius = 1.;
	sim.settings.domain_gravity_radial = false;
	sim.settings.domain_gravity = vec3(0., 0., 0.);
	sim.settings.domain_size = vec3(100., 100., 100.);
	sim.settings.domain_bounciness = 0.4;
	sim.settings.particles_bounciness = 0.95;
	sim.settings.particle_gravity = true;
	sim.settings.speed = 0.5;

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
		cam.treat_inputs(ui.window, 0.013);
		render.frame();
		ui.render_stop();

		n_frame++;
	}

	return ui.end();
}
