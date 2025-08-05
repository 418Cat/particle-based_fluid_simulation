#include "render.hpp"
#include "simulation.hpp"
#include "sim_ui.hpp"

int main()
{
	SimUI ui = SimUI();

	Simulation sim = Simulation(3000);
	Render render = Render(ui.window, &sim);

	ui.sim = &sim;
	ui.render = &render;

	// Settings =============================
	sim.settings.hertz = 1300;
	sim.settings.n_threads = 8;
	sim.settings.n_bounding_boxes_x = 200;
	sim.settings.n_bounding_boxes_y = 200;
	sim.settings.particle_radius = 1.;
	sim.settings.domain_gravity_radial = true;
	sim.settings.domain_gravity = vec2(20., 20.);
	sim.settings.domain_size = vec2(400., 400.);

	render.show_boxes = false;
	render.show_vel = true;
	// ======================================
	
	const unsigned int frames_benchmark = 1000;
	unsigned int n_frame = 0;
	const bool benchmark = false;

	sim.run();

	while(ui.render_start() && (n_frame < frames_benchmark || !benchmark))
	{
		render.frame();
		ui.render_stop();
		n_frame++;
	}

	return ui.end();
}
