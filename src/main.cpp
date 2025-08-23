#include "render.hpp"
#include "simulation.hpp"
#include "sim_ui.hpp"
#include "camera.hpp"

int main()
{
	SimUI ui = SimUI();

	Simulation sim = Simulation(200);
	Camera cam = Camera();
	Render render = Render(ui.window, &sim, &cam);

	ui.sim = &sim;
	ui.render = &render;
	ui.camera = &cam;

	// Settings =============================
	sim.settings.domain_gravity = vec3(0., -9., 0.);
	sim.settings.domain_size = vec3(100., 100., 2.);

	render.show_vel = true;

	cam.pos = sim.settings.domain_size * vec3(0.5) - vec3(20., 0., 80.);
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
