#include "geometric.hpp"
#include "imgui.h"
#include "imgui_internal.h"
#include "render.hpp"
#include "simulation.hpp"
#include "fluid_ui.hpp"
#include <chrono>

#define _now std::chrono::system_clock::now()
#define frametime std::chrono::milliseconds(0)

int main()
{
	FlUId::begin();

	Simulation sim = Simulation(
		20,
		glm::vec2(20., 20.)		
	);

	Render render = Render(FlUId::window, &sim);

	auto now = _now;

	while(FlUId::render_start())
	{
		//if(_now - now > frametime)
		{
			//FlUId::render_start();
			ImGui::SliderFloat("Sim speed", &sim.speed, 0., 2.);

			ImGui::Checkbox("Radial gravity", &sim.domain.radial_gravity);
			ImGui::SliderFloat2("Domain gravity, m/s", &sim.domain.gravity.x, -20., 20.);

			ImGui::SliderFloat("Domain bounciness", &sim.domain.bounciness, 0., 1.5);
			ImGui::SliderFloat("Particles bounciness", &sim.particles_bounciness, 0., 1.5);

			if(ImGui::Button("Reset particles"))
				sim.spawn_particles_as_rect();

			ImGui::Text("Frametime: %.1f ms", sim.last_delta_t*1000.);
			//ImGui::Text("Kinetic energy: %.1f");
			if(sim.speed == 0. && ImGui::Button("Next tick"))
			{
				sim.tick();
			}

			render.frame();
			FlUId::render_stop();
		}

		sim.tick();
	}

	return FlUId::end();
}
