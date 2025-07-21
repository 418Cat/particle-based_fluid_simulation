#include "imgui.h"
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
		60,
		glm::vec2(80., 60.)		
	);

	Render render = Render(FlUId::window, &sim);

	auto now = _now;

	while(true)
	{
		if(_now - now > frametime)
		{
			FlUId::render_start();
			ImGui::SliderFloat("Sim speed", &sim.speed, 0., 5.);
			render.frame();
			FlUId::render_stop();
		}

		sim.tick();
	}

	return FlUId::end();
}
