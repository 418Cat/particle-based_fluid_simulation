#include "simulation.hpp"
#include "fluid_ui.hpp"

int main()
{
	FlUId::begin();

	Simulation sim = Simulation(
		100,
		glm::vec2(50., 50.)		
	);

	while(FlUId::render())
	{
		sim.tick();
	}

	return FlUId::end();
}
