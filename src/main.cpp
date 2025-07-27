#include "geometric.hpp"
#include "imgui.h"
#include "imgui_internal.h"
#include "render.hpp"
#include "simulation.hpp"
#include "fluid_ui.hpp"
#include <chrono>
#include <ctime>

#define _now std::chrono::system_clock::now()
#define frametime std::chrono::milliseconds(0)

int main()
{
	FlUId::begin();

	Simulation sim = Simulation(
		500,
		glm::vec2(280., 160.)		
	);

	Render render = Render(FlUId::window, &sim);

	auto now = _now;

	bool compute_momentum = false;
	bool show_all_particles = false;

	// Temp sim settings ============
	render.zoom = 2.;
	//sim.domain.radial_gravity = true;
	sim.speed = 0.;
	sim.domain.gravity.y = 0.;
	sim.domain.bounciness = 1.;
	sim.particles_bounciness = 1.;
	// ===============================

	while(FlUId::render_start())
	{
		//if(_now - now > frametime)
		{
			ImGui::BeginTabBar("Debug tab bar");

			if(ImGui::BeginTabItem("Settings"))
			{
				ImGui::SliderFloat("Sim speed", &sim.speed, 0., 2.);

				ImGui::Checkbox("Radial gravity", &sim.domain.radial_gravity);
				ImGui::SliderFloat2("Domain gravity, m/s", &sim.domain.gravity.x, -20., 20.);

				ImGui::SliderFloat("Domain bounciness", &sim.domain.bounciness, 0., 1.5);
				ImGui::SliderFloat("Particles bounciness", &sim.particles_bounciness, 0., 1.5);

				ImGui::DragFloat2("Domain size", (float*)&sim.domain.size);
				if(ImGui::Button("Reset particles"))
					sim.spawn_particles_as_rect();

				ImGui::DragFloat("zoom", &render.zoom, render.zoom/10., 0., 100.);

				ImGui::EndTabItem();
			}

			if(ImGui::BeginTabItem("Info"))
			{
				ImGui::Checkbox("Compute momentum", &compute_momentum);
				if(compute_momentum)
				{
					float total_mom = 0.;

					for(int i = 0; i < sim.n_particles; i++)
					{
						particle_t* p = &sim.particles[i];
						float mom = p->mass * length(p->velocity);
						if(mom != mom)
						{
							printf("Particle %d has NaN momentum. Mass %.1f    Vel (%.1f , %.1f): %.1f     Momentum %.1f\n",
									i, p->mass, p->velocity.x, p->velocity.y, length(p->velocity), mom);
						}

						total_mom += mom;
					}

					ImGui::Text("Momentum: %.1f kg*m/s", total_mom);
				}

				ImGui::Checkbox("Show all particles", &show_all_particles);
				if(show_all_particles)
				{
					ImGui::BeginTable("Particles", 4,
							ImGuiTableFlags_Resizable
					);
					ImGui::TableSetupColumn("i");
					ImGui::TableSetupColumn("Position (m)");
					ImGui::TableSetupColumn("Velocity (m/s)");
					ImGui::TableSetupColumn("Mass (kg)");
					ImGui::TableHeadersRow();

					for(int i = 0; i < sim.n_particles; i++)
					{
						ImGui::TableNextRow();
						particle_t* A = &sim.particles[i];

						ImGui::TableSetColumnIndex(0);
						ImGui::Text("%d", i);
						
						ImGui::TableNextColumn();
						ImGui::Text("(%.1f , %.1f)", A->position.x, A->position.y);

						ImGui::TableNextColumn();
						ImGui::Text("(%.1f , %.1f) : %.1f", A->velocity.x, A->velocity.y, length(A->velocity));

						ImGui::TableNextColumn();
						ImGui::Text("%.1f", A->mass);
					}
					ImGui::EndTable();
				}

				ImGui::Text("Frametime: %.1f ms", sim.last_delta_t*1000.);

				ImGui::EndTabItem();
			}

			ImGui::EndTabBar();

			render.frame();
			FlUId::render_stop();
		}

		sim.tick();
	}

	return FlUId::end();
}
