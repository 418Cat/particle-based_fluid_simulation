#include "geometric.hpp"
#include "imgui.h"

#include "render.hpp"
#include "simulation.hpp"
#include "fluid_ui.hpp"
#include <chrono>

int main()
{
	FlUId::begin();

	Simulation sim = Simulation(
		2000,
		glm::vec2(500., 500.)		
	);

	Render render = Render(FlUId::window, &sim);

	unsigned int max_threads = std::thread::hardware_concurrency();
	if(max_threads < 1) max_threads = 1;

	// UI Settings & states =========
	bool compute_momentum = false;
	float last_momentum = 0.;
	float captured_momentum = 0.; // User selected momentum at a certain frame

	bool always_autoresize = true;

	bool show_all_particles = false;
	// ==============================

	// Temp sim settings ============
	render.zoom = 2.;
	//sim.domain.radial_gravity = true;
	sim.speed = 1.;
	sim.domain.bounciness = .9;
	sim.particles_bounciness = .9;
	// ===============================

	sim.run();

	while(FlUId::render_start())
	{
		ImGui::BeginTabBar("Debug tab bar");

		if(ImGui::BeginTabItem("Settings"))
		{
			ImGui::SeparatorText("Simulation");
			if(ImGui::Button("Pause") && sim.should_run) sim.should_run = false;
			ImGui::SameLine();
			if(ImGui::Button("Resume") && !sim.should_run)
			{
				sim.domain.last_update = std::chrono::high_resolution_clock::now();
				sim.run();
			}

			ImGui::SliderFloat("Sim speed", &sim.speed, 0., 2.);
			ImGui::SliderInt("Sim rate", (int*)&sim.sim_hertz, 1, 10000, "%d Hz");
			ImGui::SliderInt("Number of threads", (int*)&sim.n_threads, 1, max_threads);

			ImGui::Spacing();
			ImGui::Checkbox("Radial gravity", &sim.domain.radial_gravity);
			ImGui::SliderFloat2("Domain gravity, m/s", &sim.domain.gravity.x, -20., 20.);

			ImGui::Spacing();
			ImGui::SliderFloat("Domain bounciness", &sim.domain.bounciness, 0., 1.5);
			ImGui::SliderFloat("Particles bounciness", &sim.particles_bounciness, 0., 1.5);
			ImGui::SliderFloat("Particle radius", &sim.p_radius, 0.01, 20.);

			ImGui::Spacing();
			ImGui::DragFloat2("Domain size", (float*)&sim.domain.size);
			if(ImGui::Button("Reset particles"))
				sim.spawn_particles_as_rect();

			ImGui::SeparatorText("Display");
			ImGui::DragFloat("zoom", &render.zoom, render.zoom/10., 0., 100.);

			if(ImGui::Button("Auto resize") || always_autoresize)
			{
				render.zoom = std::min(
						render.win_x / sim.domain.size.x,
						render.win_y / sim.domain.size.y 
				);
			}
			ImGui::SameLine();
			ImGui::Checkbox("Always", &always_autoresize);

			ImGui::DragFloat("Arrow velocity", &render.arrow_max_vel, render.arrow_max_vel/10., 0.);

			ImGui::EndTabItem();
		}

		if(ImGui::BeginTabItem("Info"))
		{
			ImGui::Text("Number of particles: %d", sim.n_particles);

			ImGui::Separator();
			ImGui::Text("Sim time: %.1f µs", sim.last_delta_t*1.e6);
			ImGui::Text("Sim goal time: %.1f µs", 1.e6 / (float)sim.sim_hertz);
			ImGui::Text("Performance:");
			ImGui::ProgressBar(1./ (sim.last_delta_t * (float)sim.sim_hertz));

			ImGui::Separator();

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

				ImGui::Text("Momentum: %.1f kg*m/s\nRelative deviation: %.1f%%\nRelative to saved: %.1f%%",
						total_mom,
						(last_momentum - total_mom)/total_mom * 100.,
						(last_momentum - captured_momentum)/captured_momentum * 100.
				);

				last_momentum = total_mom;

				if(ImGui::Button("Save momentum"))
					captured_momentum = total_mom;

				ImGui::Separator();
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

			ImGui::EndTabItem();
		}

		ImGui::EndTabBar();

		render.frame();
		FlUId::render_stop();
	}

	std::cout << "Exiting" << std::endl;

	return FlUId::end();
}
