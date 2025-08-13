#ifndef SIM_UI_H
#define SIM_UI_H

#include <thread>

#include "imgui.h"
#include "render.hpp"
#include "simulation.hpp"

#include "ui.h"

// Funny name with Fluid and UI :D
class SimUI
{
	private:
		unsigned int max_threads = std::thread::hardware_concurrency();
		const float side_window_size = 0.3;

		// UI Settings & states =========
		bool compute_momentum = false;
		float last_momentum = 0.;
		float captured_momentum = 0.; // Reference momentum at a user chosen frame
		bool show_all_particles = false;
		int n_particles = 0;
		bool clear_screen = true;
		// ==============================

	public:
		GLFWwindow* window;
		Simulation* sim = NULL;
		Render* render = NULL;
		Camera* camera = NULL;
	
		SimUI()
		{
			// Couldn't get max number of threads, defaulting to 1.
			// (Assuming there's at least one thread, otherwise
			// they'd have bigger issues than this)
			if(max_threads < 1) max_threads = 1; 

			UI::ui_init();
			UI::show_ui = true;

			SimUI::window = UI::window;
		}

		// Render imgui components
		void side_window()
		{
			// Make it a left bar window
			ImGui::SetNextWindowPos(ImVec2(
						render->win_x * (1. - side_window_size),
						0.
			));
			ImGui::SetNextWindowSize(ImVec2(
						render->win_x * side_window_size,
						render->win_y	
			));
			
			ImGui::Begin(
					"side_window",
					NULL,
						ImGuiWindowFlags_NoDecoration 	|
						ImGuiWindowFlags_NoResize		|
						ImGuiWindowFlags_NoMove
			);

			ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(!sim->settings.run/2., sim->settings.run/2., 0., 1.));
			if(ImGui::Button("Pause") && sim->settings.run) sim->settings.run = false;
			ImGui::PopStyleColor();
			ImGui::SameLine();

			ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(sim->settings.run/2., !sim->settings.run/2., 0., 1.));
			if(ImGui::Button("Resume") && !sim->settings.run)
				sim->run();

			ImGui::PopStyleColor();
			ImGui::SameLine();
			ImGui::SliderFloat("Sim speed", &sim->settings.speed, 0., 2.);

			ImGui::InputInt("N_particles", &n_particles);
			if(ImGui::Button("Reset sim"))
			{
				sim->reset(n_particles);
				render->particles_buffers();
			}
			ImGui::Spacing();

			ImGui::BeginTabBar("Debug tab bar");

			if(ImGui::BeginTabItem("Settings"))
			{
				ImGui::SeparatorText("Simulation");
				ImGui::SliderInt("Sim rate", (int*)&sim->settings.hertz, 1, 10000, "%d Hz");
				ImGui::SliderInt("Number of threads", (int*)&sim->settings.n_threads, 1, max_threads);

				ImGui::Spacing();
				ImGui::Text("Bounding boxes XY");
				ImGui::SliderInt("##BboxesX", (int*)(&sim->settings.n_bounding_boxes_x), 1, (int)floor(sim->settings.domain_size.x/(sim->settings.particle_radius*2.)));
				ImGui::SliderInt("##BboxesY", (int*)(&sim->settings.n_bounding_boxes_y), 1, (int)floor(sim->settings.domain_size.y/(sim->settings.particle_radius*2.)));

				ImGui::Spacing();
				ImGui::Checkbox("Radial gravity", &sim->settings.domain_gravity_radial);
				ImGui::SliderFloat3("Domain gravity, m/s", &sim->settings.domain_gravity.x, -40., 40.);
				ImGui::Checkbox("Particles gravity", &sim->settings.particle_gravity);
				ImGui::SameLine();
				ImGui::Checkbox("Inverse", &sim->settings.particles_gravity_inverse);

				ImGui::Spacing();
				ImGui::SliderFloat("Domain bounciness", &sim->settings.domain_bounciness, 0., 1.5);
				ImGui::SliderFloat("Particles bounciness", &sim->settings.particles_bounciness, 0., 1.5);
				ImGui::SliderFloat("Particle radius", &sim->settings.particle_radius, 0.01, 20.);

				ImGui::Spacing();
				ImGui::DragFloat3("Domain size", (float*)&sim->settings.domain_size, 1., 0.);

				if(ImGui::Button("Reset particles"))
				{
					bool was_running = sim->settings.run;
					sim->end();
					sim->spawn_particles_as_rect();
					if(was_running) sim->run();
				}

				ImGui::SeparatorText("Display");
				ImGui::DragFloat("FOV", &camera->fov, camera->fov/10., 1., 100.);

				ImGui::Checkbox("Show##vel", &render->show_vel);
				ImGui::SameLine();
				ImGui::DragFloat("Velocity arrow", &render->arrow_max_vel, render->arrow_max_vel/10. + 0.0001, 0.);

				ImGui::Checkbox("Show##accel", &render->show_accel);
				ImGui::SameLine();
				ImGui::DragFloat("Acceleration arrow", &render->arrow_max_accel, render->arrow_max_accel/10. + 0.0001, 0.);

				ImGui::Spacing();
				ImGui::Checkbox("Show", &render->show_borders);
				ImGui::SameLine();
				ImGui::DragInt("Domain borders thickness", &render->border_size, 1, 1, 30);
				
				ImGui::Checkbox("Show##boxes lines", &render->show_boxes);
				ImGui::SameLine();
				ImGui::DragInt("Bounding boxes lines thickness", &render->box_line_size, 1, 1, 30);

				ImGui::Spacing();
				ImGui::Checkbox("Clear screen", &clear_screen);

				ImGui::EndTabItem();
			}

			if(ImGui::BeginTabItem("Info"))
			{
				ImGui::Text("Number of particles: %d", sim->n_particles());

				ImGui::Separator();
				ImGui::Text("Sim time: %.1f µs", sim->sim_tick_time()*1.e6);
				ImGui::Text("Sim goal time: %.1f µs", 1.e6 / (float)sim->settings.hertz);

				float ratio = 1./ (sim->sim_tick_time() * (float)sim->settings.hertz);
				ImGui::Text("Performance: %.1f %%", ratio*100.);
				ImGui::ProgressBar(ratio);

				ImGui::Separator();

				ImGui::Checkbox("Compute momentum", &compute_momentum);
				if(compute_momentum)
				{
					float total_mom = 0.;

					for(int i = 0; i < sim->n_particles(); i++)
					{
						particle_t* p = &sim->particles[i];
						float mom = p->mass * length(p->velocity);
						if(mom != mom)
						{
							printf("Particle %d has NaN momentum. Mass %.1f    Vel (%.1f , %.1f): %.1f     Momentum %.1f\n",
									i, p->mass, p->velocity.x, p->velocity.y, length(p->velocity), mom);
						}

						total_mom += mom;
					}

					ImGui::Text("Momentum: %.1f kg*m/s\nRelative to last frame: %.1f%%\nRelative to saved: %.1f%%",
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
				//ImGui::InputText();
				if(show_all_particles)
				{
					ImGui::BeginTable("Particles", 5,
							ImGuiTableFlags_Resizable
					);
					ImGui::TableSetupColumn("i");
					ImGui::TableSetupColumn("Position (m)");
					ImGui::TableSetupColumn("Velocity (m/s)");
					ImGui::TableSetupColumn("Acceleration (m/s²)");
					ImGui::TableSetupColumn("Mass (kg)");
					ImGui::TableHeadersRow();

					for(int i = 0; i < sim->n_particles(); i++)
					{
						ImGui::TableNextRow();
						particle_t* A = &sim->particles[i];

						ImGui::TableSetColumnIndex(0);
						ImGui::Text("%d", i);
						
						ImGui::TableNextColumn();
						ImGui::Text("(%.1f , %.1f, %.1f)", A->position.x, A->position.y, A->position.z);

						ImGui::TableNextColumn();
						ImGui::Text("(%.1f , %.1f, %.1f) : %.1f", A->velocity.x, A->velocity.y, A->velocity.z, length(A->velocity));

						ImGui::TableNextColumn();
						ImGui::Text("(%.1f , %.1f, %.1f) : %.1f", A->acceleration.x, A->acceleration.y, A->acceleration.z, length(A->acceleration));

						ImGui::TableNextColumn();
						ImGui::Text("%.1f", A->mass);
					}
					ImGui::EndTable();
				}

				ImGui::EndTabItem();
			}

			ImGui::EndTabBar();

			ImGui::End();
		}

		// Only function to call every frame
		bool render_start()
		{
			if(UI::ui_is_shown()) 
			{
				UI::ui_render_start(!clear_screen & UI_RENDER_FLAG_NO_CLEAR);
				SimUI::side_window();
			}

			return UI::ui_is_shown();
		}

		void render_stop()
		{
			UI::ui_render_stop();
		}


		int end()
		{
			UI::ui_cleanup();
			return 0;
		}

};

#endif
