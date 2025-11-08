#ifndef SIM_UI_H
#define SIM_UI_H

#include <thread>

// render.hpp header must be placed before ui.h to not include opengl twice
#include "util/render/render.hpp"

#include "sim/simulation.hpp"
#include "sim/particle.hpp"

#include <imgui.h>
#include <ui.h>



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
		int col_type_selected = 0;
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
						ImGuiWindowFlags_NoDecoration	|
						ImGuiWindowFlags_NoResize				|
						ImGuiWindowFlags_NoMove
			);

			ImGui::BeginChild("Settings child",
					ImVec2(
						ImGui::GetContentRegionAvail().x*1.,
						ImGui::GetContentRegionAvail().y*0.5));
			ImGui::BeginTabBar("Components settings");

			if(ImGui::BeginTabItem("Sim settings"))
			{
				ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(!sim->settings.run/2., sim->settings.run/2., 0., 1.));
				if(ImGui::Button("Pause") && sim->settings.run) sim->settings.run = false;
				ImGui::PopStyleColor();
				ImGui::SameLine();

				ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(sim->settings.run/2., !sim->settings.run/2., 0., 1.));
				if(ImGui::Button("Resume") && !sim->settings.run)
					sim->run();

				ImGui::PopStyleColor();
				ImGui::SameLine();
				ImGui::DragFloat("Sim speed", &sim->settings.speed, sim->settings.speed/20., 0., 2.);
				//ImGui::InputDouble("Sim speed", &sim->settings.speed, 0.05); // DOUBLE

				ImGui::Spacing();
				//ImGui::DragFloat3("Domain size", (float*)&sim->settings.domain_size, 1., 0.);
				ImGui::Text("Domain size");
				ImGui::InputDouble("X##domain size", &sim->settings.domain_size.x, 0.1);
				ImGui::InputDouble("Y##domain size", &sim->settings.domain_size.y, 0.1);
				ImGui::InputDouble("Z##domain size", &sim->settings.domain_size.z, 0.1);

				ImGui::Spacing();
				if(ImGui::Button("Reset particles"))
				{
					bool was_running = sim->settings.run;
					sim->end();
					sim->spawn_particles_as_rect(false);
					if(was_running) sim->run();
				}

				ImGui::Spacing();
				ImGui::InputInt("N_particles", &n_particles);
				if(ImGui::Button("Reset sim"))
				{
					sim->reset(n_particles);
					render->particles_buffers();
				}

				ImGui::Spacing();

				ImGui::SliderInt("Sim rate", (int*)&sim->settings.hertz, 1, 10000, "%d Hz");
				ImGui::SliderInt("Number of threads", (int*)&sim->settings.n_threads, 1, max_threads);

				ImGui::EndTabItem();
			}

			if(ImGui::BeginTabItem("Collisions"))
			{
				ImGui::Checkbox("Enable", &sim->settings.particles_collisions);

				const char* col_types[] = {"Velocity", "Acceleration"};

				if(ImGui::Combo("Collision type", &col_type_selected, col_types, 2))
					sim->settings.collision_type = sim_state_t::p_collision_type_t(col_type_selected);

				ImGui::Spacing();
				ImGui::SeparatorText("Bounding boxes XYZ");
				ImGui::SliderInt("##BboxesX", (int*)(&sim->settings.n_bounding_boxes_x), 1, (int)floor(sim->settings.domain_size.x/(sim->settings.particle_radius*2.)));
				ImGui::SliderInt("##BboxesY", (int*)(&sim->settings.n_bounding_boxes_y), 1, (int)floor(sim->settings.domain_size.y/(sim->settings.particle_radius*2.)));
				ImGui::SliderInt("##BboxesZ", (int*)(&sim->settings.n_bounding_boxes_z), 1, (int)floor(sim->settings.domain_size.z/(sim->settings.particle_radius*2.)));

				ImGui::Spacing();
				ImGui::Separator();
				//ImGui::SliderFloat("Domain bounciness", &sim->settings.domain_bounciness, 0., 1.5);
				//ImGui::SliderFloat("Particles bounciness", &sim->settings.particles_bounciness, 0., 1.5);
				//ImGui::SliderFloat("Particle radius", &sim->settings.particle_radius, 0.01, 20.);
				ImGui::InputDouble("Domain bounciness", &sim->settings.domain_bounciness, 0.05);
				ImGui::InputDouble("Particles bounciness", &sim->settings.particles_bounciness, 0.05);
				ImGui::InputDouble("Particles radius", &sim->settings.particle_radius, 0.05);

				ImGui::EndTabItem();
			}

			if(ImGui::BeginTabItem("Liquid"))
			{
				ImGui::Checkbox("Enable##liquid_sim", &sim->settings.liquid_sim);

				ImGui::InputDouble("Rest density", &sim->settings.liquid_rest_density, 0.05);
				ImGui::InputDouble("Smoothing radius", &sim->settings.smoothing_length_h, 0.05);
				ImGui::InputDouble("Stiffness constant", &sim->settings.stiffness_constant_k, 0.05);
				ImGui::InputDouble("Kinematic Viscosity", &sim->settings.kinematic_viscosity, 0.05);

				ImGui::EndTabItem();
			}

			if(ImGui::BeginTabItem("Gravity"))
			{
				ImGui::SeparatorText("Domain gravity");
				ImGui::Checkbox("X", &sim->settings.domain_gravity_axis[0]); ImGui::SameLine();
				ImGui::Checkbox("Y", &sim->settings.domain_gravity_axis[1]); ImGui::SameLine();
				ImGui::Checkbox("Z", &sim->settings.domain_gravity_axis[2]);
				//ImGui::SliderFloat3("Domain gravity, m/s", &sim->settings.domain_gravity.x, -40., 40.);
				ImGui::InputDouble("X##domain gravity", &sim->settings.domain_gravity.x, 0.1);
				ImGui::InputDouble("Y##domain gravity", &sim->settings.domain_gravity.y, 0.1);
				ImGui::InputDouble("Z##domain gravity", &sim->settings.domain_gravity.z, 0.1);

				ImGui::Checkbox("Radial", &sim->settings.domain_gravity_radial);

				ImGui::Spacing();
				ImGui::SeparatorText("Particles gravity [Experimental]");
				ImGui::Checkbox("Enabled", &sim->settings.particle_gravity);
				ImGui::SameLine();
				ImGui::Checkbox("Inverse", &sim->settings.particles_gravity_inverse);

				ImGui::Spacing();
				ImGui::SeparatorText("Octree [Experimental]");
				//ImGui::DragFloat("Threshold", &sim->settings.volume_to_distance_threshold, sim->settings.volume_to_distance_threshold/10., 0., 1.e3, "%.5f");
				ImGui::SliderInt("Max depth", (int*)&sim->settings.octree_max_depth, 0, 100);
				ImGui::SliderInt("Min particles in node", (int*)&sim->settings.octree_min_particles_in_node, 1, 200);
				//ImGui::DragFloat("Gravity factor", &sim->settings.particles_gravity_factor, sim->settings.particles_gravity_factor/10., 0., 1.e20);


				ImGui::EndTabItem();
			}

			if(ImGui::BeginTabItem("Display"))
			{
				ImGui::DragFloat("FOV", &camera->fov, camera->fov/10., 1., 100.);

				ImGui::SeparatorText("Particles display");

				ImGui::BeginDisabled(render->show_boxes);
				ImGui::Checkbox("Show##vel", &render->show_vel);
				ImGui::DragFloat("Max Velocity", &render->color_max_vel, render->color_max_vel/10. + 0.0001, 0.);

				ImGui::Spacing();ImGui::Spacing();
				ImGui::Checkbox("Show##accel", &render->show_accel);
				ImGui::DragFloat("Max Acceleration", &render->color_max_accel, render->color_max_accel/10. + 0.0001, 0.);
				ImGui::EndDisabled();

				ImGui::Spacing();ImGui::Spacing();
				ImGui::Checkbox("Show bbox", &render->show_boxes);

				ImGui::Spacing();ImGui::Spacing();
				ImGui::Checkbox("Show density", &render->show_density);


				ImGui::SeparatorText("Domain display");
				ImGui::Checkbox("Show", &render->show_borders);
				ImGui::SameLine();
				ImGui::DragFloat("Domain borders size", &render->border_size, render->border_size/10., 0., 5., "%.1f %%");
				ImGui::Spacing();ImGui::Spacing();

				ImGui::Separator();
				ImGui::Spacing();
				ImGui::Checkbox("Clear screen", &clear_screen);


				ImGui::EndTabItem();
			}
			ImGui::EndTabBar();
			ImGui::EndChild();


			ImGui::BeginTabBar("Info");
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
						Particle* p = &sim->particles[i];
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
						Particle* A = &sim->particles[i];

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
			if(ImGui::BeginTabItem("How to use"))
			{
				ImGui::SeparatorText("Camera movement");
				ImGui::Text("Z: Forward");
				ImGui::Text("S: Backwards");
				ImGui::Text("Q: Left");
				ImGui::Text("D: Right");
				ImGui::Text("A: Up");
				ImGui::Text("E: Down");

				ImGui::Text("Left SHIFT to increase movement speed");
				ImGui::Spacing();
				ImGui::Text("Click + drag to rotate camera [Small bugfix needed]");

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
