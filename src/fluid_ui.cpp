#include "ui.h"

#include "fluid_ui.hpp"

ImVec2 FlUId::window_size = ImVec2(0.f, 0.f);
float FlUId::settings_window_ratio = 1.f/5.f;
GLFWwindow* FlUId::window = UI::window;

void FlUId::side_window()
{
	// Make it a left bar window
	//ImGui::SetNextWindowPos(ImVec2(0., 0.));
	ImGui::SetNextWindowSize(ImVec2(
				FlUId::window_size.x * FlUId::settings_window_ratio,
				FlUId::window_size.y
	));
	
	ImGui::Begin(
			"side_window",
			NULL/*,
				ImGuiWindowFlags_NoDecoration 	|
				ImGuiWindowFlags_NoResize		|
				ImGuiWindowFlags_NoMove*/
	);

	ImGui::Text("Hello");
	ImGui::DragFloat("im a float :D", &FlUId::settings_window_ratio, 0.1f, 0.f, 1.f);

	ImGui::End();

}

void FlUId::rendering_window()
{
	ImGui::SetNextWindowPos(ImVec2(
				FlUId::window_size.x * FlUId::settings_window_ratio,
				0.
	)); 
	
	ImGui::SetNextWindowSize(ImVec2(
				FlUId::window_size.x * (1.-FlUId::settings_window_ratio),
				FlUId::window_size.y
	));

	ImGui::Begin(
			"rendering_window",
			NULL,
				ImGuiWindowFlags_NoDecoration 	| 
				ImGuiWindowFlags_NoResize 		|
				ImGuiWindowFlags_NoMove
	);

	ImGui::Text("This is the rendering window");

	ImGui::End();
}

void FlUId::update_state()
{
	FlUId::window_size = ImGui::GetIO().DisplaySize;
}

bool FlUId::render_start()
{
	if(!UI::ui_is_shown()) return false;

	UI::ui_render_start();

	FlUId::update_state();

	//FlUId::side_window();
	//FlUId::rendering_window();

	return true;
}

void FlUId::render_stop()
{
	UI::ui_render_stop();
}

void FlUId::begin()
{
	UI::ui_init();
	UI::show_ui = true;

	FlUId::window = UI::window;
}

int FlUId::end()
{
	UI::ui_cleanup();
	return 0;
}
