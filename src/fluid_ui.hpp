#ifndef FLUID_UI_H
#define FLUID_UI_H

#include "imgui.h"

// Funny name with Fluid and UI :D
class FlUId
{
	public:
		static ImVec2 window_size;
		static float settings_window_ratio;

		// To call once before the main loop
		static void begin();

		// Update attributes
		static void update_state();
		
		// Render imgui components
		static void side_window();
		static void rendering_window();

		// Only function to call every frame
		static bool render();

		static int end();

};

#endif
