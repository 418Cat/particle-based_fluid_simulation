#ifndef RENDER_H
#define RENDER_H

#include "/home/cat/Documents/Programmation/C++/lagrangian-fluid-simulation/libs/glad/include/glad/glad.h"
#include <GLFW/glfw3.h>

class Render
{
	private:


	public:
		Render()
		{
			if(!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
			{
				throw("Failed to load Glad GL Loader");
			}
		}

};

#endif
