add_subdirectory(utils/UI)

add_library(
	imgui
	imgui/imgui.cpp
	imgui/imgui_demo.cpp
	imgui/imgui_draw.cpp
	imgui/imgui_tables.cpp
	imgui/imgui_widgets.cpp

	imgui/imgui.h

	imgui/backends/imgui_impl_opengl3.cpp
	imgui/backends/imgui_impl_opengl3.h
	imgui/backends/imgui_impl_opengl3_loader.h

	imgui/backends/imgui_impl_glfw.cpp
	imgui/backends/imgui_impl_glfw.h
)

target_include_directories(
	imgui PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/imgui
	${CMAKE_CURRENT_SOURCE_DIR}/imgui/backends
)

target_link_libraries(
	imgui PUBLIC
	libglfw.so.3
	libGL.so
)

