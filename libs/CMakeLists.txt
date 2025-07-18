add_subdirectory(utils/UI)

# _____ BUILDING  IMGUI _____
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
# ____________________


# Making a convenient lib for GLM
# Doesn't actually need to build, it's header only
add_library(
	glm
	glm/glm/glm.hpp
)

# It's some .h files so cmake
# doesn't know if it's C or C++
set_target_properties(
	glm PROPERTIES
	LINKER_LANGUAGE CXX # Set to C++
)

target_include_directories(
	glm PUBLIC
	glm/glm
)
# ____________________


# _____ Glad loader _____
add_library(
	glad
	glad/src/glad.c
	glad/include/glad/glad.h
	glad/include/KHR/khrplatform.h
)

target_include_directories(
	glad PUBLIC
	glad/include
)

target_link_libraries(
	glad PUBLIC
	libglfw.so.3
	libGL.so
)
# ____________________
