# Building UI part
add_library(
	sim_ui	
	util/ui/sim_ui.hpp
)

target_link_libraries(
	sim_ui PUBLIC
	ui
	imgui
)

target_include_directories(
	sim_ui PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}
)

# Building Simulation part
add_library(
	simulation
	sim/simulation.hpp
	sim/particle.hpp
	sim/collision_info.hpp
)

target_link_libraries(
	simulation PUBLIC
	glm
)

set_target_properties(
	simulation PROPERTIES
	LINKER_LANGUAGE CXX
)

target_include_directories(
	simulation PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}
)

# Building render engine
file( # Copy shaders to corresponding bin dir
	COPY
	${CMAKE_CURRENT_SOURCE_DIR}/shaders/
	DESTINATION
	${CMAKE_CURRENT_BINARY_DIR}/shaders/
)

add_library(
	render
	util/render/render.hpp
	util/render/camera.hpp
)

target_include_directories(
	render PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(
	render PUBLIC
	glad
	shader
)

set_target_properties(
	render PROPERTIES
	LINKER_LANGUAGE CXX
)

# Main program using both libraries above
add_executable(
	main
	main.cpp
)

target_link_libraries(
	main PUBLIC
	sim_ui
	simulation
	render
)
