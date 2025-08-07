# Building UI part
add_library(
	sim_ui	
	sim_ui.hpp
)

target_link_libraries(
	sim_ui PUBLIC
	ui
	imgui
)

# Building Simulation part
add_library(
	simulation
	simulation.hpp
)

target_link_libraries(
	simulation PUBLIC
	glm
)

set_target_properties(
	simulation PROPERTIES
	LINKER_LANGUAGE CXX
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
	render.hpp
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
