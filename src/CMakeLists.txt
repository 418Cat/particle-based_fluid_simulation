# Building UI part
add_library(
	fluid
	fluid_ui.hpp
	fluid_ui.cpp
)

target_link_libraries(
	fluid PUBLIC
	ui
	imgui
)



# Building Simulation part
add_library(
	simulation
	simulation.cpp
	simulation.hpp
)

target_link_libraries(
	simulation PUBLIC
	glm
)

# Building render engine
add_library(
	render
	render.hpp
)

target_link_libraries(
	render PUBLIC
	glad
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
	fluid
	simulation
	render
)
