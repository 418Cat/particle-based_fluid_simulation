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
	simulation.hpp
)

target_link_libraries(
	simulation PUBLIC
	glm
	shader
)

set_target_properties(
	simulation PROPERTIES
	LINKER_LANGUAGE CXX
)

target_include_directories(
	simulation PUBLIC
	${CMAKE_SOURCE_DIR}/libs/glad/include
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
	shader
)

target_include_directories(
	render PUBLIC
	${GLAD_INCLUDE_DIR}
	${CMAKE_SOURCE_DIR}/libs/utils/shader/
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
