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

add_executable(
	main
	main.cpp
)

target_link_libraries(
	main PUBLIC
	fluid
)
