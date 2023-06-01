cmake_minimum_required(VERSION 3.10)

project(Bachelorprojekt)

add_subdirectory(sources)

add_executable(main
	${PROJECT_SOURCE_DIR}/main.cpp
)

target_link_libraries(main PRIVATE MolekylBibliotek)

