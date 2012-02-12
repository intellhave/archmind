set(SPHERE_MAPPING_SRC
	src/OpenCL.cpp
	src/OpenCLSolver.cpp
	src/Utils.cpp
	src/Main.cpp
	)

include_directories(Boost_INCLUDE_DIRS)
include_directories(src)
include_directories(../../src)

add_executable(sphere_mapping ${SPHERE_MAPPING_SRC})
target_link_libraries(sphere_mapping ${OPENCL_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
