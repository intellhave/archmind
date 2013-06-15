set(CONSTRAINED_MAPPING_SRC
    src/OpenCLSolver.cpp
    src/OpenCLBlas.cpp
    src/UtilsCL.cpp
	  src/Main.cpp
	)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(src)
include_directories(../../src)
include_directories(../../include)

if(MSVC)
    link_directories(${Boost_LIBRARY_DIRS})
    add_executable(constrained_mapping ${CONSTRAINED_MAPPING_SRC})
    target_link_libraries(constrained_mapping ${OPENCL_LIBRARIES})
else()
    add_executable(constrained_mapping ${CONSTRAINED_MAPPING_SRC})
    target_link_libraries(constrained_mapping ${OPENCL_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
endif()

