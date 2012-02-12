set(ARCH_PYTHON_SRC
	PyGeometry.cpp
	PyInterface.cpp
	PyLoader.cpp
	PyTypes.cpp
	)

include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(../)

if(MSVC)
    link_directories(${Boost_LIBRARY_DIRS})
    add_library(archmind ${ARCH_PYTHON_SRC} )
else()
    add_library(archmind ${ARCH_PYTHON_SRC} )
    target_link_libraries(sphere_mapping ${Boost_PYTHON_LIBRARY})
endif()

