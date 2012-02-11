set(ARCH_PYTHON_SRC
	PyGeometry.cpp
	PyInterface.cpp
	PyLoader.cpp
	PyTypes.cpp
	)

include_directories(${PYTHON_INCLUDE_PATH})
include_directories(../)

add_library(archmind ${ARCH_PYTHON_SRC})
target_link_libraries(archmind ${Boost_LIBRARIES})
