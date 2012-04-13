set(ARCH_PYTHON_SRC
	PyGeometry.cpp
	PyInterface.cpp
	PyLoader.cpp
	PyTypes.cpp
	)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(../)

if(MSVC)
  link_directories(${Boost_LIBRARY_DIRS})
  add_library(archmind STATIC ${ARCH_PYTHON_SRC} )
  target_link_libraries(archmind ${Boost_PYTHON_LIBRARY})
else()
  add_library(archmind ${ARCH_PYTHON_SRC} )
  target_link_libraries(archmind ${Boost_PYTHON_LIBRARY})
endif()

