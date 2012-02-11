#Project archmind
cmake_minimum_required(VERSION 2.6)

if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE})
else()
  set(CMAKE_BUILD_TYPE Release)
endif()

project(archmind)

set(
    OPENCL_INC_SEARCH_PATH
    ${OPENCL_INCLUDE_DIR}
    ${OPENCL_DIR}/OpenCL/common/inc
    $ENV{OPENCL_INCLUDE_DIR}
    /usr/local/include
    /usr/include
    #append additional search directories for OpenCL include here
    )

set(
    OPENCL_LIB_SEARCH_PATH
    ${OPENCL_LIBRARY_DIR}
    ${OPENCL_DIR}/lib
    ${OPENCL_DIR}/lib/x86
    $ENV{OPENCL_LIBRARY_DIR}
    $ENV{OPENCL_DIR}/lib
    $ENV{OPENCL_DIR}/lib/x86
    $ENV{LD_LIBRARY_PATH}
    /usr/local/lib64
    /usr/local/lib
    /usr/lib64
    /usr/lib
    #append additional search directories for OpenCL libraries here.
)

if(WIN64)
    set(
        OPENCL_LIB_SEARCH_PATH
        ${OPENCL_LIB_SEARCH_PATH}
        )

    #find_path(OPENCL_INCLUDE_DIR NAMES CL/cl.h PATHS ${OPENCL_INC_SEARCH_PATH})
    find_library(OPENCL_LIBRARIES NAMES OpenCL PATHS ${OPENCL_LIB_SEARCH_PATH})
else(WIN64)
    set(
        OPENCL_LIB_SEARCH_PATH
        ${OPENCL_LIB_SEARCH_PATH}
        )
    
    #find_path(OPENCL_INCLUDE_DIR NAMES CL/cl.h PATHS ${OPENCL_INC_SEARCH_PATH})
    find_library(OPENCL_LIBRARIES NAMES OpenCL PATHS ${OPENCL_LIB_SEARCH_PATH})
endif(WIN64)

find_package(PythonLibs 3.00)

#set(Boost_USE_STATIC_LIBS ONN)
#set(Boost_USE_MULTITHREADED ON)
#set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.40 COMPONENTS program_options python)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG")

if(CMAKE_BUILD_TYPE MATCHES "Release")
  message("Release build!")
  #add_subdirectory(src/Python)
  add_subdirectory(examples/cpp "Release/examples/cpp")
	add_subdirectory(examples/sphere_mapping "Release/examples/sphere_mapping")
else(CMAKE_BUILD_TYPE MATCHES "Release")
  #add_subdirectory(src/Python)
  if(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    add_subdirectory(examples/cpp "RelWithDebInfo/examples/cpp")
    add_subdirectory(examples/sphere_mapping "RelWithDebInfo/examples/sphere_mapping")
  else(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    if(CMAKE_BUILD_TYPE MATCHES "Debug")
    add_subdirectory(examples/cpp "Debug/examples/cpp")
    add_subdirectory(examples/sphere_mapping "Debug/examples/sphere_mapping")
  else(CMAKE_BUILD_TYPE MATCHES "Debug")
    add_subdirectory(examples/cpp "MinSizeRel/examples/cpp")
    add_subdirectory(examples/sphere_mapping "MinSizeRel/examples/sphere_mapping")
  endif(CMAKE_BUILD_TYPE MATCHES "Debug")
  endif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
endif(CMAKE_BUILD_TYPE MATCHES "Release")
