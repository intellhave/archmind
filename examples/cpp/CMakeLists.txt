include_directories(Boost_INCLUDE_DIRS)
include_directories(../../src)

add_executable(smooth smooth.cpp)
add_executable(stl_export stl_export.cpp)
add_executable(tetra tetra.cpp)

