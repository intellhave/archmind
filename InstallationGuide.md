# Installing the Blender Tools #
In order to install the Blender tools, the scripts and the precompiled python library should be copied into the Blender scripts folder located at .blender/scripts.
Alternatively, if there is a working python installation in the system, the python library may be copied into the python library directory.

**Note:** _Most of the Blender scripts do not require a full python installation._

# Using the c++ Interface #

The c++ interface is header-only at the moment, except from the python bindings.
Therefore, only the corresponding header files should be included in your project.
Nevertheless, the [boost library](http://www.boost.org/) is required.

# Windows Python Library Compilation #

In order to build the python library, a working python installation is required along with the boost library.
There is a Microsoft Visual Studio solution in the msvc folder and a cmake configuration file.

# Linux Python Library Compilation #

To build the project and the examples, create a build directory and run cmake and then make from there. For example:
```
mkdir build
cd build
cmake .. -G "Unix Makefiles" && make
```