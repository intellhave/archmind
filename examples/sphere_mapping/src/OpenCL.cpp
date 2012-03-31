/*
  Parallel Computation of Spherical Parameterizations for Mesh Analysis
  Copyright (C) 2011 Athanasiadis Theodoros and Fudos Ioannis

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
*/

#include "OpenCL.h"
#include <iostream>

#pragma comment (lib, "OpenCL.lib") 

using namespace cl;

// Round Up Division function
int cl::round_up(int group_size, int global_size) 
{
    int r = global_size % group_size;
    if(r == 0) 
    {
        return global_size;
    } else 
    {
        return global_size + group_size - r;
    }
}

unsigned int cl::next_pow2(unsigned int x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

std::string cl::buildlog(cl_program program, cl_device_id device)
{
	std::size_t strsize;

	clGetProgramBuildInfo(program,device,CL_PROGRAM_BUILD_LOG,0,NULL,&strsize); 
	char *log = new char[strsize+1];
	clGetProgramBuildInfo(program,device,CL_PROGRAM_BUILD_LOG,strsize,log,&strsize);  
	std::string retstr(log);
	delete [] log;
	return retstr;
}

cl_int cl::platform_info(cl_platform_id id, cl_platform_info info, std::string &str)
{
	 char chBuffer[1024];
	 if( clGetPlatformInfo (id, info, 1024, &chBuffer, NULL) != CL_SUCCESS )
		 return -1000;
	 else
	 {
		 str.assign(chBuffer);
		 return CL_SUCCESS;
	 }
}
