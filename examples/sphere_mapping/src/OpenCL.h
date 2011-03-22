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

#ifndef OPENCL_HEADER_H
#define OPENCL_HEADER_H

#include <string>
#include <CL/cl.h>

namespace cl
{

// Round Up Division function
std::size_t roundUp(int group_size, int global_size);

std::string buildlog(cl_program program, cl_device_id device);

cl_int platform_info(cl_platform_id id, cl_platform_info info, std::string &str);

template<typename ListType>
cl_int get_platforms(ListType &platforms)
{
    cl_uint num_platforms; 
    cl_platform_id* clPlatformIDs;
    cl_int ciErrNum;
  
    // Get OpenCL platform count
    ciErrNum = clGetPlatformIDs (0, NULL, &num_platforms);
    if (ciErrNum != CL_SUCCESS)
        return -1000;
    else 
    {
        if(num_platforms == 0)
        {
			std::cerr << "No OpenCL platform found!\n\n";
            return -2000;
        }
        else 
        {
            // if there's a platform or more, make space for ID's
            if ((clPlatformIDs = (cl_platform_id*)malloc(num_platforms * sizeof(cl_platform_id))) == NULL)
            {
				std::cerr << "Failed to allocate memory for cl_platform ID's!\n\n";
                return -3000;
            }

            // get platform info for each platform and trap the NVIDIA platform if found
            ciErrNum = clGetPlatformIDs (num_platforms, clPlatformIDs, NULL);

			for(cl_uint i = 0; i < num_platforms; ++i)
				platforms.push_back( clPlatformIDs[i] );
       
            free(clPlatformIDs);
        }
    }

    return CL_SUCCESS;
}

};

#endif