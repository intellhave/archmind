/*
Archmind Non-manifold Geometric Kernel
Copyright (C) 2010 Athanasiadis Theodoros

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

#include "UtilsCL.h"
#include <fstream>
#include <iostream>

// Round Up Division function
int arch::clutils::round_up(int group_size, int global_size) 
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

unsigned int arch::clutils::next_pow2(unsigned int x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}


void arch::clutils::init_best_device(int &platform, int &device, int device_type)
{
    using namespace cl;
    static int sel_platform = -1;
    static int sel_device = -1;
    int sel_compute_units = 0;

    if( sel_platform == -1 && sel_device == -1 ) {
        std::vector<Platform> platforms;
        std::vector<cl::Device> devices;
        Platform::get(&platforms);

        for( std::size_t i = 0; i < platforms.size(); ++i ) {
            cl_context_properties properties[] = 
            { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[i])(), 0};

            cl::Context context(CL_DEVICE_TYPE_ALL, properties); 
            devices = context.getInfo<CL_CONTEXT_DEVICES>();

            for( std::size_t j = 0; j < devices.size(); ++j ) 
            {
                
                if( devices[j].getInfo<CL_DEVICE_TYPE>() & device_type )
                {
                    int compute_units = devices[j].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
                    if( compute_units > sel_compute_units ) {
                        sel_platform = i;
                        sel_device = j;
                        sel_compute_units = compute_units;
                       
                        std::cout << "Platform name       : " << platforms[i].getInfo<CL_PLATFORM_NAME>() << "\n";
                        std::cout << "Platform vendor     : " << platforms[i].getInfo<CL_PLATFORM_VENDOR>() << "\n";
                        std::cout << "Platform version    : " << platforms[i].getInfo<CL_PLATFORM_VERSION>() << "\n";
                        std::cout << "Device name         : " << devices[j].getInfo< CL_DEVICE_NAME>() << "\n";
                        std::cout << "Device compute units: " << compute_units << "\n";
                        std::cout << "Device global mem   : " << devices[j].getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << " bytes \n";
                        std::cout << "Device max workgroup: " << devices[j].getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << "\n";
                        std::cout << "Device clock freq   : " << devices[j].getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << " Mhz \n";
                        std::cout << "Platform extensions : " << platforms[i].getInfo<CL_PLATFORM_EXTENSIONS>() << "\n";
                        
                    }
                }
            }
        }

        if( sel_platform == -1 ) throw std::runtime_error("Failed to find OpenCL platform\n");
    } 
   
    platform = sel_platform;
    device = sel_device;
}

std::string arch::clutils::error(cl_int err)
{
    switch(err) {                                                                            
        case CL_INVALID_VALUE                   : return std::string( "CL_INVALID_VALUE"              );
        case CL_INVALID_DEVICE_TYPE             : return std::string( "CL_INVALID_DEVICE_TYPE"        );
        case CL_INVALID_PLATFORM                : return std::string( "CL_INVALID_PLATFORM"           );
        case CL_INVALID_DEVICE                  : return std::string( "CL_INVALID_DEVICE"             );
        case CL_INVALID_CONTEXT                 : return std::string( "CL_INVALID_CONTEXT"            );
        case CL_INVALID_QUEUE_PROPERTIES        : return std::string( "CL_INVALID_QUEUE_PROPERTIES"   );
        case CL_INVALID_COMMAND_QUEUE           : return std::string( "CL_INVALID_COMMAND_QUEUE"      );
        case CL_INVALID_HOST_PTR                : return std::string( "CL_INVALID_HOST_PTR"           );
        case CL_INVALID_MEM_OBJECT              : return std::string( "CL_INVALID_MEM_OBJECT"         );
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR : return std::string( "CL_INVALID_IMAGE_FORMAT_DESCR" );
        case CL_INVALID_IMAGE_SIZE              : return std::string( "CL_INVALID_IMAGE_SIZE"         );
        case CL_INVALID_SAMPLER                 : return std::string( "CL_INVALID_SAMPLER"            );
        case CL_INVALID_BINARY                  : return std::string( "CL_INVALID_BINARY"             );
        case CL_INVALID_BUILD_OPTIONS           : return std::string( "CL_INVALID_BUILD_OPTIONS"      );
        case CL_INVALID_PROGRAM                 : return std::string( "CL_INVALID_PROGRAM"            );
        case CL_INVALID_PROGRAM_EXECUTABLE      : return std::string( "CL_INVALID_PROGRAM_EXECUTABLE" );
        case CL_INVALID_KERNEL_NAME             : return std::string( "CL_INVALID_KERNEL_NAME"        );
        case CL_INVALID_KERNEL_DEFINITION       : return std::string( "CL_INVALID_KERNEL_DEFINITION"  );
        case CL_INVALID_KERNEL                  : return std::string( "CL_INVALID_KERNEL"             );
        case CL_INVALID_ARG_INDEX               : return std::string( "CL_INVALID_ARG_INDEX"          );
        case CL_INVALID_ARG_VALUE               : return std::string( "CL_INVALID_ARG_VALUE"          );
        case CL_INVALID_ARG_SIZE                : return std::string( "CL_INVALID_ARG_SIZE"           );
        case CL_INVALID_KERNEL_ARGS             : return std::string( "CL_INVALID_KERNEL_ARGS"        );
        case CL_INVALID_WORK_DIMENSION          : return std::string( "CL_INVALID_WORK_DIMENSION"     );
        case CL_INVALID_WORK_GROUP_SIZE         : return std::string( "CL_INVALID_WORK_GROUP_SIZE"    );
        case CL_INVALID_WORK_ITEM_SIZE          : return std::string( "CL_INVALID_WORK_ITEM_SIZE"     );
        case CL_INVALID_GLOBAL_OFFSET           : return std::string( "CL_INVALID_GLOBAL_OFFSET"      );
        case CL_INVALID_EVENT_WAIT_LIST         : return std::string( "CL_INVALID_EVENT_WAIT_LIST"    );
        case CL_INVALID_EVENT                   : return std::string( "CL_INVALID_EVENT"              );
        case CL_INVALID_OPERATION               : return std::string( "CL_INVALID_OPERATION"          );
        case CL_INVALID_GL_OBJECT               : return std::string( "CL_INVALID_GL_OBJECT"          );
        case CL_INVALID_BUFFER_SIZE             : return std::string( "CL_INVALID_BUFFER_SIZE"        );
        case CL_INVALID_MIP_LEVEL               : return std::string( "CL_INVALID_MIP_LEVEL"          );
        case CL_INVALID_GLOBAL_WORK_SIZE        : return std::string( "CL_INVALID_GLOBAL_WORK_SIZE"   );
        case CL_INVALID_PROPERTY                : return std::string( "CL_INVALID_PROPERTY"           );
        default : return std::string("Unknown");
    }
}
