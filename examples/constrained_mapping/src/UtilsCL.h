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

#ifndef UTILS_CL_H
#define UTILS_CL_H

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

namespace arch
{

namespace clutils
{
void init_best_device(int &platform, int &device, int device_type = CL_DEVICE_TYPE_ALL);
int round_up(int group_size, int global_size);
unsigned int next_pow2(unsigned int x);
std::string error( cl_int code );
}

}

#endif
