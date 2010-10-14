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

#ifndef IO_WAVEFRONT_H
#define IO_WAVEFRONT_H


#include "Loader.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

namespace arch
{

namespace io 
{

template<typename mesh_t>
class WaveFront 
{
public:
    WaveFront(const std::string &fileName);

    virtual bool load(mesh_t &m);
    virtual bool save(mesh_t &m);

    static bool can_read(const std::string &filename); 
private:
    std::string m_FileName;
};

#include "WaveFront.inl"

}

}

#endif
