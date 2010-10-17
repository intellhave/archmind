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

#ifndef IO_IO_H
#define IO_IO_H

#include "Base.h"
#include "WaveFront.h"      //Wavefront file importer/exporter
#include "OFF.h"            //OFF file exporter

namespace arch
{

namespace io
{

template<typename mesh_t>
struct reader_traits
{
    typedef boost::mpl::vector< 
        WaveFront<mesh_t>,          //register Wavefront importer
        OFF<mesh_t>                 //register OFF importer
    > type;
};

template<typename mesh_t>
struct writer_traits
{
    typedef boost::mpl::vector< 
        WaveFront<mesh_t>,          //register Wavefront exporter
        OFF<mesh_t>                 //register OFF exporter
    > type;
};

}

}

#endif
