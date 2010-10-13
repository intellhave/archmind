/*
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

#include "../Geometry/Geometry.h"
#include "../Geometry/Algorithms.h"
#include "../Python/PyInterface.h"

typedef geometry::mesh<> mesh_t;
typedef mesh_t::vertex_t vertex_t;
typedef mesh_t::vertex_ptr_t vertex_ptr_t;
typedef mesh_t::face_t face_t;
typedef mesh_t::face_ptr_t face_ptr_t;

#include <deque>

int main(int argc, char **argv)
{
	using boost::shared_ptr;
    using namespace geometry;

	mesh_t mesh;

	if(argc > 1)
    {
        python::Interface pyEngine(argc,argv);

        if( !pyEngine.run_script(argv[1]) )
        {
            std::cout << "Python err : " << pyEngine.last_error() << std::endl;
        }
    }

	return 0;
}
