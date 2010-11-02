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

template<typename mesh_t>
class reader_dispatch
{
public:
    reader_dispatch(const std::string &filename, mesh_t &mesh, bool &finish) : 
        m_filename(filename), m_mesh(mesh), m_finish(finish) {}
    
    template< typename U > void operator()(U &x)
    {
        if( !m_finish && x.can_read(m_filename) )
            m_finish = x.read(m_filename,m_mesh);
    }

private:
    std::string m_filename;
    mesh_t &m_mesh;
    bool &m_finish;
};

template<typename mesh_t>
class writer_dispatch
{
public:
    writer_dispatch(const std::string &filename, mesh_t &mesh, bool &finish) : 
        m_filename(filename), m_mesh(mesh), m_finish(finish) {}
    
    template< typename U > void operator()(U &x)
    {
        if( !m_finish && x.can_read(m_filename) )
            m_finish = x.write(m_filename,m_mesh);
    }

private:
    std::string m_filename;
    mesh_t &m_mesh;
    bool &m_finish;
};

template<typename mesh_t>
struct reader_traits;

template<typename mesh_t>
bool load_from_file(const std::string &filename, mesh_t &mesh)
{
    bool success = false;
    boost::mpl::for_each<typename reader_traits<mesh_t>::type>( reader_dispatch<mesh_t>( filename,mesh,success) );
    return success;
}

template<typename mesh_t>
struct writer_traits;

template<typename mesh_t>
bool save_to_file(const std::string &filename, mesh_t &mesh)
{
    bool success = false;
    boost::mpl::for_each<typename writer_traits<mesh_t>::type>( writer_dispatch<mesh_t>( filename,mesh,success) );
    return success;
}

