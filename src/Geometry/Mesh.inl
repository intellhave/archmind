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

template<typename Traits>
mesh<Traits>::mesh()
{
#ifndef NDEBUG
    //std::clog << "Mesh()" << std::endl;
#endif
}

template<typename Traits>
mesh<Traits>::~mesh()
{
#ifndef NDEBUG
    //std::clog << "~Mesh()" << std::endl;
#endif

#if 0
    //Erase edge cyclic references
    for( face_iterator_t f = faces_begin(); f != faces_end(); ++f )
        (*f)->Edges.clear();

    for( edge_iterator_t e = edges_begin(); e != edges_end(); ++e )
    {
        (*e)->Faces.clear();
        (*e)->Vertices[0] = (*e)->Vertices[1] = vertex_ptr_t();
    }

    for( vertex_iterator_t v = verts_begin(); v != verts_end(); ++v )
        (*v)->Edges.clear();
#endif   
}

//Iterators
template<typename Traits>
typename mesh<Traits>::vertex_iterator_t mesh<Traits>::verts_begin()const
{
    return Vertices.begin();
}

template<typename Traits>
typename mesh<Traits>::vertex_iterator_t mesh<Traits>::verts_end()const
{
    return Vertices.end();
}

template<typename Traits>
typename mesh<Traits>::vertex_range_t mesh<Traits>::verts()const
{
    return vertex_range_t(Vertices.begin(), Vertices.end());
}

template<typename Traits>
std::size_t mesh<Traits>::verts_size()const
{
    return Vertices.size();
}

template<typename Traits>
typename mesh<Traits>::edge_range_t mesh<Traits>::edges()const
{
    return edge_range_t(Edges.begin(), Edges.end());
}
    
template<typename Traits>
typename mesh<Traits>::edge_iterator_t mesh<Traits>::edges_begin()const
{
    return Edges.begin();
}

template<typename Traits>
typename mesh<Traits>::edge_iterator_t mesh<Traits>::edges_end()const
{
    return Edges.end();
}

template<typename Traits>
std::size_t mesh<Traits>::edges_size()const
{
    return Edges.size();
}

template<typename Traits>
typename mesh<Traits>::face_range_t mesh<Traits>::faces()const
{
    return face_range_t(Faces.begin(), Faces.end());
}
    
template<typename Traits>
typename mesh<Traits>::face_iterator_t mesh<Traits>::faces_begin()const
{
    return Faces.begin();
}

template<typename Traits>
typename mesh<Traits>::face_iterator_t mesh<Traits>::faces_end()const
{
    return Faces.end();
}

template<typename Traits>
std::size_t mesh<Traits>::faces_size()const
{
    return Faces.size();
}

template<typename Traits>
bool mesh<Traits>::remove_vertex( vertex_ptr_t v )
{
    
    if( v->id() == NO_ID )      //check if this is a valid mesh vertex
    {
        //print a message ?
        return true;
    }
    else
    {
        delete_vertex( v );
    }

    return false;
}

template<typename Traits>
void mesh<Traits>::swap_vertex( vertex_ptr_t v0, vertex_ptr_t v1 )
{
	uid_t id0 = v0->id();
	uid_t id1 = v1->id();

	std::swap( Vertices[id0], Vertices[id1] );

	v0->set_id(id1);
	v1->set_id(id0);
}

template<typename Traits>
bool mesh<Traits>::remove_edge( edge_ptr_t e )
{
    if( e->id() == NO_ID )
    {
        return false;
    }
    else
    {
        typename edge_t::vertex_iterator_t v;

        //remove the references from the vertices
        for( v = e->verts_begin(); v != e->verts_end(); ++v )
        {
            for( std::size_t i = 0; i < (*v)->Edges.size(); )
                if( (*v)->Edges[i].lock() == e )
                {
                    std::swap( (*v)->Edges[i], (*v)->Edges.back() );
                    (*v)->Edges.pop_back();
                }
                else
                    ++i;

            //check if we also have to delete the vertex
            if( (*v)->Edges.empty() )
                delete_vertex( *v );
        }

        //printf("Deleting edge (%d)\n", e->unique_id() );
        delete_edge( e );
    }

    return true;
}

template<typename Traits>
bool mesh<Traits>::remove_face( face_ptr_t f )
{
    if( f->id() == NO_ID )
    {
        return false;
    }
    else
    {
        typename face_t::edge_iterator_t e;

        //Remove the face from the edges that reference it
        for( e = f->edges_begin(); e != f->edges_end(); ++e )
        {
            for( std::size_t i = 0; i < (*e)->Faces.size(); )
                if( (*e)->Faces[i].lock() == f )
                {
                    std::swap( (*e)->Faces[i], (*e)->Faces.back() );
                    (*e)->Faces.pop_back();
                }
                else
                    ++i;

            //check if we also have to delete the edge
            if( (*e)->Faces.empty() )
                remove_edge( *e );
        }

        delete_face( f );
    }

    return true;
}

template<typename Traits>
void mesh<Traits>::delete_vertex( vertex_ptr_t v )
{
    uid_t ID = v->id();
    v->set_id( NO_ID );

    //check if the ID is valid
    assert(ID < Vertices.size());

    //reassign the vertex id
    Vertices.back()->set_id( ID );

    //O(1) removal from vector
    std::swap( Vertices[ID], Vertices.back() );
    Vertices.pop_back();
}

template<typename Traits>
void mesh<Traits>::delete_edge( edge_ptr_t e )
{
    uid_t ID = e->id();
    e->set_id( NO_ID );
    
    //check if the ID is valid
    assert(ID < Edges.size());

    //reassign the edge id
    Edges.back()->set_id( ID );

    std::swap( Edges[ID], Edges.back() );
    Edges.pop_back();
}

template<typename Traits>
void mesh<Traits>::delete_face( face_ptr_t f )
{
    uid_t ID = f->id();
    f->set_id( NO_ID );

    //check if the ID is valid
    assert(ID < Faces.size());

    //reassign the face id
    Faces.back()->set_id( ID );
    
    std::swap( Faces[ID], Faces.back() );
    Faces.pop_back();
}

template<typename Traits>
uid_t mesh<Traits>::add_vertex( vertex_ptr_t v )
{
    //Check if the vertex has allready been inserted in the mesh
    if( v->id() != NO_ID  )
    {
        if( v->id() < Vertices.size() && Vertices[v->id()] == v )
            return v->id();
    }

    v->set_id( Vertices.size() );
    Vertices.push_back( v );
    return v->id();
}

template<typename Traits>
uid_t mesh<Traits>::add_face( face_ptr_t f )
{
    typename face_t::vertex_iterator_t v;
    typename vertex_t::edge_iterator_t ve;

    //dont add duplicate faces
    if( f->id() != NO_ID )
        return f->id();

    std::size_t e_id = 0;

    //Insert the polygon vertices
    for( v = f->verts_begin(); v != f->verts_end(); ++e_id, ++v )
    {
        add_vertex( *v );

        //Check if we add a new edge or reference an old one
        for( ve = (*v)->edges_begin(); ve != (*v)->edges_end(); ++ve )
        {
            if( (*ve)->v0() == f->Edges[e_id]->v0() &&
                (*ve)->v1() == f->Edges[e_id]->v1() )
            {
                f->Edges[e_id] = *ve;       //assign the correct edge
                (*ve)->Faces.push_back( f );        //add a reference to the new face
                break;
            }
        }

        //check if we failed to find an old edge
        if( ve == (*v)->edges_end() )
        {
            //Add a face reference to the edge
            f->Edges[e_id]->Faces.push_back( f );
            f->Edges[e_id]->set_id( Edges.size() );
            
            Edges.push_back( f->Edges[e_id] );  //add the new edge
            
            //add the edge to the vertices
            f->Edges[e_id]->Vertices[0]->Edges.push_back( f->Edges[e_id] );
            f->Edges[e_id]->Vertices[1]->Edges.push_back( f->Edges[e_id] );

            //add_vertex( f->Edges[e_id]->Vertices[0] );
            //add_vertex( f->Edges[e_id]->Vertices[1] );
        }
    }

    //Add the face at the end
    f->set_id( Faces.size() );
    Faces.push_back( f );

    return f->id();
}

template<typename Traits>
typename mesh<Traits>::vertex_ptr_t 
mesh<Traits>::split_edge( edge_ptr_t e, real_t t, bool triangulate )
{
    if( e->id() == NO_ID )      //check if this edge is part of the mesh
        return e->verts()[0];

    typename vertex_t::point_t v0 = e->verts()[0]->point();      //edge first point
    typename vertex_t::point_t v1 = e->verts()[1]->point();      //edge second point

    vertex_ptr_t sv( new vertex_t(v0 + (v1 - v0) * t) );

    add_vertex( sv );

    //store the old list of faces incident to this edge
    std::vector< face_ptr_t > faces( e->faces_begin(), e->faces_end() );
    
    //split all faces incident to this edge
    for( std::size_t i = 0; i < faces.size(); ++i )
    {
        if( !triangulate )
        {
            //copy old poly vertices
            std::vector< vertex_ptr_t > verts( faces[i]->verts_begin(), faces[i]->verts_end() );
    
            //find the edge position on the face
            std::size_t pos = std::distance( 
                faces[i]->edges_begin(),
                std::find( faces[i]->edges_begin(), faces[i]->edges_end(), e ) );

            //insert the new vertex between the vertices of the edge
            verts.insert( verts.begin() + pos + 1, sv );

            //face_ptr_t nf(new face_t( *faces[i] ));
            //nf->assign( verts.begin(), verts.end() );
            
            face_ptr_t nf(new face_t( verts.begin(), verts.end() ));
            add_face( nf );
        }
        else
        {
            boost::array<vertex_ptr_t,3> verts; 

            for( typename face_t::edge_iterator_t ei = faces[i]->edges_begin(); ei != faces[i]->edges_end(); ++ei )
            {
                if( (*ei) != e )
                {
                    if( faces[i]->edge_cw( *ei ) )
                    {
                        verts[0] = (*ei)->verts()[0];
                        verts[1] = (*ei)->verts()[1];
                        verts[2] = sv; 
                    }
                    else
                    {
                        verts[0] = (*ei)->verts()[1];
                        verts[1] = (*ei)->verts()[0];
                        verts[2] = sv; 
                    }
                 
                    //face_ptr_t nf(new face_t( *faces[i] ));
                    //nf->assign( verts.begin(), verts.end() );
                    face_ptr_t nf(new face_t( verts.begin(), verts.end() ));
                    
                    add_face( nf );
                }
            }
        }

        remove_face( faces[i] );
    }

    return sv;
}

template<typename Traits>
bool mesh<Traits>::join_edge( edge_ptr_t e, vertex_ptr_t v )
{
    if( e->id() == NO_ID )
        return false;

    //store the other vertex of the edge
    vertex_ptr_t ov = (e->verts()[0] == v) ? e->verts()[1] : e->verts()[0];

    //keep a copy of old edges faces
    std::vector< face_ptr_t > e_faces( e->faces_begin(), e->faces_end() );
    std::vector< face_ptr_t > v_faces( ov->faces_begin(), ov->faces_end() );

    for( std::size_t i = 0; i < v_faces.size(); ++i )
    {
        //check if the face is incident to the edge collapsed
        if( std::find( e_faces.begin(), e_faces.end(), v_faces[i] ) == e_faces.end() )
        {
            //build a new face with the new vertex instead of the old one
            std::vector< vertex_ptr_t > verts;
            for( typename face_t::vertex_iterator_t vi = v_faces[i]->verts_begin(); vi != v_faces[i]->verts_end(); ++vi )
                verts.push_back( (*vi == ov) ? v : *vi ); 

            //face_ptr_t nf(new face_t( *v_faces[i] ));
            //nf->assign( verts.begin(), verts.end() );
            face_ptr_t nf(new face_t( verts.begin(), verts.end() ));
            add_face( nf ); 
        }

        //delete the old face
        remove_face( v_faces[i] );
    }

    return true;
}

template<typename Traits>
typename mesh<Traits>::edge_ptr_t 
mesh<Traits>::split_face( face_ptr_t f, vertex_ptr_t v0, vertex_ptr_t v1 )
{
    if( f->id() == NO_ID )
        //return a dummy edge
        return edge_ptr_t( new edge_t(v0,v1) );

    std::vector< vertex_ptr_t > verts( f->verts_begin(), f->verts_end() );
    std::vector< vertex_ptr_t > poly0;
    std::vector< vertex_ptr_t > poly1;
    std::size_t pos0 = 0, pos1 = 0, i;

    //find the positions of the vertices
    for( i = 0; i < verts.size(); ++i )
    {
        if( verts[i] == v0 ) pos0 = i;
        else if( verts[i] == v1 ) pos1 = i;
    }

    for( i = pos0; i != (pos1+1)%verts.size(); i = (i+1) % verts.size() )
        poly0.push_back( verts[i] );
    
    for( i = pos1; i != (pos0+1)%verts.size(); i = (i+1) % verts.size() )
        poly1.push_back( verts[i] );
    
    //build new faces and copy the original face data
    //face_ptr_t nf0( new face_t(*f) );
    //nf0->assign( poly0.begin(), poly0.end() );
    face_ptr_t nf0( new face_t( poly0.begin(), poly0.end() ) );
    add_face( nf0 );

    //face_ptr_t nf1( new face_t(*f) );
    //nf1->assign( poly1.begin(), poly1.end() );
    face_ptr_t nf1( new face_t( poly1.begin(), poly1.end() ) );
    add_face( nf1 );

#ifndef NDEBUG
    using std::abs;
    if( nf0->verts_size() < 3 || nf1->verts_size() < 3 )
        std::clog << "split_face : warning : degenerated face created\n";
#endif

    //remove the original face
    remove_face( f );

    //the new edge is the last edge of each new poly
    return nf1->Edges.back();
}

template<typename Traits>
typename mesh<Traits>::face_ptr_t 
mesh<Traits>::join_face( face_ptr_t f0, face_ptr_t f1 )
{
    using std::abs;

    //copy the vertices of both faces
    std::vector< vertex_ptr_t > verts0(f0->verts_begin(),f0->verts_end());
    std::vector< vertex_ptr_t > verts1(f1->verts_begin(),f1->verts_end());
    std::deque< vertex_ptr_t > diff_verts;

    typedef std::pair< std::size_t, std::size_t > v_t;
    std::vector< v_t > pairs;

    //find the common vertices
    for( std::size_t i = 0; i < verts0.size(); ++i )
        for( std::size_t j = 0; j < verts1.size(); ++j )
            if( verts0[i] == verts1[j] ) 
                pairs.push_back( v_t(i,j) );

    //check if the polygons share only one edge
    if( pairs.size() != 2 ) 
    {
#ifndef NDEBUG
        std::clog << "join_face : warning : faces share more than one edge\n";
#endif
        return f0;
    }

    //the vertices of the common edge (v0,v1) are visited either a) ...v0v1.. or b) v1....v0 in the first poly
    //the vertices of the other poly must be inserted after v0  
    if( pairs[1].first != pairs[0].first + 1 )      //check for case b)
        std::swap( pairs[1].first, pairs[0].first );
    
    /*
         v3
        /  \
       / <- \
   v1 /______\v0   
      \  ->  /
       \    / 
        \  /
         v2

    v1 v3 v0
    |     |
    |    / 
    v1 v0 v2

    correct orientation scheme where either both of the faces are clockwise or counter-clockwise
    the vertices on the second faces are visited in the reversed order
    */
    
    //the vertices of the common edge (v0,v1) in the second poly are visited either 
    //a) in the same order as in the first poly (different orientations)
    //b) in the reverse order (correct orientation)
    bool diff_orient = pairs[1].second == (pairs[0].second+1)%verts1.size();     //check same order
    
    if( !diff_orient )
    {
        for( std::size_t i = (pairs[0].second+1)%verts1.size(); i != pairs[1].second; i = (i+1) % verts1.size() )
            diff_verts.push_back( verts1[i] );
    }
    else
    {
        for( std::size_t i = (pairs[1].second+1)%verts1.size(); i != pairs[0].second; i = (i+1) % verts1.size() )
            diff_verts.push_front( verts1[i] );     //reverse the order of traversal with the stack
    }

    //insert into the first polygons vertices the seconds one
    verts0.insert( verts0.begin() + pairs[0].first + 1, diff_verts.begin(), diff_verts.end() );

    //build a new face with the data of f0
    //face_ptr_t nf( new face_t( *f0 ) );
    //nf->assign( verts0.begin(), verts0.end() );
    face_ptr_t nf( new face_t( verts0.begin(), verts0.end() ) );
    add_face( nf );

    remove_face( f0 );
    remove_face( f1 );

    return nf;
}

template<typename Traits>
bool mesh<Traits>::flip_face( face_ptr_t f )
{
    //just flip the orientation of edges
    f->EdgesOrientation.flip();

    return true;
}

template<typename Traits>
void mesh<Traits>::set_point( vertex_ptr_t v, const point_t &point )
{
    v->Point = point;
}
