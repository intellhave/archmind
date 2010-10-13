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

#ifndef GEOMETRY_ITERATORS_H
#define GEOMETRY_ITERATORS_H

#include <vector>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>
#include <iterator>

namespace geometry
{

template<typename Value, typename IterType>
class face_point_iterator : public boost::iterator_facade<
	face_point_iterator<Value,IterType>, 
	Value, 
	boost::forward_traversal_tag,
    Value                           >
{
public:
    face_point_iterator() {}

    face_point_iterator(IterType iter) : m_Iter(iter)
    {
    }
    
    template<typename OtherValue,typename OtherIterType> 
    face_point_iterator(const face_point_iterator<OtherValue,OtherIterType> &other) :
        m_Iter(other.m_Iter)
    {}

private:
    friend class boost::iterator_core_access;
    template<typename,typename> friend class face_point_iterator;

    void increment() 
	{
		m_Iter++; 
	}

    template<typename OtherValue,typename OtherIterType>
    bool equal(const face_point_iterator<OtherValue,OtherIterType> &other)const
    {
        return m_Iter == other.m_Iter;
    }

	Value dereference() const
    {
        return (*m_Iter)->point();
    }

    IterType m_Iter;
};

template<typename Value, typename IterType, typename BitSet >
class face_vertex_iterator : public boost::iterator_facade<
	face_vertex_iterator<Value,IterType,BitSet>, 
	Value, 
	boost::forward_traversal_tag>
{
public:
    face_vertex_iterator() {}

    face_vertex_iterator(IterType iter, BitSet orientation) : m_Iter(iter), m_Pos(0), m_Orientation(orientation)
    {
    }
    
    template<typename OtherValue,typename OtherIterType> 
    face_vertex_iterator(const face_vertex_iterator<OtherValue,OtherIterType,BitSet> &other) :
        m_Iter(other.m_Iter), m_Pos(other.m_Pos), m_Orientation(other.m_Orientation)
    {}

private:
    friend class boost::iterator_core_access;
    template<typename,typename,typename> friend class face_vertex_iterator;

    void increment() 
	{
		m_Pos++;
		m_Iter++; 
	}

    template<typename OtherValue,typename OtherIterType>
    bool equal(const face_vertex_iterator<OtherValue,OtherIterType,BitSet> &other)const
    {
        return m_Iter == other.m_Iter;
    }

	Value &dereference() const
    {
		return m_Orientation[m_Pos] ? 
			*((*m_Iter)->vertices_begin()) : *((*m_Iter)->vertices_begin()+1);
    }

    IterType m_Iter;
	std::size_t m_Pos;
	BitSet m_Orientation;
};

template<typename Value, typename FaceIter, typename EdgeIter>
class face_face_iterator : public boost::iterator_facade<
	face_face_iterator<Value,FaceIter,EdgeIter>, 
	Value, 
	boost::forward_traversal_tag>
{
public:
    face_face_iterator() {}

    face_face_iterator(EdgeIter iter, EdgeIter end, uid_t id) : 
		m_EdgeIter(iter), m_EdgeIterEnd(end), m_VID(id)
    {	
		//printf("face_face_iterator()\n");
		while( m_EdgeIter != m_EdgeIterEnd )
        {
			//printf("While\n");
			m_FaceIter = (*m_EdgeIter)->faces_begin();
			m_FaceIterEnd = (*m_EdgeIter)->faces_end();       

			while( m_FaceIter != m_FaceIterEnd )
			{
				//Check if this is not the original
				if( (*m_FaceIter)->unique_id() != m_VID)
				{
					//printf("found valid face\n");
					return;
				}
				else
				{
					//printf("++f\n");
					++m_FaceIter;	//go to the next
				}
			}

			//printf("++e\n");
			++m_EdgeIter;		//go to the next edge
		}
    }
    
    template<typename OtherValue,typename OtherIterType> 
    face_face_iterator(const face_face_iterator<Value,OtherValue,OtherIterType> &other) :
        m_EdgeIter(other.m_EdgeIter), 
		m_EdgeIterEnd(other.m_EdgeIterEnd), 
		m_FaceIter(other.m_FaceIter),
		m_FaceIterEnd(other.m_FaceIter),
		m_VID(other.m_VID)
    {}

private:
    friend class boost::iterator_core_access;
    template<typename, typename,typename> friend class face_face_iterator;

    void increment() 
	{
		do
        {
            //go to the next face
            ++m_FaceIter;
			
            //check if we must move to the next edge
            if( m_FaceIter == m_FaceIterEnd ) 
            {      
                //go to the next edge
                ++m_EdgeIter;

                //start traversing the faces of the edge
                if( m_EdgeIter != m_EdgeIterEnd ) 
                {
					m_FaceIter = (*m_EdgeIter)->faces_begin();
					m_FaceIterEnd = (*m_EdgeIter)->faces_end();            
                }
            }
			else
			{
				//Check if this is not the original
				if( (*m_FaceIter)->unique_id() != m_VID)
					break;
			}
		}
		while( m_EdgeIter != m_EdgeIterEnd );
	}

    template<typename OtherValue,typename OtherIterType>
    bool equal(const face_face_iterator<Value,OtherValue,OtherIterType> &other)const
    {
        return m_EdgeIter == other.m_EdgeIter;
    }

	Value &dereference() const
    {
		return *m_FaceIter;
    }

    EdgeIter m_EdgeIter;
	EdgeIter m_EdgeIterEnd;
	FaceIter m_FaceIter;
	FaceIter m_FaceIterEnd;

	uid_t m_VID;
};

template<typename Value, typename IterType>
class vertex_vertex_iterator : public boost::iterator_facade<
    vertex_vertex_iterator<Value,IterType>,
    Value,
	boost::forward_traversal_tag>
{
public:
    vertex_vertex_iterator() {}

    vertex_vertex_iterator(IterType iter, uid_t id ) : m_Iter(iter), m_VID(id)
    {}
    
    template<typename OtherValue,typename OtherIterValue> 
    vertex_vertex_iterator(const vertex_vertex_iterator<OtherValue,OtherIterValue> &other) :
        m_Iter(other.m_Iter), m_VID(other.m_VID)
    {}

private:
    friend class boost::iterator_core_access;
    template<typename,typename> friend class vertex_vertex_iterator;

    void increment() 
	{
		++m_Iter;
    }

    template<typename OtherValue,typename OtherIterValue>
    bool equal(const vertex_vertex_iterator<OtherValue,OtherIterValue> &other)const
    {
        return m_Iter == other.m_Iter;
    }

	Value &dereference() const
    {
        return
			(*((*m_Iter)->vertices_begin()))->unique_id() == m_VID ? 
			*((*m_Iter)->vertices_begin()+1) : *((*m_Iter)->vertices_begin());
    }

    IterType m_Iter;
	uid_t m_VID;
};

template<typename Value, typename EdgeIter, typename FaceIter>
class vertex_face_iterator : public boost::iterator_facade<
    vertex_face_iterator<Value,EdgeIter,FaceIter>, 
	Value,
    boost::forward_traversal_tag>
{
public:
    vertex_face_iterator() {}

    vertex_face_iterator(EdgeIter iter, EdgeIter iterend) : m_EdgeIter(iter), m_EdgeIterEnd(iterend)
    {
        //Find a valid face or reach the end
        while( m_EdgeIter != m_EdgeIterEnd )
        {
            m_FaceIter = (*m_EdgeIter)->faces_begin();
			m_FaceIterEnd = (*m_EdgeIter)->faces_end();
            
			//Check if can find a valid face
            if( m_FaceIter == m_FaceIterEnd )
                m_EdgeIter++;
            else
			{
				//valid face found
				m_Visited[ (*m_FaceIter)->unique_id() ] = true;		//mark as visited 
                break;
			}
        }
    }

private:
    friend class boost::iterator_core_access;

    void increment() 
	{
        do
        {
            //go to the next face
            ++m_FaceIter;
			
            //check if we must move to the next edge
            if( m_FaceIter == m_FaceIterEnd ) 
            {      
                //go to the next edge
                ++m_EdgeIter;

                //start traversing the faces of the edge
                if( m_EdgeIter != m_EdgeIterEnd ) 
                {
					m_FaceIter = (*m_EdgeIter)->faces_begin();
					m_FaceIterEnd = (*m_EdgeIter)->faces_end();            
                }
            }
			
			if( m_FaceIter != m_FaceIterEnd )
			{
				uid_t id = (*m_FaceIter)->unique_id();

				if( m_Visited.find( id ) == m_Visited.end() )
				{
					m_Visited[ id ] = true;
					break;
				}
			}
        }
        while( m_EdgeIter != m_EdgeIterEnd );
    }

    bool equal(const vertex_face_iterator &other)const
    {
        return m_EdgeIter == other.m_EdgeIter;
    }

	Value &dereference() const
    {
        return (*m_FaceIter);
    }

	boost::unordered_map< uid_t, bool > m_Visited;

    EdgeIter m_EdgeIter;
	EdgeIter m_EdgeIterEnd;

    FaceIter m_FaceIter;
	FaceIter m_FaceIterEnd;
};

};

#endif
