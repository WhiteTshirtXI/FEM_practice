#include <vector>

#include <glm/glm.hpp>

#include "Types.h"
#include "TetraMesh.h"
#include "Dynamic.h"

namespace BallonFEM
{
    //////////////////////////////////////////////////////////////////////
    /* ObjState implementation */
    ObjState::ObjState(TetraMesh* tetra)
    {
       this->input(tetra);
    }

    ObjState::ObjState(const ObjState& other)
    {
		*this = other;
    }

    ObjState& ObjState::operator=(const ObjState& other)
    {
		m_tetra = other.m_tetra;
		m_size = other.m_size;

		world_space_pos.assign(other.world_space_pos.begin(), other.world_space_pos.end());
		m_pos.assign(other.m_pos.begin(), other.m_pos.end());

		return *this;
    }

    void ObjState::input(TetraMesh* tetra)
    {
        m_tetra = tetra;
        m_size = tetra->vertices.size();

        world_space_pos.assign( m_size, Vec3(0));

        for (size_t i = 0; i < m_size; i++)
        {
           world_space_pos[i] = m_tetra->vertices[i].m_pos; 
        }

        m_pos.assign( world_space_pos.begin(), world_space_pos.end() );
    }

    void ObjState::output()
    {
        this->project();

        for (size_t i = 0; i < m_size; i++)
        {
           m_tetra->vertices[i].m_pos = world_space_pos[i]; 
        }
    }

    void ObjState::project()
    {
        world_space_pos.assign(m_pos.begin(), m_pos.end());
    }


    void ObjState::update(DeltaState& dpos)
    {
        for (size_t i = 0; i < m_size; i++)
        {
            this->m_pos[i] += dpos.m_pos[i];
        }
    }

    //////////////////////////////////////////////////////////////////////////
    /* DeltaStae implementation */
    DeltaState::DeltaState(const ObjState& other)
    {
        m_tetra = other.m_tetra;
        m_size = other.m_size;

        m_pos.assign( m_size, Vec3(0));
		world_space_pos.assign(m_size, Vec3(0));
    }

    void DeltaState::assign(const double alpha, const DeltaState& other)
    {
        for (size_t i = 0; i < m_size; i++)
            m_pos[i] = alpha * other.m_pos[i];
    }

    void DeltaState::addup( const double alpha, const DeltaState& other)
    {
        for (size_t i = 0; i < m_size; i++)
            m_pos[i] += alpha * other.m_pos[i];
    }

    void DeltaState::multiply(const double alpha)
    {
        for (size_t i = 0; i < m_size; i++)
            m_pos[i] *= alpha;
    }

    void DeltaState::project()
    {
        for (size_t i = 0; i < m_size; i++)
            world_space_pos[i] = m_pos[i]; 
    }

    void DeltaState::conterProject()
    {
        for (size_t i = 0; i < m_size; i++)
            m_pos[i] = world_space_pos[i];

        std::vector<size_t>::iterator i;
        for(i = m_tetra->fixed_ids.begin(); i != m_tetra->fixed_ids.end(); i++)
        {
            m_pos[*i] = Vec3(0);
        }
    }

    double DeltaState::dot(const DeltaState& other)
    {
        double count = 0;

        for (size_t i = 0; i < m_size; i++)
			count += glm::dot(m_pos[i], other.m_pos[i]);
		
		return count;
    }

}
