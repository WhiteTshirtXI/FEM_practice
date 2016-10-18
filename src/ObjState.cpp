#include <vector>

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

#include "Types.h"
#include "TetraMesh.h"
#include "ObjState.h"

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
		m_r_size = other.m_r_size;

		world_space_pos.assign(other.world_space_pos.begin(), other.world_space_pos.end());
		m_pos.assign(other.m_pos.begin(), other.m_pos.end());
		m_r_pos.assign(other.m_r_pos.begin(), other.m_r_pos.end());
		m_r_rot.assign(other.m_r_rot.begin(), other.m_r_rot.end());

		return *this;
    }

    void ObjState::input(TetraMesh* tetra)
    {
        m_tetra = tetra;
        m_size = tetra->vertices.size();
        m_r_size = tetra->rigid_bodies.size();

        /* initialize world_space_pos */
        world_space_pos.assign( m_size, Vec3(0));
        for (size_t i = 0; i < m_size; i++)
        {
           world_space_pos[i] = m_tetra->vertices[i].m_pos; 
        }

        /* initialize m_r_pos and m_r_rot */
        m_r_pos.assign( m_r_size, Vec3(0));
        m_r_rot.assign( m_r_size, Quat(1,0,0,0));
        for (size_t i = 0; i < m_r_size; i++)
        {
            Rigid &r = m_tetra->rigid_bodies[i];
            m_r_pos[i] = r.m_pos;
            m_r_rot[i] = r.m_rot;
        }

        /* initialize m_pos */
        m_pos.assign( world_space_pos.begin(), world_space_pos.end() );
    }

    void ObjState::output()
    {
        this->project();

        for (size_t i = 0; i < m_size; i++)
        {
           m_tetra->vertices[i].m_pos = world_space_pos[i]; 
        }

        for (size_t i = 0; i < m_r_size; i++)
        {
            Rigid &r = m_tetra->rigid_bodies[i];
            r.m_pos = m_r_pos[i];
            r.m_rot = m_r_rot[i];
        }
    }

    void ObjState::project()
    {
        /* first assign those not bounded vertex */
        world_space_pos.assign(m_pos.begin(), m_pos.end());
        
        /* calculate rigid body elements world space pos */
        for (size_t i = 0; i < m_r_size; i++)
        {
            /* rotation matrix of current rigid body i */
            Mat3 R = glm::toMat3(m_r_rot[i]);

            /* current mass center position of rigid body */
            Vec3 rc = m_r_pos[i];
            
            /* current rigid body object */
            Rigid &r = m_tetra->rigid_bodies[i];
            
            /* rigid body object mass center pos in material space */
            Vec3 cordc = r.m_cord;
            
            for (size_t j = 0; j < r.elements.size(); j++)
            {
                size_t v_id = r.elements[j];
                /* r = R * (X - Xc) + r_ c*/
                world_space_pos[v_id] = rc +
                        R * (m_tetra->vertices[v_id].m_cord - cordc);
            }
        }
    }

    Quat ObjState::omegaToQuat(Vec3 delta_phi)
    {
        double half_angle = delta_phi.length() / 2.;
        double cos = std::cos(half_angle);
        double sin = std::sin(half_angle);
        Vec3 w = delta_phi;
        glm::normalize(w);

        return Quat(cos, sin * w);
    }
        
    void ObjState::update(DeltaState& dpos)
    {
        for (size_t i = 0; i < m_size; i++)
        {
            this->m_pos[i] += dpos.dm_pos[i];
        }

        for (size_t i = 0; i < m_r_size; i++)
        {
            this->m_r_pos[i] += dpos.dm_r_pos[i];
            this->m_r_rot[i] = omegaToQuat(dpos.dm_r_rot[i]) * this->m_r_rot[i]; 
        }
    }

    //////////////////////////////////////////////////////////////////////////
    /* DeltaStae implementation */
    DeltaState::DeltaState(ObjState& other)
    {
        m_tetra = other.m_tetra;
        m_state = &other;
        m_size = other.m_size;
        m_r_size = other.m_r_size;

		world_space_pos.assign(m_size, Vec3(0));

        dm_pos.assign( m_size, Vec3(0));
        dm_r_pos.assign( m_r_size, Vec3(0));
        dm_r_rot.assign( m_r_size, Vec3(0));
    }

    void DeltaState::assign(const double alpha, const DeltaState& other)
    {
        for (size_t i = 0; i < m_size; i++)
            dm_pos[i] = alpha * other.dm_pos[i];

        for (size_t i = 0; i < m_r_size; i++)
        {
            dm_r_pos[i] = alpha * other.dm_r_pos[i];
            dm_r_rot[i] = alpha * other.dm_r_rot[i];
        }
    }

    void DeltaState::addup( const double alpha, const DeltaState& other)
    {
        for (size_t i = 0; i < m_size; i++)
            dm_pos[i] += alpha * other.dm_pos[i];

        for (size_t i = 0; i < m_r_size; i++)
        {
            dm_r_pos[i] += alpha * other.dm_r_pos[i];
            dm_r_rot[i] += alpha * other.dm_r_rot[i];
        }
    }

    void DeltaState::multiply(const double alpha)
    {
         for (size_t i = 0; i < m_size; i++)
            dm_pos[i] *= alpha;

        for (size_t i = 0; i < m_r_size; i++)
        {
            dm_r_pos[i] *= alpha; 
            dm_r_rot[i] *= alpha; 
        }
    }

    void DeltaState::project()
    {
       /* first assign those not bounded vertex */
        world_space_pos.assign(dm_pos.begin(), dm_pos.end());
        
        /* calculate rigid body elements world space pos */
        for (size_t i = 0; i < m_r_size; i++)
        {
            /* rotation matrix of current rigid body i */
            Vec3 dR = dm_r_rot[i];

            /* current mass center position of rigid body */
            Vec3 drc = dm_r_pos[i];
            Vec3 rc = m_state->m_r_pos[i];
            
            /* current rigid body object */
            Rigid &r = m_tetra->rigid_bodies[i];
            
            for (size_t j = 0; j < r.elements.size(); j++)
            {
                size_t v_id = r.elements[j];
                /* dr = drc + dR x (r - rc) */
                world_space_pos[v_id] = drc + 
                    glm::cross(dR, m_state->m_pos[v_id] - rc);
            }
        }
		
		/* for those fixed point world_spce_pos shold be set to 0 */
		std::vector<size_t>::iterator i;
		for (i = m_tetra->fixed_ids.begin(); i != m_tetra->fixed_ids.end(); i++)
		{
			world_space_pos[*i] = Vec3(0);
		}
    }

    void DeltaState::conterProject()
    {
        /* copy data of free verteics from world_space_pos*/
        dm_pos.assign(world_space_pos.begin(), world_space_pos.end());

        /* make fixed vertex delta to be 0 */
        std::vector<size_t>::iterator i;
        for(i = m_tetra->fixed_ids.begin(); i != m_tetra->fixed_ids.end(); i++)
        {
            dm_pos[*i] = Vec3(0);
        }

        /* calculate rigid body delta */
        for (size_t i = 0; i < m_r_size; i++)
        {
            /* rotation matrix of current rigid body i */
            Vec3 &drot = dm_r_rot[i];
            drot = Vec3(0);

            /* current mass center position of rigid body */
            Vec3 &drc = dm_r_pos[i];
            drc = Vec3(0);

            /* current rigid body mass center pos */
            Vec3 rc = m_state->m_r_pos[i];

            /* current rigid body object */
            Rigid &r = m_tetra->rigid_bodies[i];
            
            for (size_t j = 0; j < r.elements.size(); j++)
            {
                size_t v_id = r.elements[j];
				
				/* clear dm_pos */
				dm_pos[v_id] = Vec3(0);

                /* drc = sum(dr) */
                drc += world_space_pos[v_id];

                /* drot = sum(r x dr) */
                drot += glm::cross(m_state->m_pos[v_id] - rc, 
                                   world_space_pos[v_id]);
            }
        } 
    }

    double DeltaState::dot(const DeltaState& other)
    {
        double count = 0;

        for (size_t i = 0; i < m_size; i++)
			count += glm::dot(dm_pos[i], other.dm_pos[i]);

        for (size_t i = 0; i <  m_r_size; i++)
        {
            count += glm::dot(dm_r_pos[i], other.dm_r_pos[i]);
            count += glm::dot(dm_r_rot[i], other.dm_r_rot[i]);
        }
		
		return count;
    }

}
