#include <vector>

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

#include "Types.h"
#include "TetraMesh.h"
#include "ObjState.h"

namespace BalloonFEM
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
        m_h_size = other.m_h_size;

		world_space_pos.assign(other.world_space_pos.begin(), other.world_space_pos.end());
        volume_gradient.assign(other.volume_gradient.begin(), other.volume_gradient.end());
        hole_volume.assign(other.hole_volume.begin(), other.hole_volume.end());

		m_pos.assign(other.m_pos.begin(), other.m_pos.end());
		m_r_pos.assign(other.m_r_pos.begin(), other.m_r_pos.end());
		m_r_rot.assign(other.m_r_rot.begin(), other.m_r_rot.end());

		return *this;
    }

    void ObjState::input(TetraMesh* tetra)
    {
        m_tetra = tetra;
        m_size = tetra->vertices.size();
        m_r_size = tetra->rigids.size();
        m_h_size = tetra->holes.size();

        /* initialize m_pos */
        m_pos.assign( m_size, Vec3(0));
        for (size_t i = 0; i < m_size; i++)
        {
           m_pos[i] = m_tetra->vertices[i].m_pos; 
        }


        /* initialize m_r_pos and m_r_rot */
        m_r_pos.assign( m_r_size, Vec3(0));
        m_r_rot.assign( m_r_size, Quat(1,0,0,0));
        for (size_t i = 0; i < m_r_size; i++)
        {
            Rigid &r = m_tetra->rigids[i];
            m_r_pos[i] = r.m_pos;
            m_r_rot[i] = r.m_rot;
        }

        /* initialize world_space_pos, volume gradients and hole volume */
        world_space_pos.assign( m_pos.begin(), m_pos.end() );
        volume_gradient.assign( m_size, Vec3(0) );
        hole_volume.assign( m_h_size, 0 );
        this->project();
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
            Rigid &r = m_tetra->rigids[i];
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
            Rigid &r = m_tetra->rigids[i];
            
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

        this->holeVolume();
        this->volumeGradient();
    }

    void ObjState::holeVolume()
    {
        Vvec3 &pos = world_space_pos;

        for(size_t i = 0; i < m_h_size; i++)
        {
            double &V = hole_volume[i];
            Hole &h = m_tetra->holes[i];

            V = 0;
            for(FCIter f = h.holeface.begin(); f != h.holeface.end(); f++)
            {
                iVec3 v = f->v_id;
                V += glm::dot(glm::cross(pos[v[0]], pos[v[1]]), pos[v[2]]);
            }

        /* since the hole surface normal towards inner, the V is negtive*/
            V /= - 6.0;
        }
    }

    void ObjState::volumeGradient()
    {   
        Vvec3 &pos = world_space_pos;
        volume_gradient.assign( m_size, Vec3(0));

        for(size_t i = 0; i < m_h_size; i++)
        {
            Hole &h = m_tetra->holes[i];
            for(FCIter f = h.holeface.begin(); f != h.holeface.end(); f++)
            {
                iVec3 v = f->v_id;
                /* 2An = (r1 - r3) x (r2 - r3) */
                Vec3 n = glm::cross(pos[v[0]] - pos[v[2]], pos[v[1]] - pos[v[2]]);
                volume_gradient[v[0]] += n;
                volume_gradient[v[1]] += n;
                volume_gradient[v[2]] += n;
            }
        }

        /* since the hole surface normal towards inner, the V is negtive*/
        for(size_t i = 0; i < m_size; i++)
            volume_gradient[i] /= - 6.0;
    }

    void ObjState::volumeGradientDiff(Vvec3& dr, Vvec3& dg)
    {
        Vvec3 &pos = world_space_pos;
        dg.assign( m_size, Vec3(0));

        for(size_t i = 0; i < m_h_size; i++)
        {
            Hole &h = m_tetra->holes[i];
            for(FCIter f = h.holeface.begin(); f != h.holeface.end(); f++)
            {
                iVec3 v = f->v_id;
                /* d(2An) = (dr1 - dr3) x (r2 -r3) + (r1 - r3) x (dr2 - dr3) */
                Vec3 dn = glm::cross(dr[v[0]] - dr[v[2]], pos[v[1]] - pos[v[2]])
                        + glm::cross(pos[v[0]] - pos[v[2]], dr[v[1]] - dr[v[2]]);
                dg[v[0]] += dn;
                dg[v[1]] += dn;
                dg[v[2]] += dn;
            }
        }

        /* since the hole surface normal towards inner, the V is negtive*/
        for(size_t i = 0; i < m_size; i++)
            dg[i] /= - 6.0;
    }

    Quat ObjState::omegaToQuat(Vec3 delta_phi)
    {
        double angle = length(delta_phi) ;
        double cos = std::cos(angle / 2.0);
        double sin = std::sin(angle / 2.0 ) / angle;

        return Quat(cos, sin * delta_phi);
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

	void DeltaState::clear()
	{
		world_space_pos.assign(m_size, Vec3(0));

		dm_pos.assign(m_size, Vec3(0));
		dm_r_pos.assign(m_r_size, Vec3(0));
		dm_r_rot.assign(m_r_size, Vec3(0));
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
            Vec3 drot = dm_r_rot[i];

            /* current mass center position of rigid body */
            Vec3 drc = dm_r_pos[i];
            Vec3 rc = m_state->m_r_pos[i];
            
            /* current rigid body object */
            Rigid &r = m_tetra->rigids[i];
            
            for (size_t j = 0; j < r.elements.size(); j++)
            {
                size_t v_id = r.elements[j];
                /* dr = drc + dR x (r - rc) */
                world_space_pos[v_id] = drc + 
                    glm::cross(drot, m_state->world_space_pos[v_id] - rc);
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
            Rigid &r = m_tetra->rigids[i];
            
            for (size_t j = 0; j < r.elements.size(); j++)
            {
                size_t v_id = r.elements[j];
				
				/* clear dm_pos */
				dm_pos[v_id] = Vec3(0);

                /* drc = sum(dr) */
                drc += world_space_pos[v_id];

                /* drot = sum(r x dr) */
                drot += glm::cross(m_state->world_space_pos[v_id] - rc, 
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
