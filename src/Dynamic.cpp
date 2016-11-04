#include <vector>

#include <glm/glm.hpp>

#include <Eigen/Sparse>

#include "Types.h"
#include "Viewer.h"
#include "controls.h"
#include "Dynamic.h"

extern int shadFlag;
extern View::Viewer *p_viewer;

namespace BalloonFEM
{
    Engine::Engine(TetraMesh* tetra, ElasticModel* model, AirModel* a_model)
    {
        m_tetra = tetra;
        m_size = tetra->vertices.size();
        m_model = model;
        m_a_model = a_model;

		f_ext.assign(m_size, Vec3(0));

		this->inputData();
    }

	void Engine::inputData()
	{
        cur_state.input(m_tetra);

		/* load data */
		for (size_t i = 0; i < m_size; i++)
		{
			Vertex &v = m_tetra->vertices[i];
			f_ext[i] = v.m_f_ext;
		}
	}

    void Engine::outputData()
    {
        cur_state.project();
        cur_state.output();
    }

    void Engine::stepToNext()
    {
        std::swap(cur_state, next_state);
    }
    
    void Engine::computeElasticForces(ObjState &state, Vvec3 &f_elas)
    {
        f_elas.assign( m_size, Vec3(0.0));

        Vvec3 &pos = state.world_space_pos;

        for(TIter t = m_tetra->tetrahedrons.begin();
                t != m_tetra->tetrahedrons.end(); t++)
        {
            iVec4 &id = t->v_id;
            Vec3 &v0 = pos[id[0]];
            Vec3 &v1 = pos[id[1]];
            Vec3 &v2 = pos[id[2]];
            Vec3 &v3 = pos[id[3]];
            
            /* calculate deformation in world space */
            Mat3 Ds = Mat3(v0 - v3, v1 - v3, v2 - v3);

            /* calculate deformation gradient */
            Mat3 F = Ds * t->Bm;

            /* calculate Piola for this tetra */
            Mat3 P = m_model->Piola(F);

            /* calculate forces contributed from this tetra */
            Mat3 H = - t->W * P * transpose(t->Bm);

            f_elas[id[0]] += H[0];
            f_elas[id[1]] += H[1];
            f_elas[id[2]] += H[2];
            f_elas[id[3]] -= H[0] + H[1] + H[2];
        }

        /* add air pressure force */
        for (size_t i = 0; i < m_tetra->holes.size(); i++)
        {
            double p = m_a_model->pressure(state.hole_volume[i]);
            Hole &h = m_tetra->holes[i];
            std::vector<size_t>::iterator j;

            for (j = h.vertices.begin(); j != h.vertices.end(); j++)
            {
                f_elas[*j] += p * state.volume_gradient[*j];
            }
        }
    }

    void Engine::forceTest()
    {
        Vvec3 f_elas;
        f_elas.assign( m_size, Vec3(0) );

        /* compute nodal force for each vertex */
        computeElasticForces(cur_state, f_elas);

        /* output data */
        for (size_t i = 0; i < m_size; i++)
        {
            m_tetra->vertices[i].m_pos = cur_state.world_space_pos[i];
            m_tetra->vertices[i].m_velocity = f_elas[i];
        }
    }

}
