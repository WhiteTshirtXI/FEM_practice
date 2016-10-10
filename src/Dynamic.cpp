#include <vector>

#include <glm/glm.hpp>
#include "ElasticModel.h"
#include "Dynamic.h"

namespace BallonFEM
{
    Engine::Engine(TetraMesh* tetra, ElasticModel* model)
    {
        m_tetra = tetra;
        m_size = tetra->vertices.size();
        m_model = model;

        /* initialize */
        f_elas.assign       ( m_size, Vec3(0));
        df_elas.assign      ( m_size, Vec3(0));
        v_pos.assign        ( m_size, Vec3(0));
        v_velocity.assign   ( m_size, Vec3(0));
        v_pos_next.assign   ( m_size, Vec3(0));
        v_velo_next.assign  ( m_size, Vec3(0));
        dv_pos_next.assign  ( m_size, Vec3(0));

        /* load data */
        for (size_t i = 0; i < m_size; i++)
        {
            Vec3 &pos = tetra->vertices[i].m_pos;
            v_pos[i] = pos;
            v_pos_next[i] = pos;

            Vec3 &velocity = tetra->vertices[i].m_velocity;
            v_velocity[i] = velocity;
            v_velo_next[i] = velocity;
        }
    }

    void Engine::outputData()
    {
        for (size_t i = 0; i < m_size; i++)
        {
            m_tetra->vertices[i].m_pos = v_pos[i];
            m_tetra->vertices[i].m_velocity = v_velocity[i];
        }
    }

    void Engine::computeElasticForces(Vvec3 &pos)
    {
        f_elas.assign( m_size, Vec3(0.0));

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
    }

    void Engine::computeForceDifferentials(Vvec3 &pos, Vvec3 &dpos)
    {
        df_elas.assign( m_size, Vec3(0.0));

        for(TIter t = m_tetra->tetrahedrons.begin();
                t != m_tetra->tetrahedrons.end(); t++)
        {
            /* assgin world space position */
            iVec4 &id = t->v_id;
            Vec3 &v0 = pos[id[0]];
            Vec3 &v1 = pos[id[1]];
            Vec3 &v2 = pos[id[2]];
            Vec3 &v3 = pos[id[3]];
            
            /* calculate deformation in world space */
            Mat3 Ds = Mat3(v0 - v3, v1 - v3, v2 - v3);

            /* calculate deformation gradient */
            Mat3 F = Ds * t->Bm;
            
            /* assign world space delta position */
            Vec3 &dv0 = dpos[id[0]];
            Vec3 &dv1 = dpos[id[1]];
            Vec3 &dv2 = dpos[id[2]];
            Vec3 &dv3 = dpos[id[3]];
            
            /* calculate delta deformation in world space */
            Mat3 dDs = Mat3(dv0 - dv3, dv1 - dv3, dv2 - dv3);

            /* calculate delta deformation gradient */
            Mat3 dF = dDs * t->Bm;

            /* calculate delta Piola */
            Mat3 dP = m_model->StressDiff(F, dF);
 
            /* calculate forces contributed from this tetra */
            Mat3 dH = - t->W * dP * transpose(t->Bm);

            df_elas[id[0]] += dH[0];
            df_elas[id[1]] += dH[1];
            df_elas[id[2]] += dH[2];
            df_elas[id[3]] -= dH[0] + dH[1] + dH[2];     
        }
    }

    void Engine::solveNextTimestep(float timestep)
    {
    }

    void Engine::forceTest()
    {
        /* manually deformation along x axis */
        for (size_t i = 0; i < m_size; i++)
        {
            v_pos[i].x += v_pos[i].z * 0.01f;
        }

        /* compute nodal force for each vertex */
        computeElasticForces(v_pos);

        /* coutput data */
        for (size_t i = 0; i < m_size; i++)
        {
            m_tetra->vertices[i].m_pos = v_pos[i];
            m_tetra->vertices[i].m_velocity = f_elas[i];
        }

        m_tetra->recomputeSurfaceNorm();
    }

}
