#include <vector>

#include <glm/glm.hpp>
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

    void Engine::computeElasticForces(Vvec3 &pos, Vvec3 &f_elas)
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

    void Engine::computeForceDifferentials(Vvec3 &pos, Vvec3 &dpos, Vvec3 &df_elas)
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

    void solveStaticPos()
    {
        /* initialize v_pos_next */
        for(size_t i = 0; i < m_size; i++)
            v_pos_next[i] = v_pos[i];

        /* initialize temp variable for iterative implicit solving */
        Vvec3 f_elas;
        computeElasticForces(v_pos_next, f_elas); 
        for(size_t i = 0; i < m_size; i++)
            f_elas[i] += f_ext[i];
        
        Vvec3 dv_pos( m_size, Vec3(0) );

        /*  conjugate gradient begin
         *  initialize temp variable
         */
        Vvec3 d( m_size, Vec3(0) );
        Vvec3 r( m_size, Vec3(0) );
        Vvec3 Ad( m_size, Vec3(0) );
        Vvec3 &x = dv_pos;
        Vvec3 &b = f_elas;

        /* while not converge f == 0, iterate */
        while ( vvec3Dot(f_elas, f_elas) > 1e-5 )
        {
            /* conjugate gradient method 
             * solves K dv_pos = f_elas */

            /* d0 = r0 = b - Ax0*/
            computeForceDifferentials(v_pos_next, dv_pos, r);
            for(size_t i = 0; i < m_size; i++)
            {
                r[i] = b[i] - r[i];
                d[i] = r[i];
            }

            float riTri = vvec3Dot(r, r);  /* compute dot( r_i, r_i )*/
            while (riTri > 1e-5)
            {
                computeForceDifferentials(v_pos_next, d, Ad); /* compute A * d_i */
                float alpha = riTri / vvec3Dot(d, Ad);

                for(size_t i = 0; i < m_size; i++)
                {
                    x[i] += alpha * d[i];       /* update x */
                    r[i] -= alpha * Ad[i];      /* update r */
                }

                float riTri_next = vvec3Dot(r, r);  /* compute dot( r_i+1, r_i+1 )*/
                float beta = riTri_next / riTri; /* compute beta */
                riTri = riTri_next;             /* update riTri */

                for(size_t i = 0; i < m_size; i++)
                    d[i] += r[i] + beta * d[i]; /* update d */
            }
            /* conjugate gradient end */

            /* update v_pos_next and f_elas */
            for(size_t i = 0; i < m_size; i++)
                v_pos_next[i] += dv_pos[i];

            /* update f_elas */
            computeElasticForces(v_pos_next, f_elas); 
            for(size_t i = 0; i < m_size; i++)
                f_elas[i] += f_ext[i];
        }

    }

    void Engine::forceTest()
    {
        Vvec3 f_elas;
        f_elas.assign( m_size, Vec3(0) );

        /* compute nodal force for each vertex */
        computeElasticForces(v_pos, f_elas);

        /* coutput data */
        for (size_t i = 0; i < m_size; i++)
        {
            m_tetra->vertices[i].m_pos = v_pos[i];
            m_tetra->vertices[i].m_velocity = f_elas[i];
        }
    }

}
