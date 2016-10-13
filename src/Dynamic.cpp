#include <vector>

#include <glm/glm.hpp>

#include "Types.h"
#include "Dynamic.h"

namespace BallonFEM
{
    Engine::Engine(TetraMesh* tetra, ElasticModel* model)
    {
        m_tetra = tetra;
        m_size = tetra->vertices.size();
        m_model = model;

        /* initialize */
        v_pos.assign        ( m_size, Vec3(0));
        v_velocity.assign   ( m_size, Vec3(0));
        f_ext.assign        ( m_size, Vec3(0));
        v_pos_next.assign   ( m_size, Vec3(0));
        v_velo_next.assign  ( m_size, Vec3(0));

        /* load data */
        for (size_t i = 0; i < m_size; i++)
        {
            Vertex &v = tetra->vertices[i];

            v_pos[i] = v.m_pos;

            v_velocity[i] = v.m_velocity;

            f_ext[i] = v.m_f_ext; 
        }
    }

    void Engine::labelFixedId()
    {
        fixed_id.clear();
        for (size_t i = 0; i < m_size; i++)
        {
            if (m_tetra->vertices[i].m_fixed)
                fixed_id.push_back(i);
        }
    }

	void Engine::inputData()
	{
		for (size_t i = 0; i < m_size; i++)
		{
			Vertex &v = m_tetra->vertices[i];
			v_pos[i] = v.m_pos;
			v_velocity[i] = v.m_velocity;
			f_ext[i] = v.m_f_ext;
		}
	}

    void Engine::outputData()
    {
        for (size_t i = 0; i < m_size; i++)
        {
			Vertex &v = m_tetra->vertices[i];
            v.m_pos = v_pos[i];
            v.m_velocity = v_velocity[i];
        }
    }

    void Engine::stepToNext()
    {
        for (size_t i = 0; i < m_size; i++)
        {
            v_pos[i] = v_pos_next[i];
            v_velocity[i] = v_velo_next[i];
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

    void Engine::solveNextTimestep(double timestep)
    {
    }


	#define CONVERGE_ERROR_RATE 0.01
    void Engine::solveStaticPos()
    {
        /* initialize v_pos_next */
        for(size_t i = 0; i < m_size; i++)
            v_pos_next[i] = v_pos[i];

        /* initialize temp variable for iterative implicit solving */
        Vvec3 f_elas;
        computeElasticForces(v_pos_next, f_elas); 
        for(size_t i = 0; i < m_size; i++)
            f_elas[i] += f_ext[i];
		purifyFixedId(f_elas);

        Vvec3 dv_pos( m_size, Vec3(0) );

        /*  conjugate gradient begin
         *  initialize temp variable
         */
        Vvec3 d( m_size, Vec3(0) );
        Vvec3 r( m_size, Vec3(0) );
        Vvec3 Ad( m_size, Vec3(0) );
        Vvec3 &x = dv_pos;
        Vvec3 &b = f_elas;
		
		double err_felas = vvec3Dot(f_elas, f_elas);
		double err_begin = err_felas;
		int count_iter = 0;		/* K dx = f iter count */
        /* while not converge f == 0, iterate */
		while ((err_felas > CONVERGE_ERROR_RATE * err_begin) && (err_felas > 1e-5))
        {
            /* debug use */
			count_iter++;
            printf("%d iter of K dv = f , err_felas = %f \n", count_iter, err_felas);

            /* conjugate gradient method 
             * solves K dv_pos = f_elas 
			 * since here A == -K, the
			 * eq becomes A dx = - f_elas
			 */

            /* d0 = r0 = b - Ax0 */
            computeForceDifferentials(v_pos_next, dv_pos, Ad);
            for(size_t i = 0; i < m_size; i++)
                r[i] = - b[i] - Ad[i];
			purifyFixedId(r);

			for (size_t i = 0; i < m_size; i++)
                d[i] = r[i];

			int count_iter_cg = 0;		   /* conjugate gradient iter count */
            double riTri = vvec3Dot(r, r);  /* compute dot( r_i, r_i )*/
			double riTri_begin = riTri;
			while ( (riTri > CONVERGE_ERROR_RATE * riTri_begin ) && ( count_iter_cg < 100 ) )
            {
                /* debug use */
				count_iter_cg++;
                printf("%d iter: %d iteration of CG , riTri = %f \n", count_iter, count_iter_cg, riTri);

                computeForceDifferentials(v_pos_next, d, Ad); /* compute A * d_i */
				purifyFixedId(Ad);
                double alpha = riTri / vvec3Dot(d, Ad);

				/* update x */
                for(size_t i = 0; i < m_size; i++)
                    x[i] += alpha * d[i];       

				/* update r */
				if (count_iter_cg % 50 == 0){
					/* re-calculate r */
					computeForceDifferentials(v_pos_next, dv_pos, Ad);
					for (size_t i = 0; i < m_size; i++)
						r[i] = - b[i] - Ad[i];
					purifyFixedId(r);
				}
				else{
					for (size_t i = 0; i < m_size; i++)
						r[i] -= alpha * Ad[i];      
				}

                double riTri_next = vvec3Dot(r, r);  /* compute dot( r_i+1, r_i+1 )*/
                double beta = riTri_next / riTri; /* compute beta */
                riTri = riTri_next;             /* update riTri */

                for(size_t i = 0; i < m_size; i++)
                    d[i] = r[i] + beta * d[i];  /* update d */
            }
            /* conjugate gradient end */

            /* update v_pos_next and f_elas */
            for(size_t i = 0; i < m_size; i++)
                v_pos_next[i] += dv_pos[i];

            /* update f_elas */
            computeElasticForces(v_pos_next, f_elas); 
            for(size_t i = 0; i < m_size; i++)
                f_elas[i] += f_ext[i];
			purifyFixedId(f_elas);

			err_felas = vvec3Dot(f_elas, f_elas);
        }
		/* debug use */
		printf("f_sum error %f \n", err_felas);
		printf("finish solving \n");
    }

    void Engine::forceTest()
    {
        Vvec3 f_elas;
        f_elas.assign( m_size, Vec3(0) );

        /* compute nodal force for each vertex */
        computeElasticForces(v_pos, f_elas);

        /* output data */
        for (size_t i = 0; i < m_size; i++)
        {
            m_tetra->vertices[i].m_pos = v_pos[i];
            m_tetra->vertices[i].m_velocity = f_elas[i];
        }
    }

}
