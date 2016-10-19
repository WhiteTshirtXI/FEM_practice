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
        for (size_t i = 0; i < m_size; i++)
        {
			Vertex &v = m_tetra->vertices[i];
            v.m_pos = cur_state.world_space_pos[i];
        }
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
    }

    void Engine::computeForceDifferentials(ObjState &state, DeltaState &dstate, Vvec3 &df_elas)
    {
        df_elas.assign( m_size, Vec3(0.0));
        
        /* project from constrained freedom state to world space */
        Vvec3 &pos = state.world_space_pos;

        dstate.project();
        Vvec3 &dpos = dstate.world_space_pos;

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


	#define CONVERGE_ERROR_RATE 1e-4
    void Engine::solveStaticPos()
    {
        /* initialize next_state */
        next_state = cur_state;
		next_state.project();

        /* initialize temp variable for iterative implicit solving */
        DeltaState f_elas(next_state);
        computeElasticForces(next_state, f_elas.world_space_pos); 
        for(size_t i = 0; i < m_size; i++)
            f_elas.world_space_pos[i] += f_ext[i];
        f_elas.conterProject();

        DeltaState dstate(next_state);

        /*  conjugate gradient begin
         *  initialize temp variable
         */
        DeltaState d( next_state );
        DeltaState r( next_state );
        DeltaState Ad(next_state );
        DeltaState &x = dstate;
        DeltaState &b = f_elas;
		
		double err_felas = f_elas.dot(f_elas);
		double err_begin = err_felas;
		int count_iter = 0;		/* K dx = f iter count */
        /* while not converge f == 0, iterate */
		while ((err_felas > CONVERGE_ERROR_RATE * err_begin) && (err_felas > 1e-10))
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
            computeForceDifferentials(next_state, x, Ad.world_space_pos);
            Ad.conterProject();
            r.assign(-1, b);
            r.addup(-1, Ad);

            d.assign(1, r);

			int count_iter_cg = 0;		   /* conjugate gradient iter count */
            double riTri = r.dot(r);  /* compute dot( r_i, r_i )*/
			double riTri_begin = riTri;
			while ( (riTri > CONVERGE_ERROR_RATE * riTri_begin ) && ( count_iter_cg < 100 ) )
            {
                /* debug use */
				count_iter_cg++;
                printf("%d iter: %d iteration of CG , riTri = %f \n", count_iter, count_iter_cg, riTri);

				x.project();
				for (size_t i = 0; i < x.world_space_pos.size(); i++)
					printf("v%d %.4f %.4f %.4f, ", i, x.world_space_pos[i].x, x.world_space_pos[i].y, x.world_space_pos[i].z);
				printf("\n");

                computeForceDifferentials(next_state, d, Ad.world_space_pos); /* compute A * d_i */
                Ad.conterProject();
				double gamma = d.dot(Ad);
                double alpha = riTri / gamma;

				/* update x */
                    x.addup(alpha, d);       

				/* update r */
				if (count_iter_cg % 50 == 0){
					/* re-calculate r */
                    computeForceDifferentials(next_state, x, Ad.world_space_pos);
                    Ad.conterProject();
                    r.assign(-1, b);
                    r.addup(-1, Ad);
				}
				else{
						r.addup(-alpha, Ad);      
				}

                double riTri_next = r.dot(r);  /* compute dot( r_i+1, r_i+1 )*/
                double beta = riTri_next / riTri; /* compute beta */
                riTri = riTri_next;             /* update riTri */

                d.multiply(beta);
                d.addup(1, r);  /* update d */
            }
            /* conjugate gradient end */

            /* update v_pos_next and f_elas */
            next_state.update(dstate);
			next_state.project();

            /* update f_elas */
            computeElasticForces(next_state, f_elas.world_space_pos); 
            for(size_t i = 0; i < m_size; i++)
                f_elas.world_space_pos[i] += f_ext[i];
            f_elas.conterProject();

			err_felas = f_elas.dot(f_elas);
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
        computeElasticForces(cur_state, f_elas);

        /* output data */
        for (size_t i = 0; i < m_size; i++)
        {
            m_tetra->vertices[i].m_pos = cur_state.world_space_pos[i];
            m_tetra->vertices[i].m_velocity = f_elas[i];
        }
    }

}
