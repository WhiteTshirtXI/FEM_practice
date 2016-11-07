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
    Engine::Engine(TetraMesh* tetra, ElasticModel* volume_model, AirModel* air_model, FilmModel* film_model)
    {
        m_tetra = tetra;
        m_size = tetra->vertices.size();
        m_volume_model = volume_model;
        m_air_model = air_model;
		m_film_model = film_model;

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
    
    SpMat Engine::computeForceAndGradient(ObjState &state, Vvec3 &f_sum)
    {
        f_sum.assign( m_size, Vec3(0.0));

        /* compute elastic forces by tetrahedrons */
        computeElasticForces(state, f_sum); 

		/* compute film forces by pieces */
		//computeFilmForces(state, f_sum);
        
        /* compute forces by air pressure */
		computeAirForces(state, f_sum);

        for(size_t i = 0; i < m_size; i++)
            f_sum[i] += f_ext[i];

		/* compute forces diff by air pressure */
        SpMat K = computeAirDiffMat(state);

		/* compute film forces by pieces */
		//K += computeFilmDiffMat(state);

		/* compute elastic forces diff by tetrahedrons */
		K += computeElasticDiffMat(state);

        /* compute bending force and gradient */
        K -= bendingForceAndGradient(state, f_sum);

        return K;
    }

    #define CONVERGE_ERROR_RATE 1e-4
    void Engine::solveStaticPos()
    {
      /* initialize next_state */
        next_state = cur_state;
		next_state.project();
        DeltaState f_sum(next_state);

        /* initialize temp variable for iterative implicit solving */
        /* f is total force on each vertex */
        /* K = - df/dr, here Force Diff Mat compute df/dr */
        SpMat K = - computeForceAndGradient(next_state, f_sum.world_space_pos);
        f_sum.conterProject();
        SpVec b = f_sum.toSpVec();


        DeltaState dstate(next_state);

        /* convert K to \tilde K with W transfer. The restricted vertices
         * has all 0 colume and raw so we add 1 to its diagnal */
        SpMat W = dstate.projectMat();
        SpMat R = dstate.restrictedMat();
        K = W.transpose() * K * W + R;
        
        /* solver */
        Eigen::SimplicialLDLT<SpMat> solver;
		
		double err_f = f_sum.dot(f_sum);
		double err_begin = err_f;
		int count_iter = 0;		/* K dx = f iter count */
        /* while not converge f == 0, iterate */
		while ((err_f > CONVERGE_ERROR_RATE * err_begin) && (err_f > 1e-10))
        {
            /* debug use */
			count_iter++;
            printf("%d iter of K dv = f , err_felas = %f \n", count_iter, err_f);

            /* r0 = b - Ax0 */
            SpVec r = b - K * dstate.toSpVec();

            printf("building solver\n");
            solver.compute(K);
            if (solver.info() != Eigen::Success)
            {
                printf("decomposition failed!\n");
				printf("Number of non zeros: %d \n", K.nonZeros());
                return;
            }
			printf("solve delta_x \n");
            SpVec x = solver.solve(r);
            dstate.readSpVec(x);

            /* update v_pos_next and f_sum */
            next_state.update(dstate);
			next_state.project();
            dstate.clear();

			/* debug watch use*/
			next_state.output();
			shadFlag = 1;
			p_viewer->refresh();
			//Control::mOutput();
			
            /* update K and f*/
            K = - computeForceAndGradient(next_state, f_sum.world_space_pos);
            K = W.transpose() * K * W + R;
            f_sum.conterProject();
            b = f_sum.toSpVec();

			err_f = f_sum.dot(f_sum);
        }

		printf("f_sum error %f \n", err_f);
		printf("finish solving \n");
    }


    void Engine::forceTest()
    {
        Vvec3 f_sum;
        f_sum.assign( m_size, Vec3(0) );

        /* compute nodal force for each vertex */
        computeElasticForces(cur_state, f_sum);

        /* output data */
        for (size_t i = 0; i < m_size; i++)
        {
            m_tetra->vertices[i].m_pos = cur_state.world_space_pos[i];
            m_tetra->vertices[i].m_velocity = f_sum[i];
        }
    }

}
