#include <vector>

#include <glm/glm.hpp>

#include <Eigen/Sparse>

#include "Types.h"
#include "Viewer.h"
#include "Dynamic.h"

extern View::Viewer *p_viewer;

namespace BalloonFEM
{
    Engine::Engine(TetraMesh* tetra, 
            ElasticModel*   volume_model, 
            AirModel*       air_model, 
            FilmModel*      film_model,
            BendingModel*   bend_model)
    {
        m_tetra = tetra;
        m_size = tetra->vertices.size();
        m_volume_model = volume_model;
        m_air_model = air_model;
		m_film_model = film_model;
        m_bend_model = bend_model;

		f_ext.assign(m_size, Vec3(0));

		cur_state = new ObjState();
		next_state = new ObjState();

		this->inputData();
    }

	void Engine::inputData()
	{
		cur_state->input(m_tetra);

		/* load data */
		for (size_t i = 0; i < m_size; i++)
		{
			Vertex &v = m_tetra->vertices[i];
			f_ext[i] = v.m_f_ext;
		}
	}

    void Engine::outputData()
    {
        cur_state->output();
    }

    void Engine::stepToNext()
    {
        std::swap(cur_state, next_state);
    }
    
    void Engine::computeForceAndGradient(ObjState &state, SpVec &f, SpMat &A)
    {
		Vvec3 f_sum;
		f_sum.assign(m_size, Vec3(0.0));

		SpMat Tri(3 * m_tetra->num_vertex, m_tetra->num_pieces);
		SpMat Sigma(3 * m_tetra->num_vertex, 2 * m_tetra->num_pieces);

		computeElasticForces(state, f_sum);      /* compute elastic forces by tetrahedrons */
		computeFilmForces(state, f_sum);	     /* compute film forces by pieces */
		computeAirForces(state, f_sum);          /* compute forces by air pressure */

		for (size_t i = 0; i < m_size; i++) f_sum[i] += f_ext[i];   /* add external force */

		SpMat K = computeAirDiffMat(state);     /* compute forces diff by air pressure */
		K += computeFilmDiffMat(state);         /* compute film forces by pieces */
		K += computeElasticDiffMat(state);      /* compute elastic forces diff by tetrahedrons */
		//K -= bendingForceAndGradient(state, f_sum); /* compute bending force and gradient */
        
        /* convert force to SpVec */
        SpVec f_real = SpVec::Zero(3 * m_size);
        for (size_t i = 0; i < f_sum.size(); i++)
        {
            f_real( 3 * i ) = f_sum[i].x;
            f_real( 3 * i + 1) = f_sum[i].y;
            f_real( 3 * i + 2) = f_sum[i].z;
        }

        /* convert K to \tilde K with W transfer. The restricted vertices
         * has all 0 colume and raw so we add 1 to its diagnal */
        A = - state.projectMat().transpose() * K * state.projectMat() + state.restrictedMat().transpose() * state.restrictedMat();

        f = state.projectMat().transpose() * f_real;
    }

    #define CONVERGE_ERROR_RATE 1e-4
    void Engine::solveStaticPos()
    {
		/* initialize next_state */
        *next_state = *cur_state;
		next_state->project();
        SpVec f_sum = SpVec::Zero(next_state->freedomDegree());

        /* initialize temp variable for iterative implicit solving */
        /* f is total force on each vertex */
        /* K = - df/dr, here Force Diff Mat compute df/dr */
        SpMat K;
        computeForceAndGradient(*next_state, f_sum, K);
        
        /* solver */
        Eigen::SimplicialLDLT<SpMat> solver;
        SpVec &b = f_sum;
		
		double err_f = f_sum.dot(f_sum);
		double err_begin = err_f;
		int count_iter = 0;		/* K dx = f iter count */
        /* while not converge f == 0, iterate */
		while ((err_f > CONVERGE_ERROR_RATE * err_begin) && (err_f > 1e-10))
        {
            /* debug use */
			count_iter++;
            printf("%d iter of K dv = f , err_felas = %.4e \n", count_iter, err_f);

            printf("building solver\n");
            solver.compute(K);
            if (solver.info() != Eigen::Success)
            {
                printf("decomposition failed!\n");
				printf("Number of non zeros: %d \n", K.nonZeros());
                return;
            }
			printf("solve delta_x \n");
            SpVec dstate = solver.solve(b);

            /* update v_pos_next and f_sum */
            next_state->update(dstate);

			/* debug watch use*/
			next_state->output();
			p_viewer->refresh(1);
			//Control::mOutput();
			
            /* update K and f*/
            computeForceAndGradient(*next_state, f_sum, K);

			/* control step length */
			dstate = - dstate;
			int cut_count = 0;
			double err_next = f_sum.dot(f_sum);
			while (err_f < err_next && cut_count < 5)
			{
				printf("cutting down dstate by half.\n");
				dstate /= 2.0;
				next_state->update( dstate);
				computeForceAndGradient(*next_state, f_sum, K);
				err_next = f_sum.dot(f_sum);
				cut_count++;
			}
			err_f = err_next;
        }

		printf("f_sum error %f \n", err_f);
		printf("finish solving \n");
    }


  //  void Engine::forceTest(Vvec3 &f_sum, Vvec3 &f_loc)
  //  {
  //      f_sum.assign( m_size, Vec3(0) );
		//f_loc = cur_state.world_space_pos;

		//SpVec f = SpVec::Zero(cur_state.freedomDegree());
  //      /* compute nodal force for each vertex */
  //      return computeForceAndGradient(cur_state, f_sum);
  //  }

}
