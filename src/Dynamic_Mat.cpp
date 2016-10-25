#include <vector>

#include <glm/glm.hpp>

#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "Types.h"
#include "Viewer.h"
#include "controls.h"
#include "Dynamic.h"

extern int shadFlag;
extern View::Viewer *p_viewer;

namespace{
    using namespace BallonFEM;
    const Vec3 v[4] = {Vec3(0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    const Mat3 m[4][3] = {
        { Mat3(v[1], v[0], v[0]), Mat3(v[2], v[0], v[0]), Mat3(v[3], v[0], v[0])},
        { Mat3(v[0], v[1], v[0]), Mat3(v[0], v[2], v[0]), Mat3(v[0], v[3], v[0])},
        { Mat3(v[0], v[0], v[1]), Mat3(v[0], v[0], v[2]), Mat3(v[0], v[0], v[3])},
        { -Mat3(v[1], v[1], v[1]), -Mat3(v[2], v[2], v[2]), -Mat3(v[3], v[3], v[3])}
    };

}

namespace BallonFEM
{

	void Engine::computeForceDiffMat(ObjState &state, SpMat &K)
	{
		printf("building force differential matrix \n");
		/* project from constrained freedom state to world space */
		Vvec3 &pos = state.world_space_pos;

		///////////////////////////////////////////////////////////////////
		/* compute air pressure force differential matrix */
		state.volumeGradientDiffMat(K);

		std::vector<T> pressure;
		pressure.reserve(pos.size());
		for (size_t i = 0; i < m_tetra->holes.size(); i++)
		{
			double p = m_a_model->pressure(state.hole_volume[i]);
			Hole &h = m_tetra->holes[i];
			std::vector<size_t>::iterator j;
            
			double dV = 0;
			for (j = h.vertices.begin(); j != h.vertices.end(); j++)
			{
				/* p(V) * dG */
				pressure.push_back( T(3 * *j    , 3 * *j    , p) );
				pressure.push_back( T(3 * *j + 1, 3 * *j + 1, p) );
				pressure.push_back( T(3 * *j + 2, 3 * *j + 2, p) );
			}
		}
		SpMat P(3 * pos.size(), 3 * pos.size());
		P.setFromTriplets(pressure.begin(), pressure.end());
		K = K * P;

		//////////////////////////////////////////////////////////////////
		/* compute elastic force Differentials */
		std::vector<T> coefficients;
		coefficients.clear();
		coefficients.reserve( 12 * 12 * pos.size());

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
            
			/* i is index of vertex, j is index of dimention */
			for (size_t i = 0; i < 4; i++)
				for(size_t j = 0; j < 3; j++)
			{
				/* calculate delta deformation in world space */
				Mat3 dDs = m[i][j]; 

				/* calculate delta deformation gradient */
				Mat3 dF = dDs * t->Bm;

				/* calculate delta Piola */
				Mat3 dP = m_model->StressDiff(F, dF);
 
				/* calculate forces contributed from this tetra */
				Mat3 dH = - t->W * dP * transpose(t->Bm);

				for(size_t w = 0; w < 3; w++)
					for(size_t l = 0; l < 3; l++)
						coefficients.push_back( T(3*id[w] + l, 3*id[i] + j,  dH[w][l]));

				Vec3 df_4 = - dH[0] - dH[1] - dH[2];
				for(size_t l = 0; l < 3; l++)
					coefficients.push_back( T(3*id[3] + l, 3*id[i] + j,  df_4[l]));
			}
		}
		SpMat E( 3 * pos.size(), 3 * pos.size());
		E.setFromTriplets(coefficients.begin(), coefficients.end());

		K += E;
    }

#define CONVERGE_ERROR_RATE 1e-4
  void Engine::solveStaticPosMat()
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
        SpVec b = f_elas.toSpVec();

        DeltaState dstate(next_state);

        /* K = - df/dr, here Force Diff Mat compute df/dr */
        SpMat K;
        computeForceDiffMat(next_state, K);
        K = -K;

        /* convert K to \tilde K with W transfer. The restricted vertices
         * has all 0 colume and raw so we add 1 to its diagnal */
        SpMat W = dstate.projectMat();
        K = W.transpose() * K * W + dstate.restrictedMat();
        
        /* solver */
        Eigen::SimplicialLDLT<SpMat> solver;
		
		double err_felas = f_elas.dot(f_elas);
		double err_begin = err_felas;
		int count_iter = 0;		/* K dx = f iter count */
        /* while not converge f == 0, iterate */
		while ((err_felas > CONVERGE_ERROR_RATE * err_begin) && (err_felas > 1e-10))
        {
            /* debug use */
			count_iter++;
            printf("%d iter of K dv = f , err_felas = %f \n", count_iter, err_felas);

            /* r0 = b - Ax0 */
            SpVec r = b - K * dstate.toSpVec();


            printf("building solver\n");
            solver.compute(K);
            if (solver.info() != Eigen::Success)
            {
                printf("decomposition failed!\n");
                return;
            }
			printf("solve delta_x \n");
            SpVec x = solver.solve(r);
            dstate.readSpVec(x);

            /* update v_pos_next and f_elas */
            next_state.update(dstate);
			next_state.project();
            dstate.clear();

			/* debug watch use*/
			next_state.output();
			shadFlag = 1;
			p_viewer->refresh();
			//Control::mOutput();
			
            /* update K */
            computeForceDiffMat(next_state, K);
            K = -K;
            dstate.projectMat();
            K = W.transpose() * K * W + dstate.restrictedMat();

            /* update f_elas */
            computeElasticForces(next_state, f_elas.world_space_pos); 
            for(size_t i = 0; i < m_size; i++)
                f_elas.world_space_pos[i] += f_ext[i];
            f_elas.conterProject();
            b = f_elas.toSpVec();

			err_felas = f_elas.dot(f_elas);
        }

		printf("f_sum error %f \n", err_felas);
		printf("finish solving \n");
  }

}
