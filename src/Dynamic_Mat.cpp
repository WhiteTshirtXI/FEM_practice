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


namespace BalloonFEM
{

	void Engine::computeForceDiffMat(ObjState &state, SpMat &K)
	{
		printf("building force differential matrix \n");
        
        K = computeAirDiffMat(state);

		K += computeElasticDiffMat(state);
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
