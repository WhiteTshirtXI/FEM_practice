#include <vector>

#include <glm/glm.hpp>

#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "Types.h"
#include "Dynamic.h"

namespace BalloonFEM
{
    void Engine::computeAirForces(ObjState &state, Vvec3 &f_elas)
    {
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
        
    SpMat Engine::computeAirDiffMat(ObjState &state)
    {
	    printf("building air pressure differential matrix \n");
		/* project from constrained freedom state to world space */
		Vvec3 &pos = state.world_space_pos;

		/* compute air pressure force differential matrix */
		SpMat K = state.volumeGradientDiffMat();

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

        return K;
    }
}
