#include <vector>

#include <glm/glm.hpp>

#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "Types.h"
#include "Dynamic.h"

namespace{
    using namespace BalloonFEM;
    const Vec3 v[4] = {Vec3(0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    const Mat3 m[4][3] = {
        { Mat3(v[1], v[0], v[0]), Mat3(v[2], v[0], v[0]), Mat3(v[3], v[0], v[0])},
        { Mat3(v[0], v[1], v[0]), Mat3(v[0], v[2], v[0]), Mat3(v[0], v[3], v[0])},
        { Mat3(v[0], v[0], v[1]), Mat3(v[0], v[0], v[2]), Mat3(v[0], v[0], v[3])},
        { -Mat3(v[1], v[1], v[1]), -Mat3(v[2], v[2], v[2]), -Mat3(v[3], v[3], v[3])}
    };

}

namespace BalloonFEM
{
    void Engine::computeElasticForces(ObjState &state, Vvec3 &f_sum)
    {
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
            Mat3 P = m_volume_model->Piola(F);

            /* calculate forces contributed from this tetra */
            Mat3 H = - t->W * P * transpose(t->Bm);

            f_sum[id[0]] += H[0];
            f_sum[id[1]] += H[1];
            f_sum[id[2]] += H[2];
            f_sum[id[3]] -= H[0] + H[1] + H[2];
        }
    }

    SpMat Engine::computeElasticDiffMat(ObjState &state)
    {	
        printf("building elastic differential matrix \n");
        /* project from constrained freedom state to world space */
		Vvec3 &pos = state.world_space_pos;
	    
        /* compute elastic force Differentials */
		std::vector<T> coefficients;
		coefficients.clear();
		coefficients.reserve( 12 * 12 * m_tetra->num_vertex);

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
				Mat3 dP = m_volume_model->StressDiff(F, dF);
 
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

        return E;
    }
}

