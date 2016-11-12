#include <vector>

#include <glm/glm.hpp>

#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "Types.h"
#include "Dynamic.h"

namespace{
    using namespace BalloonFEM;
    const Vec3 v[4] = {Vec3(0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    const Mat3x2 m[3][3] = {
        { Mat3x2(v[1], v[0]), Mat3x2(v[2], v[0]), Mat3x2(v[3], v[0])},
        { Mat3x2(v[0], v[1]), Mat3x2(v[0], v[2]), Mat3x2(v[0], v[3])},
        { -Mat3x2(v[1], v[1]), -Mat3x2(v[2], v[2]), -Mat3x2(v[3], v[3])}
    };

    void pushMat3x2(std::vector<T> &coeff, iVec3 v_id, int f_id, Mat3x2 H)
    {
        Mat3 F(H[0], H[1], -H[0]-H[1]);
        for (size_t i = 0; i < 3; i++)
            for(size_t j = 0; j < 3; j++)
            {
                coeff.push_back( T( 3 * v_id[i] + j, f_id, F[i][j]) );
            }
    }
}

namespace BalloonFEM
{
    void Engine::computeFilmForces(ObjState &state, Vvec3 &f_sum)
    {
        Vvec3 &pos = state.world_space_pos;

		/* compute film elastic force */

		for (MIter f = m_tetra->films.begin(); f != m_tetra->films.end(); f++)
		{
			for (PIter p = f->pieces.begin(); p != f->pieces.end(); p++)
			{
				iVec3 &id = p->v_id;
				Vec3 &v0 = pos[id[0]];
				Vec3 &v1 = pos[id[1]];
				Vec3 &v2 = pos[id[2]];

				/* calculate deformation in world space */
				Mat3x2 Ds = Mat3x2(v0 - v2, v1 - v2);

				/* calculate deformation gradient */
				Mat3x2 F = Ds * p->Bm;

				/* calculate Piola for this tetra */
				Mat3x2 P = m_film_model->Piola(F);

				/* calculate forces contributed from this tetra */
				Mat3x2 H = - p->volume() * P * transpose(p->Bm);

				f_sum[id[0]] += H[0];
				f_sum[id[1]] += H[1];
				f_sum[id[2]] -= H[0] + H[1];
			}
		}
    }
     
    SpMat Engine::computeFilmDiffMat(ObjState &state)
    {
       	printf("building film force differential matrix \n");
		/* project from constrained freedom state to world space */
		Vvec3 &pos = state.world_space_pos; 

	    std::vector<T> coefficients;
		coefficients.clear();
		size_t count_film = 0;
		for (size_t i = 0; i < m_tetra->films.size(); i++)
			count_film += m_tetra->films[i].pieces.size();
		coefficients.reserve( 9 * 9 * count_film);

		for (MIter f = m_tetra->films.begin(); f != m_tetra->films.end(); f++)
        {
			for (PIter p = f->pieces.begin(); p != f->pieces.end(); p++)
			{
				/* assgin world space position */
				iVec3 &id = p->v_id;
				Vec3 &v0 = pos[id[0]];
				Vec3 &v1 = pos[id[1]];
				Vec3 &v2 = pos[id[2]];

				/* calculate deformation in world space */
				Mat3x2 Ds = Mat3x2(v0 - v2, v1 - v2);

				/* calculate deformation gradient */
				Mat3x2 F = Ds * p->Bm;

				/* i is index of vertex, j is index of dimention */
				for (size_t i = 0; i < 3; i++)
				for (size_t j = 0; j < 3; j++)
				{
					/* calculate delta deformation in world space */
					Mat3x2 dDs = m[i][j];

					/* calculate delta deformation gradient */
					Mat3x2 dF = dDs * p->Bm;

					/* calculate delta Piola */
					Mat3x2 dP = m_film_model->StressDiff(F, dF);

					/* calculate forces contributed from this tetra */
					Mat3x2 dH = - p->volume() * dP * transpose(p->Bm);

                    pushMat3x2(coefficients, id, id[i] + j, dH);
                    /*
					for (size_t w = 0; w < 2; w++)
					for (size_t l = 0; l < 3; l++)
						coefficients.push_back(T(3 * id[w] + l, 3 * id[i] + j, dH[w][l]));

					Vec3 df_3 = - dH[0] - dH[1];
					for (size_t l = 0; l < 3; l++)
						coefficients.push_back(T(3 * id[2] + l, 3 * id[i] + j, df_3[l]));
                    */
				}
			}
		}
		
		SpMat E( 3 * pos.size(), 3 * pos.size());
		E.setFromTriplets(coefficients.begin(), coefficients.end());

        return E;
    }
}
