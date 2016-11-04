#include <vector>
#include <string>
#include <sstream>

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

#include "Types.h"
#include "TetraMesh.h"
#include "ObjState.h"

namespace{
    using namespace BalloonFEM;
    const Vec3 v[3] = {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    /* coeff of dr[i] in
     * d(2An) = (dr1 - dr3) x (r2 -r3) + (r1 - r3) x (dr2 - dr3) 
     */
    const double coeff_v[3][2] = {
        { 1, 0},
        { 0, 1},
        { -1, -1},
    };

}

namespace BalloonFEM
{
    void ObjState::volumeGradientDiffMat(SpMat &M)
    {
        /* clear Matrix */
        M.setZero();
        M.resize( 3 * m_size, 3 * m_size);

        std::vector<T> coefficients;
        coefficients.clear();
        size_t count_holeface = 0;
        for (size_t i = 0; i < m_h_size; i++)
            count_holeface += m_tetra->holes[i].holeface.size();
        coefficients.reserve( 9 * 9 * count_holeface);

        /* compute contributions triangle by triangle */
        Vvec3 &pos = world_space_pos;
        for(size_t i = 0; i < m_h_size; i++)
        {
            Hole &h = m_tetra->holes[i];
            for(FCIter f = h.holeface.begin(); f != h.holeface.end(); f++)
            {
                iVec3 v_id = f->v_id;
                /* j is index of vertex, k is index of dimention */
                for (size_t j = 0; j < 3; j++)
                for (size_t k = 0; k < 3; k++)
                {
                    /* d(2An) = (dr1 - dr3) x (r2 -r3) + (r1 - r3) x (dr2 - dr3) */
                    Vec3 dn = coeff_v[j][0] * glm::cross(v[k], pos[v_id[1]] - pos[v_id[2]])
                        + coeff_v[j][1] * glm::cross(pos[v_id[0]] - pos[v_id[2]], v[k]);
                    
                    for(size_t w = 0; w < 3; w++)
                        for(size_t l = 0; l < 3; l++)
                            coefficients.push_back( T( 3*v_id[w] + l, 3*v_id[j] + k, dn[l]));
                }
            }
        }

        /* since the hole surface normal towards inner, the V is negtive*/
        M.setFromTriplets(coefficients.begin(), coefficients.end());
        M /= -6.0; 
    }


    SpMat DeltaState::projectMat()
    {
        /* present only consider fixed points */
        SpMat W(3 * m_size, 3 * m_size + 6 * m_r_size);

        std::vector<T> coefficients;
        coefficients.clear();
        size_t count_rigid = 0;
        for (size_t i = 0; i < m_r_size; i++)
            count_rigid += m_tetra->rigids[i].elements.size();
        /* each vertex contribute 3 and rigid element contribute 9 */
        coefficients.reserve(3 * m_size + 9 * count_rigid);

        /* vertices, fixed and rigid elements get 0 influence */
        for(VIter vi = m_tetra->vertices.begin(); vi != m_tetra->vertices.end(); vi++)
        {
            if (!vi->m_fixed  &&  vi->rigid == -1)
            {
                coefficients.push_back(T( 3 * vi->id    , 3 * vi->id    , 1));
                coefficients.push_back(T( 3 * vi->id + 1, 3 * vi->id + 1, 1));
                coefficients.push_back(T( 3 * vi->id + 2, 3 * vi->id + 2, 1));
            }
        }

        /* rigid bodies, begin at 3 * m_size */
        size_t count = 3 * m_size;
        for(size_t i = 0; i < m_r_size; i++)
        {
            Rigid &R = m_tetra->rigids[i];
            for(size_t j = 0; j < R.elements.size(); j++)
            {
                size_t id = R.elements[j];
                /* dr = drc + dw x (r - rc) */
                coefficients.push_back(T( 3 * id    , count    , 1));
                coefficients.push_back(T( 3 * id + 1, count + 1, 1));
                coefficients.push_back(T( 3 * id + 2, count + 2, 1));

                Vec3 r = m_state->m_r_pos[i] - m_state->world_space_pos[id];
                coefficients.push_back(T( 3 * id + 1, count + 3, r.z ));
                coefficients.push_back(T( 3 * id + 2, count + 3, -r.y));
                coefficients.push_back(T( 3 * id    , count + 4, -r.z));
                coefficients.push_back(T( 3 * id + 2, count + 4, r.x ));
                coefficients.push_back(T( 3 * id    , count + 5, r.y ));
                coefficients.push_back(T( 3 * id + 1, count + 5, -r.x));
            }
        }

        W.setFromTriplets(coefficients.begin(), coefficients.end());

        return W;
    }

    SpMat DeltaState::restrictedMat()
    {
        SpMat I(3 * m_size + 6 * m_r_size, 3 * m_size + 6 * m_r_size);

        std::vector<T> coefficients;
        coefficients.clear();
        coefficients.reserve(3 * m_size);

        for(VIter vi = m_tetra->vertices.begin(); vi != m_tetra->vertices.end(); vi++)
        {
            if (vi->m_fixed  ||  vi->rigid != -1)
            {
                coefficients.push_back(T( 3 * vi->id    , 3 * vi->id    , 1));
                coefficients.push_back(T( 3 * vi->id + 1, 3 * vi->id + 1, 1));
                coefficients.push_back(T( 3 * vi->id + 2, 3 * vi->id + 2, 1));
            }
        }
        I.setFromTriplets(coefficients.begin(), coefficients.end());

        return I;
    }

    SpVec DeltaState::toSpVec()
    {
        SpVec stateVec( 3 * m_size + 6 * m_r_size );

		size_t count = 0;
        /* output vertices state, fixed and rigid part should be zero */
		for (size_t i = 0; i < m_size; i++)
		{
			stateVec(count    ) = dm_pos[i].x;
			stateVec(count + 1) = dm_pos[i].y;
			stateVec(count + 2) = dm_pos[i].z;
			count += 3;
		}

        /* output rigid body state */
        for (size_t i = 0; i < m_r_size; i++)
        {
			stateVec(count	  ) = dm_r_pos[i].x;
			stateVec(count + 1) = dm_r_pos[i].y;
			stateVec(count + 2) = dm_r_pos[i].z;
			stateVec(count + 3) = dm_r_rot[i].x;
			stateVec(count + 4) = dm_r_rot[i].y;
			stateVec(count + 5) = dm_r_rot[i].z;
			count += 6;
        }

        return stateVec;
    }

    void DeltaState::readSpVec(SpVec& stateVec)
    {
		size_t count = 0;
		/* output vertices state, fixed and rigid part should be zero */
		for (size_t i = 0; i < m_size; i++)
		{
			dm_pos[i].x = stateVec(count);
			dm_pos[i].y = stateVec(count + 1);
			dm_pos[i].z = stateVec(count + 2);
			count += 3;
		}

		/* output rigid body state */
		for (size_t i = 0; i < m_r_size; i++)
		{
			dm_r_pos[i].x = stateVec(count);
			dm_r_pos[i].y = stateVec(count + 1);
			dm_r_pos[i].z = stateVec(count + 2);
			dm_r_rot[i].x = stateVec(count + 3);
			dm_r_rot[i].y = stateVec(count + 4);
			dm_r_rot[i].z = stateVec(count + 5);
			count += 6;
		}
    }
}
