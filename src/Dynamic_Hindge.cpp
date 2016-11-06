#include <vector>

#include <glm/glm.hpp>

#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "Types.h"
#include "Dynamic.h"

namespace BalloonFEM
{
    void Engine::computeBendingForces(ObjState &state, Vvec3 &f_elas)
    {
        /* compute theta */
        Vvec3 &pos = state.world_space_pos;

        size_t count_hindge = 0;
		for (size_t i = 0; i < m_tetra->films.size(); i++)
			count_hindge += m_tetra->films[i].hindges.size();
        SpVec dphi(count_hindge);

        /* compute detrivation of hindge angle theta */
        int offset = 0;
        for(MIter f = m_tetra->films.begin(); f != m_tetra->films.end(); f++)
        {
            for(EIter h = f->hindges.begin(); h != f->hindges.end(); h++)
            {
                /* normal of peice_info[0] */
                iVec3 &id0 = f->peices[h->peice_info[0].x].v_id;
                Vec3 n0 = glm::cross( pos[id0[0]] - pos[id0[2]], pos[id0[1]] - pos[id0[2]] );
                n0 /= glm::length(n0);

                /* normal of peice_info[1] */
                iVec3 &id1 = f->peices[h->peice_info[1].x].v_id;
                Vec3 n1 = glm::cross( pos[id1[0]] - pos[id1[2]], pos[id1[1]] - pos[id1[2]] );
                n1 /= glm::length(n1);

                /* edge direction */
                Vec3 e = pos[id0[h->peice_info[0].y]] -  pos[id0[h->peice_info[0].y]];
                e /= glm::length(e);

                /* energy is 2*sin(x/2)^2 */
                dphi(offset) = glm::dot(e, glm::cross(n0, n1));
                offset ++;
            }
        }

        SpVec bendforce = dphi.transpose() * state.hindgeAngleGradient();

    }

    SpMat Engine::computeBendingDiffMat(ObjState &state)
    {
        SpMat K(m_size, m_size);

        return K;

    }
}
