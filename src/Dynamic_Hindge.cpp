#include <vector>

#include <glm/glm.hpp>

#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "Types.h"
#include "Dynamic.h"

namespace BalloonFEM
{
    SpMat Engine::bendingForceAndGradient(ObjState &state, Vvec3 &f_sum)
    {
        /* compute theta */
        Vvec3 &pos = state.world_space_pos;

        size_t count_hindge = 0;
		for (size_t i = 0; i < m_tetra->films.size(); i++)
			count_hindge += m_tetra->films[i].hindges.size();
        SpVec dphi(count_hindge);
        SpMat ddphi(count_hindge, count_hindge);

        ///////////////////////////////////////////////////////////////////////
        /* compute hindge angle theta and phi, dphi, ddphi */
        int offset = 0;
        for(MIter f = m_tetra->films.begin(); f != m_tetra->films.end(); f++)
        {
            for(EIter h = f->hindges.begin(); h != f->hindges.end(); h++)
            {
                /* normal of piece_info[0] */
                iVec3 &id0 = f->pieces[h->piece_info[0].x].v_id;
                Vec3 n0 = glm::cross( pos[id0[0]] - pos[id0[2]], pos[id0[1]] - pos[id0[2]] );
                n0 /= glm::length(n0);

                /* normal of piece_info[1] */
                iVec3 &id1 = f->pieces[h->piece_info[1].x].v_id;
                Vec3 n1 = glm::cross( pos[id1[0]] - pos[id1[2]], pos[id1[1]] - pos[id1[2]] );
                n1 /= glm::length(n1);

                /* edge direction */
                Vec3 e = pos[id0[h->piece_info[0].y]] -  pos[id0[h->piece_info[0].y]];
                e /= glm::length(e);

                /* energy is 2*sin(x/2)^2 */
                dphi(offset) = glm::dot(e, glm::cross(n0, n1));
                ddphi.coeffRef(offset, offset) = glm::dot(n0, n1);
                offset ++;
            }
        }
        
        //////////////////////////////////////////////////////////////////////
        /* compute heissen matrix and theta gradient */
        std::vector<T> theta_coeff;
        theta_coeff.clear();
		theta_coeff.reserve( 2 * 9 * count_hindge);

        /* compute detrivation of hindge angle theta */
        offset = 0;     /* since there may be more than one film obj */
        for(MIter f = m_tetra->films.begin(); f != m_tetra->films.end(); f++)
        {
            for(PIter p = f->pieces.begin(); p != f->pieces.end(); p++)
            {
                iVec3 &v_id = p->v_id;
                /* edge dire , norm and Area */
                Vec3 e[3] = { 
                    pos[v_id[2]] - pos[v_id[1]],
                    pos[v_id[0]] - pos[v_id[2]],
                    pos[v_id[1]] - pos[v_id[0]],
                };

                Vec3 norm = - glm::cross( e[0], e[1] );
                double A2 = glm::length(norm);
                norm /= A2;

                /* height of vertex */
                double h_inv[3] = { 
                    glm::length(e[0]) / A2,
                    glm::length(e[1]) / A2,
                    glm::length(e[2]) / A2,
                };

                /* cosine of vertex */
                for(int i = 0; i < 3; i++)
                    e[i] /= glm::length(e[i]);

                double cos[3] = { 
                    - glm::dot(e[1], e[2]),
                    - glm::dot(e[2], e[0]),
                    - glm::dot(e[0], e[1]),
                };

                /* edge normal m */
                /* Vec3 m[3] = { 
                    glm::cross(e[0], norm),
                    glm::cross(e[1], norm),
                    glm::cross(e[2], norm),
                };*/

                /* for each vertex if its correspond edge is hindge
                 * compute gradient of theta 
                 */
                for(int i = 0; i < 3; i++)
                {
                    if (p->hindge_id[i] != -1)
                    {
                        int h_id = p->hindge_id[i] + offset;

                        int j = (i+1) % 3;
                        int k = (i+2) % 3;

                        Vec3 tmp;

                        /* \dev_{x0} theta = - n / h0 */
                        tmp = - norm * h_inv[i];
                        theta_coeff.push_back( T(h_id, 3 * v_id[i]    , tmp.x) );
                        theta_coeff.push_back( T(h_id, 3 * v_id[i] + 1, tmp.y) );
                        theta_coeff.push_back( T(h_id, 3 * v_id[i] + 2, tmp.z) );

                        /* \dev_{x1} theta = n * cos2 / h1 */
                        tmp = norm * cos[k] * h_inv[j];
                        theta_coeff.push_back( T(h_id, 3 * v_id[j]    , tmp.x) );
                        theta_coeff.push_back( T(h_id, 3 * v_id[j] + 1, tmp.y) );
                        theta_coeff.push_back( T(h_id, 3 * v_id[j] + 2, tmp.z) );

                        /* \dev_{x1} theta = n * cos1 / h2 */
                        tmp = norm * cos[j] * h_inv[k];
                        theta_coeff.push_back( T(h_id, 3 * v_id[k]    , tmp.x) );
                        theta_coeff.push_back( T(h_id, 3 * v_id[k] + 1, tmp.y) );
                        theta_coeff.push_back( T(h_id, 3 * v_id[k] + 2, tmp.z) );
                    }
                }

            }
            offset += f->hindges.size();
        }

        /* compute bending force */
        SpMat theta_grad(count_hindge, m_size);
        theta_grad.setFromTriplets(theta_coeff.begin(), theta_coeff.end());

        SpVec bendforce = dphi.transpose() * state.hindgeAngleGradient();

        for(size_t i = 0; i < f_sum.size(); i++)
        {
            f_sum[i] += Vec3(
                    bendforce(3*i    ), 
                    bendforce(3*i + 1), 
                    bendforce(3*i + 2));
        }

        SpMat K(m_size, m_size);

        K += theta_grad * ddphi * theta_grad.transpose();
        return K;

    }
}
