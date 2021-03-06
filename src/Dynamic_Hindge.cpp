#include <vector>
#include <iostream>

#include <glm/glm.hpp>
using namespace glm;

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "Types.h"
#include "Dynamic.h"

namespace{
    using namespace BalloonFEM;
    /* push mat3 coefficients into SpMat,
     * T(target_off + i, source_off + j, H[j][i]) 
     * pay attention that H is column major
     */
    void pushMat3(std::vector<T> &coeff, int target_off, int source_off, Mat3 H)
    {
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
			if (H[j][i] != 0)
                coeff.push_back( T(target_off + i, source_off + j, H[j][i]) );
    }

    /* push mat3 coefficients into SpMat,
     * T(target_off, source_off + i, H[i]) 
     */
    void pushVec3(std::vector<T> &coeff, int target_off, int source_off, Vec3 H)
    {
        for(int i = 0; i < 3; i++)
            coeff.push_back( T(target_off, source_off + i, H[i]) );
    }

    int sgn(double val){
        return ((double(0) < val) - (val < double(0)));
    }
}

namespace BalloonFEM
{
    SpMat Engine::bendingForceAndGradient(ObjState &state, Vvec3 &f_sum)
    {
        /* compute theta */
        Vvec3 &pos = state.world_space_pos;

        SpVec dphi = SpVec::Zero(m_tetra->num_hindges);
        SpMat ddphi(m_tetra->num_hindges, m_tetra->num_hindges);

        ///////////////////////////////////////////////////////////////////////
        /* compute hindge angle theta and phi, dphi, ddphi */
        int offset = 0;
        for(MIter f = m_tetra->films.begin(); f != m_tetra->films.end(); f++)
        {
			std::vector<Piece> &pieces = f->pieces;
            for(EIter h = f->hindges.begin(); h != f->hindges.end(); h++)
            {
                /* normal of piece_info[0] */
                iVec3 &id0 = pieces[h->piece_info[0].x].v_id;
                Vec3 n0 = cross( pos[id0[0]] - pos[id0[2]], pos[id0[1]] - pos[id0[2]] );
                n0 /= length(n0);

                /* normal of piece_info[1] */
                iVec3 &id1 = pieces[h->piece_info[1].x].v_id;
                Vec3 n1 = cross( pos[id1[0]] - pos[id1[2]], pos[id1[1]] - pos[id1[2]] );
                n1 /= length(n1);

                /* edge direction , is the positive direction of piece_info[0]*/
				int i = h->piece_info[0].y;
				int j = (i + 1) % 3, k = (i + 2) % 3;
                Vec3 e = pos[id0[k]] -  pos[id0[j]];
                e /= length(e);

				/* signed theta is defined as positive when n0 n1 point away from each other */
                /* energy is 2*sin(x/2)^2 */
				double tmp = max(min(dot(n0, n1), 1.0), -1.0);
                double theta = acos(tmp) * sgn(dot(e, cross(n0, n1)));
                dphi(offset) = m_bend_model->dphi(theta, h->theta);
                ddphi.coeffRef(offset, offset) = m_bend_model->ddphi(theta, h->theta);
                offset ++;
            }
        }
        
        //////////////////////////////////////////////////////////////////////
        /* compute hessian matrix and theta gradient */
		printf("compute bending force and hessian matrix \n");
        std::vector<T> theta_coeff;
        theta_coeff.clear();
		theta_coeff.reserve( 2 * 9 * m_tetra->num_hindges );

        std::vector<T> hessian_coeff;
        hessian_coeff.clear();
        hessian_coeff.reserve( 9 * 9 * m_tetra->num_pieces );

        /* compute detrivation of hindge angle theta */
        offset = 0;     /* since there may be more than one film obj */
        for(MIter f = m_tetra->films.begin(); f != m_tetra->films.end(); f++)
        {
            for(PIter p = f->pieces.begin(); p != f->pieces.end(); p++)
            {
                iVec3 &v_id = p->v_id;
                /* edge dire , length , norm and area*/
                Vec3 e[3] = { 
                    pos[v_id[2]] - pos[v_id[1]],
                    pos[v_id[0]] - pos[v_id[2]],
                    pos[v_id[1]] - pos[v_id[0]],
                };

                Vec3 norm = cross( e[0], e[1] );
                double A2 = length(norm);
                norm /= A2;

				Vec3 l = Vec3(length(e[0]), length(e[1]), length(e[2]));
				for (int i = 0; i < 3; i++)
					e[i] /= l[i];

                /* height of vertex */
                Vec3 h_inv = l / A2;
				Mat3 w = outerProduct(h_inv, h_inv);

                /* cosine of vertex */

                double cos[3] = { 
                    - dot(e[1], e[2]),
                    - dot(e[2], e[0]),
                    - dot(e[0], e[1]),
                };

                /* edge normal m */
                Vec3 m[3] = { 
                    cross(e[0], norm),
                    cross(e[1], norm),
                    cross(e[2], norm),
                };

                /* intermediate var Mi */
                Mat3 M[3] = { 
                    outerProduct(norm, m[0]), 
					outerProduct(norm, m[1]),
					outerProduct(norm, m[2]),
                };

                /* intermediate var Ni */
                Mat3 N[3] = {
                    M[0] / (l.x * l.x), M[1] / (l.y * l.y), M[2] / (l.z * l.z)
                };

                /* intermediate var ci */
                double c[3] = {0, 0, 0};
                Mat3 R[3] = {Mat3(0), Mat3(0), Mat3(0)};
                for(int i = 0; i < 3; i++)
                    if (p->hindge_id[i] != -1)
                    {
						c[i] = dphi(p->hindge_id[i] + offset);
                        R[i] = c[i]*N[i];
                    }

                /* intermediate var di */
                double d[3] = {0, 0, 0};
                for(int i = 0; i < 3; i++)
                {
                    int j = (i+1) % 3, k = (i+2) % 3;
                    d[i] = c[j] * cos[k] + c[k] * cos[j] - c[i];
                }

                /* for each triangle, compute contribution to Heissen */
                for(int i = 0; i < 3; i++)
                {
                    int j = (i+1) % 3, k = (i+2) % 3;

                    /* j = i */
                    Mat3 H = w[i][i] * d[i] * (M[i] + transpose(M[i])) - R[j] - R[k];

                    pushMat3(hessian_coeff, 3 * v_id[i], 3 * v_id[i], H);

                    /* j = i + 1 */
                    H = w[i][j] * (d[i] * transpose(M[j]) + d[j] * M[i]) + R[k];

                    pushMat3(hessian_coeff, 3 * v_id[i], 3 * v_id[j], H);
                    
                    /* symmetric j = i + 2 */
                    pushMat3(hessian_coeff, 3 * v_id[j], 3 * v_id[i], transpose(H));
                }

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
                        pushVec3( theta_coeff, h_id, 3 * v_id[i], tmp);

                        /* \dev_{x1} theta = n * cos2 / h1 */
                        tmp = norm * cos[k] * h_inv[j];
                        pushVec3( theta_coeff, h_id, 3 * v_id[j], tmp);

                        /* \dev_{x1} theta = n * cos1 / h2 */
                        tmp = norm * cos[j] * h_inv[k];
                        pushVec3( theta_coeff, h_id, 3 * v_id[k], tmp);
                    }
                }

            }
            offset += f->hindges.size();
        }

        /* compute bending force */
        SpMat theta_grad(m_tetra->num_hindges, 3 * m_tetra->num_vertex);
        theta_grad.setFromTriplets(theta_coeff.begin(), theta_coeff.end());

        SpVec bendforce = - dphi.transpose() * theta_grad;

        for(size_t i = 0; i < f_sum.size(); i++)
        {
            f_sum[i] += Vec3(
                    bendforce(3*i    ), 
                    bendforce(3*i + 1), 
                    bendforce(3*i + 2));
        }

        SpMat H(3 * m_tetra->num_vertex, 3 * m_tetra->num_vertex);

        H.setFromTriplets(hessian_coeff.begin(), hessian_coeff.end());

		H +=  theta_grad.transpose() * ddphi * theta_grad;

        return H;
    }
}
