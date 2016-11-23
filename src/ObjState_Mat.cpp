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
    /* volume gradient of holes */
    SpMat ObjState::volumeGradientDiffMat()
    {
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
        SpMat M( m_size * 3, m_size * 3);
        M.setFromTriplets(coefficients.begin(), coefficients.end());
        M /= -6.0; 

        return M;
    }

    SpMat ObjState::hindgeAngleGradient()
    {
        Vvec3 &pos = world_space_pos;

        std::vector<T> coefficients;
        coefficients.clear();
        size_t count_hindge = 0;
		for (size_t i = 0; i < m_tetra->films.size(); i++)
			count_hindge += m_tetra->films[i].hindges.size();
		coefficients.reserve( 2 * 9 * count_hindge);

        /* compute detrivation of hindge angle theta */
        int offset = 0;
        for(MIter f = m_tetra->films.begin(); f != m_tetra->films.end(); f++)
        {
            for(PIter p = f->pieces.begin(); p != f->pieces.end(); p++)
            {
                iVec3 &id = p->v_id;
                /* edge dire , norm and Area */
                Vec3 e[3] = { 
                    pos[id[2]] - pos[id[1]],
                    pos[id[0]] - pos[id[2]],
                    pos[id[1]] - pos[id[0]],
                };

                Vec3 norm = - glm::cross( e[0], e[1] );
                double A2 = glm::length(norm);
                norm /= A2;

                /* height of vertex */
                double h[3] = { 
                    A2 / glm::length(e[0]),
                    A2 / glm::length(e[1]),
                    A2 / glm::length(e[2]),
                };

                /* cosine of vertex */
                for(int i = 0; i < 3; i++)
                    e[i] /= glm::length(e[i]);

                double cos[3] = { 
                    - glm::dot(e[1], e[2]),
                    - glm::dot(e[2], e[0]),
                    - glm::dot(e[0], e[1]),
                };

                /* for each vertex if its correspond edge is hindge */
                for(int i = 0; i < 3; i++)
                {
                    if (p->hindge_id[i] != -1)
                    {
                        int h_id = p->hindge_id[i];

                        int j = (i+1) % 3;
                        int k = (i+2) % 3;

                        /* \dev_{x0} theta = - n / h0 */
                        coefficients.push_back( T(h_id + offset, 3 * p->v_id[i]    , - norm.x/h[i]) );
                        coefficients.push_back( T(h_id + offset, 3 * p->v_id[i] + 1, - norm.y/h[i]) );
                        coefficients.push_back( T(h_id + offset, 3 * p->v_id[i] + 2, - norm.z/h[i]) );

                        /* \dev_{x1} theta = n * cos2 / h1 */
                        coefficients.push_back( T(h_id + offset, 3 * p->v_id[j]    , norm.x * cos[k]/h[j]) );
                        coefficients.push_back( T(h_id + offset, 3 * p->v_id[j] + 1, norm.y * cos[k]/h[j]) );
                        coefficients.push_back( T(h_id + offset, 3 * p->v_id[j] + 2, norm.z * cos[k]/h[j]) );

                        /* \dev_{x1} theta = n * cos1 / h2 */
                        coefficients.push_back( T(h_id + offset, 3 * p->v_id[k]    , norm.x * cos[j]/h[k]) );
                        coefficients.push_back( T(h_id + offset, 3 * p->v_id[k] + 1, norm.y * cos[j]/h[k]) );
                        coefficients.push_back( T(h_id + offset, 3 * p->v_id[k] + 2, norm.z * cos[j]/h[k]) );
                    }
                }

            }
            offset += f->hindges.size();
        }

        SpMat Theta(offset, m_size);
        Theta.setFromTriplets(coefficients.begin(), coefficients.end());

        return Theta;
    }

    void ObjState::computeProjectMat()
    {
        /* present only consider fixed points */
        m_W = SpMat(3 * m_size, this->kineticDegree());

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

                Vec3 r = m_r_pos[i] - world_space_pos[id];
                coefficients.push_back(T( 3 * id + 1, count + 3, r.z ));
                coefficients.push_back(T( 3 * id + 2, count + 3, -r.y));
                coefficients.push_back(T( 3 * id    , count + 4, -r.z));
                coefficients.push_back(T( 3 * id + 2, count + 4, r.x ));
                coefficients.push_back(T( 3 * id    , count + 5, r.y ));
                coefficients.push_back(T( 3 * id + 1, count + 5, -r.x));
            }
            count += 6;
        }

        m_W.setFromTriplets(coefficients.begin(), coefficients.end());
    }

    void ObjState::computeRestrictedMat()
    {
        m_R = SpMat(3 * m_tetra->num_vertex, this->freedomDegree());

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
        m_R.setFromTriplets(coefficients.begin(), coefficients.end());
    }
}
