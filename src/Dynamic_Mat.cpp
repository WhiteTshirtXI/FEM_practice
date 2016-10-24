#include <vector>

#include <glm/glm.hpp>

#include <Eigen/Sparse>

#include "Types.h"
#include "Viewer.h"
#include "controls.h"
#include "Dynamic.h"

namespace BallonFEM
{
  void Engine::computeForceDiffMatrix(ObjState &pos, SpMat &K)
    {
        /* project from constrained freedom state to world space */
        Vvec3 &pos = state.world_space_pos;

        typedef Eigen::Triplet<double> T;
        std::vector<T> coefficients;
        coefficients.reserve( 12 * 12 * pos.size() + );

        /* compute air pressure force differential matrix */
        state.volumeGradientDiffMatrix(K);

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
                pressure.push_back( T(*j, *j, p) );
            }
        }
        SpMat P(pos.size(), pos.size());
        P.setFromTriplets(pressure.begin(), pressure.end());
        K *= P;

        /* compute elastic force Differentials */
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
            
            /* assign world space delta position */
            Vec3 &dv0 = dpos[id[0]];
            Vec3 &dv1 = dpos[id[1]];
            Vec3 &dv2 = dpos[id[2]];
            Vec3 &dv3 = dpos[id[3]];
            
            /* calculate delta deformation in world space */
            Mat3 dDs = Mat3(dv0 - dv3, dv1 - dv3, dv2 - dv3);

            /* calculate delta deformation gradient */
            Mat3 dF = dDs * t->Bm;

            /* calculate delta Piola */
            Mat3 dP = m_model->StressDiff(F, dF);
 
            /* calculate forces contributed from this tetra */
            Mat3 dH = - t->W * dP * transpose(t->Bm);

            df_elas[id[0]] += dH[0];
            df_elas[id[1]] += dH[1];
            df_elas[id[2]] += dH[2];
            df_elas[id[3]] -= dH[0] + dH[1] + dH[2];     
        }
    }
}
