#include <vector>

#include "Dynamic.h"

namespace BallonFEM
{
    void Engine::computeElasticForces()
    {
        f_elas.assign( m_tetra->vertices.size(), Vec3(0.0));

        for(TIter t = m_tetra->tetrahedrons.begin();
                t != m_tetra->tetrahedrons.end(); t++)
        {
            iVec4 &id = t->v_id;
            Vec3 &v0 = v_pos[id[0]];
            Vec3 &v1 = v_pos[id[1]];
            Vec3 &v2 = v_pos[id[2]];
            Vec3 &v3 = v_pos[id[3]];
            
            /* calculate deformation in world space */
            Mat3 Ds = Mat3(v0 - v3, v1 - v3, v2 - v3);

            /* calculate deformation gradient */
            Mat3 F = Ds * t->Bm;

            /* calculate Piola for this tetra */
            Mat3 P = m_model.Piola(F);

            /* calculate forces contributed from this tetra */
            Mat3 H = - t->W * P * transpose(t->Bm);

            f_elas[id[0]] += H[0];
            f_elas[id[1]] += H[1];
            f_elas[id[2]] += H[2];
            f_elas[id[3]] -= H[0] + H[1] + H[2];
        }
    }

    void Engine::computeForceDifferentials()
    {
        df_elas.assign( m_tetra->vertices.size(), Vec3(0.0));

        for(TIter t = m_tetra->tetrahedrons.begin();
                t != m_tetra->tetrahedrons.end(); t++)
        {
            /* assgin world space position */
            iVec4 &id = t->v_id;
            Vec3 &v0 = v_pos[id[0]];
            Vec3 &v1 = v_pos[id[1]];
            Vec3 &v2 = v_pos[id[2]];
            Vec3 &v3 = v_pos[id[3]];
            
            /* calculate deformation in world space */
            Mat3 Ds = Mat3(v0 - v3, v1 - v3, v2 - v3);

            /* calculate deformation gradient */
            Mat3 F = Ds * t->Bm;
            
            /* assign world space delta position */
            Vec3 &dv0 = dv_pos[id[0]];
            Vec3 &dv1 = dv_pos[id[1]];
            Vec3 &dv2 = dv_pos[id[2]];
            Vec3 &dv3 = dv_pos[id[3]];
            
            /* calculate delta deformation in world space */
            Mat3 dDs = Mat3(dv0 - dv3, dv1 - dv3, dv2 - dv3);

            /* calculate delta deformation gradient */
            Mat3 dF = dDs * t->Bm;

            /* calculate delta Piola */
            Mat3 dP = m_model.StressDiff(F, dF);
 
            /* calculate forces contributed from this tetra */
            Mat3 dH = - t->W * dP * transpose(t->Bm);

            df_elas[id[0]] += dH[0];
            df_elas[id[1]] += dH[1];
            df_elas[id[2]] += dH[2];
            df_elas[id[3]] -= dH[0] + dH[1] + dH[2];     
        }
    }
   
}
