#ifndef _DYNAMIC_H_
#define _DYNAMIC_H_

#include <vector>

#include "Types.h"
#include "TetraMesh.h"
#include "ElasticModel.h"

namespace BallonFEM
{
    class Engine
    {
        public:
            Engine(){};
            Engine(TetraMesh* tetra, ElasticModel* model);

            void setModel(ElasticModel* model){ m_model = model; };

            /* output v_pos and v_velocity data to tetra mesh */
            void outputData();

            /* solve v_pos and v_velocity for next timestep */
            void solveNextTimestep(float timestep);


            /* solve v_pos for quasistatic simulation or static 
             * output static position to v_pos_next 
             */
            void solveStaticPos();

            /* test the comuteElasticForces, deforme v_pos and output f_ela 
             * to v_velocity
             * */
            void forceTest();

        private:
            TetraMesh* m_tetra;
            ElasticModel* m_model;

            /* number of tetramesh's vertices */
            size_t  m_size;

            typedef std::vector<Vec3> Vvec3;
            
            /* current position and current velocity */
            Vvec3 v_pos;
            Vvec3 v_velocity;
            Vvec3 f_ext;
            
            /* position and velocity for next time step, need to be solved
             * by backward Euler method
             * */
            Vvec3 v_pos_next;
            Vvec3 v_velo_next;

            /* compute elastic force when given positions pos*/
            void computeElasticForces(Vvec3 &pos, Vvec3 &f_elas);

            /* compute delta elastic force when given position pos and disturbe dpos */
            void computeForceDifferentials(Vvec3 &pos, Vvec3 &dpos, Vvec3 &df_elas);

            /* compute dot */
            float vvec3Dot(const Vvec3& a, const Vvec3& b)
            {
                float count = 0;
                for (size_t i = 0; i < m_size; i++)
                    count += glm::dot( a[i], b[i] );
                return count;
            }
    };

}

#endif // _DYNAMIC_H_
