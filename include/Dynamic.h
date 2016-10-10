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
            Engine(TetraMesh* tetra, ElasticModel model);

            void setModel(ElasticModel model){ m_model = model; };

            /* output v_pos and v_velocity data to tetra mesh */
            void outputData();

            /* solve v_pos and v_velocity for next timestep */
            void solveNextTimestep(float timestep);

            /* test the comuteElasticForces, deforme v_pos and output f_ela 
             * to v_velocity
             * */
            void forceTest();

        private:
            TetraMesh* m_tetra;
            ElasticModel m_model;

            /* number of tetramesh's vertices */
            size_t  m_size;

            typedef std::vector<Vec3> Vvec3;
            
            /* elastic force on pos, 
             * output of computeElasticForces(Vvec3 &pos) 
             * */
            Vvec3 f_elas;

            /* delta elastic force on pos when disturbed by dpos,
             * output of computeForceDifferentials(Vvec3 &pos, Vvec3 &dpos)
             * */
            Vvec3 df_elas;

            /* current position and current velocity */
            Vvec3 v_pos;
            Vvec3 v_velocity;
            
            /* position and velocity for next time step, need to be solved
             * by backward Euler method
             * */
            Vvec3 v_pos_next;
            Vvec3 v_velo_next;

            /* delta pos_next when iterating to solve position for next time step
             * v_pos_next(k+1) = v_pos_next(k) + dv_pos_next 
             * */
            Vvec3 dv_pos_next;

            /* compute force when given positions pos*/
            void computeElasticForces(Vvec3 &pos);

            /* compute delta force when given position pos and disturbe dpos */
            void computeForceDifferentials(Vvec3 &pos, Vvec3 &dpos);
    };

}

#endif // _DYNAMIC_H_
