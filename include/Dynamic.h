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
            Engine(TetraMesh* tetra);

            void outputPosition();

        private:
            TetraMesh* m_tetra;
            ElasticModel m_model;
            
            std::vector<Vec3> f_elas;
            std::vector<Vec3> df_elas;

            std::vector<Vec3> v_pos;
            std::vector<Vec3> dv_pos;

            void computeElasticForces();

            void computeForceDifferentials();
    };

}

#endif // _DYNAMIC_H_
