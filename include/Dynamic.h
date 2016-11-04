#ifndef _DYNAMIC_H_
#define _DYNAMIC_H_

#include <vector>

#include "Types.h"
#include "TetraMesh.h"
#include "ElasticModel.h"
#include "AirModel.h"
#include "ObjState.h"

namespace BalloonFEM
{
	/* take tetramesh and build ObjState and then solve them */
    class Engine
    {
        public:
            Engine(){};
            Engine(TetraMesh* tetra, ElasticModel* model, AirModel* a_model);

            void setElasticModel(ElasticModel* model){ m_model = model; };

            void setAirModel(AirModel* model){ m_a_model = model; };

			/* input data v_pos, v_velocity and f_ext */
			void inputData();

            /* output v_pos and v_velocity data to tetra mesh */
            void outputData();

            /* copy v_pos_next, to v_pos*/
            void stepToNext();

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

            AirModel* m_a_model;

			size_t m_size;

            typedef std::vector<Vec3> Vvec3;

            /* current position and current velocity */
            ObjState cur_state;
            Vvec3 f_ext;

            /* constrains */
            std::vector<size_t> fixed_id;
            
            /* position and velocity for next time step, need to be solved
             * by backward Euler method
             * */
            ObjState next_state;

            /* compute elastic force when given positions pos*/
            void computeForces(ObjState &state, Vvec3 &f_elas);

            void computeElasticForces(ObjState &state, Vvec3 &f_elas);
            void computeAirForces(ObjState &state, Vvec3 &f_elas);
            void computeFilmForces(ObjState &state, Vvec3 &f_elas);

            /* compute force differntial matrix */
			void computeForceDiffMat(ObjState &state, SpMat& K);

			SpMat computeElasticDiffMat(ObjState &state);
			SpMat computeAirDiffMat(ObjState &state);
			SpMat computeFilmDiffMat(ObjState &state);
    };

}

#endif // _DYNAMIC_H_
