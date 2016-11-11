#ifndef _DYNAMIC_H_
#define _DYNAMIC_H_

#include <vector>

#include "Types.h"
#include "TetraMesh.h"
#include "Models.h"
#include "ObjState.h"

namespace BalloonFEM
{
	/* take tetramesh and build ObjState and then solve them */
    class Engine
    {
        public:
            Engine(){};
            Engine(
				TetraMesh* tetra, 
				ElasticModel* model = new Elastic_neohookean(0.4, 0.4), 
				AirModel* a_model = new AirModel_Isobaric(0, 0), 
				FilmModel* film_model = new Film_neohookean_3d(0.4, 0.4),
                BendingModel* bend_model = new Bending_MeanCurvature(0.01)
				);

            void setElasticModel(ElasticModel* model){ m_volume_model = model; };
            void setAirModel(AirModel* model){ m_air_model = model; };
			void setFilmModel(FilmModel* model){ m_film_model = model; };
			void setBendingModel(BendingModel* model){ m_bend_model = model; };

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

            /* test the computeForceAndGradient, output force and Gradient
             */
			//SpMat forceTest(Vvec3 &f_sum, Vvec3 &f_loc);

        private:

            TetraMesh* m_tetra;

            ElasticModel* m_volume_model;

            AirModel* m_air_model;

			FilmModel* m_film_model;

            BendingModel* m_bend_model;

			size_t m_size;

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
            void computeForceAndGradient(ObjState &state, SpVec &f, SpMat &A);

            void computeElasticForces(ObjState &state, Vvec3 &f_sum);
            void computeAirForces(ObjState &state, Vvec3 &f_sum);
            void computeFilmForces(ObjState &state, Vvec3 &f_sum);
			void computeFilmForces(ObjState &state, Vvec3 &f_sum, SpMat& Tri);

            /* compute force differntial matrix */
			SpMat computeElasticDiffMat(ObjState &state);
			SpMat computeAirDiffMat(ObjState &state);
			SpMat computeFilmDiffMat(ObjState &state);
			SpMat bendingForceAndGradient(ObjState &state, Vvec3 &f_sum);
    };

}

#endif // _DYNAMIC_H_
