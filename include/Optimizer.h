#ifndef _OPTIMIZER_H_
#define _OPTIMIZER_H_

#include "Types.h"
#include "TetraMesh.h"
#include "Models.h"
#include "Dynamic.h"


namespace BalloonFEM
{
    class Optimizer : public Engine
    {
        public:
            Optimizer():Engine(){};
            Optimizer(
				TetraMesh* tetra, 
				ElasticModel* model = new Elastic_neohookean(0.4, 0.4), 
				AirModel* a_model = new AirModel_Isobaric(0, 0), 
				FilmModel* film_model = new Film_neohookean_3d(0.4, 0.4),
                BendingModel* bend_model = new Bending_MeanCurvature(0.01)
				):Engine(tetra, model, a_model, film_model, bend_model){};

            void solveOptimal();

        protected:

            void computeForceAndGradient(ObjState &state, SpVec &f, SpMat &A);

			void computeFilmForces(ObjState &state, Vvec3 &f_sum, SpMat& Tri);

        private:
            const double m_alpha = 1.0;
            const double m_beta = 1.0;
            const double m_gamma = 1.0;
    };


}


#endif //! _OPTIMIZER_H_
