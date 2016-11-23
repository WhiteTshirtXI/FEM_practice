#ifndef _OPTIMIZER_H_
#define _OPTIMIZER_H_

#include "Types.h"
#include "TetraMesh.h"
#include "Models.h"
#include "ObjState.h"
#include "Dynamic.h"


namespace BalloonFEM
{
	class OptState : public ObjState
	{
		public:
			OptState(){};
			OptState(TetraMesh* tetra): ObjState(tetra){};

			void update(SpVec dpos);

			size_t freedomDegree();
	};

    class Optimizer : public Engine
    {
        public:
            Optimizer():Engine(){};
            Optimizer(
				TetraMesh* tetra, 
				TetraMesh* target,
				ElasticModel* model = new Elastic_neohookean(0.4, 0.4), 
				AirModel* a_model = new AirModel_Isobaric(0, 0), 
				FilmModel* film_model = new Film_neohookean_3d(0.4, 0.4),
                BendingModel* bend_model = new Bending_MeanCurvature(0.01));

			/* solve optimal problem with newton method */
            void solveOptimal();

			/* solve optimal problem with Gaussian Newton method */
			void solveOptimalGN();

			void testFunc();
			TetraMesh* Target(){ return m_target;; }

			void setCoeff(Vec3 coeff){ m_alpha = coeff.x; m_beta = coeff.y; m_gamma = coeff.z; };
        protected:

            /* compute thickness lap m_L, which compute the difference between
             * thickness of current face and mean thickness of neighbors */
			void computeThicknessLap();

            /* compute film force and the gradient of thickness and aniso_sigma on f_real */
			void computeFilmForces(ObjState &state, Vvec3 &f_sum, SpMat& Thk, SpMat& Sigma);

			/* compute gradient and heissen matrix, return energy value */
            double computeGradientAndHessian(ObjState &state, ObjState &target, SpVec &f, SpMat &A);

			/* compute residual and jacobian matrix, return energy value */
			double computeResidualAndJacobian(ObjState &state, ObjState &target, SpVec &f, SpMat &A);

			TetraMesh* m_target;
			ObjState* target_state;

			SpMat m_L;
            double m_alpha = 10.0;
            double m_beta = 1.0;
            double m_gamma = 100.0;
			double m_penalty = 10;
    };


}


#endif //! _OPTIMIZER_H_
