#ifndef _ELASTICMODEL_H_
#define _ELASTICMODEL_H_

#include "Types.h"

namespace BalloonFEM
{
    class ElasticModel;
    class Elastic_linear;
    class Elastic_StVK;
    class Elastic_neohookean;
    
    /* base class of elasticity models */
    class ElasticModel
    {
        public:

            ElasticModel(){};

            ElasticModel(double mu, double lambda): 
                m_mu( mu ), m_lambda( lambda ) {};

            ElasticModel(const ElasticModel& other):
                m_mu( other.m_mu ), m_lambda( other.m_lambda ) {};

            void setMu(double mu){ m_mu = mu; };

            void setLambda(double lambda){ m_lambda = lambda; };

            /* Piola P = \frac{ \partial{Energy} }{ \partial{F} }*/
			virtual Mat3 Piola(Mat3 F){ return Mat3(0);  };

            /* Stress difference describe how Piola differs on F+dF */
			virtual Mat3 StressDiff(Mat3 F, Mat3 dF){ return Mat3(0);  };

        protected:
            /* elastic parameters called Lame coefficients
             * mu = k / [2 * (1 + v)]
             * lambda = k * v / [(1+v)*(1-2*v)]
             * k is young's modulus and v is Poisson's ratio
             * */
            double m_mu;
            double m_lambda;
    };


    /* class with linear elasticity model 
     * Piola P(F) = mu * ( F + transpose(F) - 2I ) + lambda * tr(F-I) * I
     * Stress Diff dP(F; dF) = mu * ( dF + transpose(dF) ) + lambda * tr(dF) * I
     * */
    class Elastic_linear : public ElasticModel
    {
        public:

            Elastic_linear(){};

            Elastic_linear(double mu, double lambda):
                ElasticModel(mu, lambda){};

            Elastic_linear(const Elastic_linear& other):
                ElasticModel(other){};

            Mat3 Piola(Mat3 F);

            Mat3 StressDiff(Mat3 F, Mat3 dF);
    };


    /* class with St. Venant-Kirchhoff elasticity model 
     * Green Strain E = [transpose(F) * F - I] / 2
     * Piola P(F) = F * [2 * mu * E + lambda * tr(E) * I]
     * Stress Diff dP(F; dF) = 
     * */
    class Elastic_StVK : public ElasticModel
    {
        public:

            Elastic_StVK(){};

            Elastic_StVK(double mu, double lambda):
                ElasticModel(mu, lambda){};

            Elastic_StVK(const Elastic_StVK& other):
                ElasticModel(other){};

            Mat3 Piola(Mat3 F);

            Mat3 StressDiff(Mat3 F, Mat3 dF);
    };


    /* class with neohooken elasticity model, assuming material isotropic,
     * getting rid of rotation freedom, both in world and material space.
     * three rotation invatiant I_1, I_2, I_3
     * I_1 = tr(trans(F) * F)
     * I_2 = tr( [trans(F) * F]^2 )
     * I_3 = det((trans(F) * F) = (det(F))^2
     * J = sqrt(I_3) = det(F)
     *
     * Piola P(F) = mu * ( F - F^{-T} ) + lambda * log(det(F)) * F^{-T}
     * Stress Diff dP(F; dF) = 
     * */
    class Elastic_neohookean : public ElasticModel
     {
        public:

            Elastic_neohookean(){};

            Elastic_neohookean(double mu, double lambda):
                ElasticModel(mu, lambda){};

            Elastic_neohookean(const Elastic_neohookean& other):
                ElasticModel(other){};

            Mat3 Piola(Mat3 F);

            Mat3 StressDiff(Mat3 F, Mat3 dF);
     };

}
#endif //!_ELASTICMODEL_H_
