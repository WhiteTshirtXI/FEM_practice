#ifndef _FILMMODEL_H_
#define _FILMMODEL_H_

#include "Types.h"

namespace BalloonFEM
{
    class FilmModel;

    class Film_neohookean;
    class Film_neohookean_3d;

    /* base class of elasticity models */
    class FilmModel
    {
        public:

            FilmModel(){};

            FilmModel(double mu, double lambda): 
                m_mu( mu ), m_lambda( lambda ) {};

            FilmModel(const FilmModel& other):
                m_mu( other.m_mu ), m_lambda( other.m_lambda ) {};

            void setMu(double mu){ m_mu = mu; };

            void setLambda(double lambda){ m_lambda = lambda; };

            /* Piola P = \frac{ \partial{Energy} }{ \partial{F} }*/
			virtual Mat3x2 Piola(Mat3x2 F){ return Mat3x2(0);  };

            /* Stress difference describe how Piola differs on F+dF */
			virtual Mat3x2 StressDiff(Mat3x2 F, Mat3x2 dF){ return Mat3x2(0);  };

        protected:
            /* elastic parameters called Lame coefficients
             * mu = k / [2 * (1 + v)]
             * lambda = k * v / [(1+v)*(1-2*v)]
             * k is young's modulus and v is Poisson's ratio
             * */
            double m_mu;
            double m_lambda;
    };    

     /* class with neohooken elasticity model, assuming material isotropic,
     * getting rid of rotation freedom, both in world and material space.
     * three rotation invatiant I_1, I_2, I_3
     * I_1 = tr(trans(F) * F)
     * I_2 = tr( [trans(F) * F]^2 )
     * I_3 = det((trans(F) * F)
     * J = sqrt(I_3) = det(F)
     * \Phi = mu/2 * (I1 - log(I3) - 3) + lambda/8 * log(I3)^2 
     * Piola = mu (F - F * inv(trans(F)F)) + lambda/2 log(I3) * F * inv(trans(F)F)
     */
    class Film_neohookean : public FilmModel
    {
        public:

            Film_neohookean(){};

            Film_neohookean(double mu, double lambda):
                FilmModel(mu, lambda){};

            Film_neohookean(const Film_neohookean& other):
                FilmModel(other){};

            Mat3x2 Piola(Mat3x2 F);

            Mat3x2 StressDiff(Mat3x2 F, Mat3x2 dF);
     };

   
    /* neohookean model, but projected to 3D assume the volume of 
     * film does not change, F is a Mat3x2
     * \Phi = mu/2 * (I1 - 3) = mu/2 (F:F + det(trans(F) * F) - 3)
     */
    class Film_neohookean_3d : public FilmModel
     {
        public:

            Film_neohookean_3d(){};

            Film_neohookean_3d(double mu, double lambda):
                FilmModel(mu, lambda){};

            Film_neohookean_3d(const Film_neohookean_3d& other):
                FilmModel(other){};

            Mat3x2 Piola(Mat3x2 F);

            Mat3x2 StressDiff(Mat3x2 F, Mat3x2 dF);
     };

}

#endif // !_FILMMODEL_H_
