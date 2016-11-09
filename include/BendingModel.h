#ifndef _BENDINGMODEL_H_
#define _BENDINGMODEL_H_

#include "Types.h"

namespace BalloonFEM
{
    class BendingModel;

    class Bending_MeanCurvature;
    class Bending_Willmore;
    class Bending_Bridson;


    class BendingModel
    {
        public:
            BendingModel(){};
            BendingModel(double a):m_a(a){};
            
			void setAmplitude(double a){ m_a = a; };

            virtual double dphi(double theta, double theta_0){return 0;};

            virtual double ddphi(double theta, double theta_0){return 0;};
    
        protected:
            double m_a;
    };

    /* model by Grinspun, mean curvature diff 
     * Phi = a * (theta - theta_0)^2 
     */
    class Bending_MeanCurvature : public BendingModel
    {
        public:
            Bending_MeanCurvature(){};
            Bending_MeanCurvature(double a):BendingModel(a){};

            double dphi(double theta, double theta_0)
            {
                return 2 * m_a * (theta - theta_0);
            }

            double ddphi(double theta, double theta_0)
            {
                return 2 * m_a;
            }
	};


    /* model by Peter Schroder, using willmore energy 
     * Phi = a * sin^2(theta/2) 
     */
    class Bending_Willmore : public BendingModel
    {
        public:
            Bending_Willmore(){};
            Bending_Willmore(double a):BendingModel(a){};

            double dphi(double theta, double theta_0 = 0 )
            {
                return 0.5 * m_a * sin(theta);
            }

            double ddphi(double theta, double theta_0 = 0)
            {
                return 0.5 * m_a * cos(theta);
            }
	};

    /* model by Bridson
     * Phi = a * (cos(theta/2) - b * theta) 
     */
    class Bending_Bridson : public BendingModel
    {
        public:
            Bending_Bridson(){};
            Bending_Bridson(double a, double b):BendingModel(a), m_b(b){};

            double dphi(double theta, double theta_0)
            {
                return - m_a * ( 0.5 * sin(0.5 * theta) + m_b);
            }

            double ddphi(double theta, double theta_0)
            {
                return -0.25 * m_a * cos(0.5 * theta);
            }

        protected:
            double m_b;
	};

}


#endif //!_BENDINGMODEL_H_
