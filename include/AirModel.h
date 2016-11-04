#ifndef _AIRMODEL_H_
#define _ELASTICMODEL_H_

#include "Types.h"

namespace BalloonFEM
{
    class AirModel;

    /* pressure remains constant */
    class AirModel_Isobaric;

    /* no heat transfer with environment */
    class AirModel_Adiabatic;

    class AirModel
    {
        public:

            AirModel(){};

            AirModel(double p, double V):
                m_p(p), m_V(V) {};

            AirModel(const AirModel& other):
                m_p(other.m_p), m_V(other.m_V) {};

            virtual double pressure(double V) {return m_p;};

            virtual double pressureDiff(double V, double dV) {return 0;};

        protected:

            double m_p;
            double m_V;
    };

    class AirModel_Isobaric : public AirModel
    {
        public:

            AirModel_Isobaric(){};

            AirModel_Isobaric(double p, double V):
                AirModel(p, V) {};

            AirModel_Isobaric(const AirModel_Isobaric& other):
                AirModel(other) {};

            double pressure(double V) {return m_p;};

            double pressureDiff(double V, double dV) {return 0;};
    };
  
    class AirModel_Adiabatic : public AirModel
    {
        public:

            AirModel_Adiabatic(){};

            AirModel_Adiabatic(double p, double V, double gamma):
                AirModel(p, V), m_gamma(gamma) {};

            AirModel_Adiabatic(const AirModel_Adiabatic& other):
                AirModel(other), m_gamma(other.m_gamma) {};

            double pressure(double V) 
            {
                return m_p * pow( V / m_V, -m_gamma );
            };

            double pressureDiff(double V, double dV) 
            {
                return - m_gamma * m_p * pow( V / m_V, -m_gamma - 1) * dV / m_V ;
            };

        protected:

            double m_gamma;
    };
}
#endif // !_AIRMODEL_H_
