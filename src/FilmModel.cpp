#include "Types.h"
#include "FilmModel.h"

#include <glm/glm.hpp>
using namespace glm;

namespace BalloonFEM
{
	double trace(Mat2 F)
	{
		return F[0][0] + F[1][1];
	}

    Mat3x2 Film_neohookean::Piola(Mat3x2 F)
    {
        Mat2 FTF = transpose(F) * F;
        /*  F * inv(trans(F)F) */
        Mat3x2 T = F * inverse( FTF );
        double I3 = determinant( FTF );
        return m_mu*( F - T ) + T * m_lambda * log(I3) / 2.0; 
    }

    Mat3x2 Film_neohookean::StressDiff(Mat3x2 F, Mat3x2 dF)
    {
        Mat2 FTF = transpose(F) * F;
        double I3 = determinant(FTF);
        Mat2 invFTF = inverse(FTF);

        /*  F * inv(trans(F)F) */
        Mat3x2 T = F * invFTF;
        Mat3x2 dT = dF * invFTF - 2.0 * F * invFTF * transpose(F) * dF * invFTF;
        /* d ( log(I3) )*/
        double dlogI3 = 2 * trace(invFTF * transpose(F) * dF);

        return m_mu*(dF - dT) + dT * m_lambda * log(I3) / 2.0 + T * m_lambda * dlogI3 / 2.0;
    }

    Mat3x2 Film_neohookean_3d::Piola(Mat3x2 F)
    {
        Mat2 FTF = transpose(F) * F;
        /*  F * inv(trans(F)F) */
        Mat3x2 T = F * inverse( FTF );
        double I3 = determinant( FTF );
        return m_mu * ( F - T / I3 );
    }

    Mat3x2 Film_neohookean_3d::StressDiff(Mat3x2 F, Mat3x2 dF)
    {
        Mat2 FTF = transpose(F) * F;
		Mat2 invFTF = inverse(FTF);
        /*  F * inv(trans(F)F) */
        Mat3x2 T = F * inverse( FTF );
        double I3 = determinant( FTF );

        Mat3x2 dT = dF * invFTF - 2.0 * F * invFTF * transpose(F) * dF * invFTF;
        /* d ( 1 / I3 ) */
        double dinvI3 = - trace( transpose( 2.0 * T / I3) * dF ); 

		Mat3x2 dP = m_mu * (dF - dT / I3 - T * dinvI3);
		return dP;
    }

    Mat3x2 Film_aniso_neohookean::Piola(Mat3x2 F)
    {
        Mat2 FTF = transpose(F) * F;
        /*  F * inv(trans(F)F) */
        Mat3x2 T = F * inverse( FTF );
        double I3 = determinant( FTF );
        return m_mu*( F * m_sigma - T ) + T * m_lambda * log(I3) / 2.0; 
    }

    Mat3x2 Film_aniso_neohookean::StressDiff(Mat3x2 F, Mat3x2 dF)
    {
        Mat2 FTF = transpose(F) * F;
        double I3 = determinant(FTF);
        Mat2 invFTF = inverse(FTF);

        /*  F * inv(trans(F)F) */
        Mat3x2 T = F * invFTF;
        Mat3x2 dT = dF * invFTF - 2.0 * F * invFTF * transpose(F) * dF * invFTF;
        /* d ( log(I3) )*/
        double dlogI3 = 2 * trace(invFTF * transpose(F) * dF);

        return m_mu*(dF * m_sigma - dT) + dT * m_lambda * log(I3) / 2.0 + T * m_lambda * dlogI3 / 2.0;
    }

    
    Mat3x2 Film_aniso_neohookean_3d::Piola(Mat3x2 F)
    {
        Mat2 FTF = transpose(F) * F;
        /*  F * inv(trans(F)F) */
        Mat3x2 T = F * inverse( FTF );
        double I3 = determinant( FTF );
        return m_mu * ( F * m_sigma - T / I3 );
    }

    Mat3x2 Film_aniso_neohookean_3d::StressDiff(Mat3x2 F, Mat3x2 dF)
    {
        Mat2 FTF = transpose(F) * F;
		Mat2 invFTF = inverse(FTF);
        /*  F * inv(trans(F)F) */
        Mat3x2 T = F * inverse( FTF );
        double I3 = determinant( FTF );

        Mat3x2 dT = dF * invFTF - 2.0 * F * invFTF * transpose(F) * dF * invFTF;
        /* d ( 1 / I3 ) */
        double dinvI3 = - trace( transpose( 2.0 * T / I3) * dF ); 

		Mat3x2 dP = m_mu * (dF * m_sigma - dT / I3 - T * dinvI3);
		return dP;
    }

}
