#include "Types.h"
#include "ElasticModel.h"

#include <glm/glm.hpp>
using namespace glm;

namespace BalloonFEM
{
	float trace(Mat3 F)
	{
		return F[0][0] + F[1][1] + F[2][2];
	}

    /* Linear Modle */
    
    Mat3 Elastic_linear::Piola(Mat3 F)
	{
		Mat3 I = Mat3(1.0);

        Mat3 P = m_mu * ( F + transpose(F) -  2.0 * I ) 
                 + m_lambda * trace(F - I) * I;

        return P;
	}

    Mat3 Elastic_linear::StressDiff(Mat3 F, Mat3 dF)
    {
        Mat3 I = Mat3(1.0);

        Mat3 dP = m_mu * ( dF + transpose(dF) ) 
                 + m_lambda * trace(dF) * I;
        return dP;
    }


    /* StVK Model */

	Mat3 Elastic_StVK::Piola(Mat3 F)
	{
		Mat3 I = Mat3(1.0);
		Mat3 E = 0.5 * (transpose(F) * F - I);

	    Mat3 P = F * ( 2 * m_mu * E + m_lambda * trace(E) * I);

        return P;
	}

    Mat3 Elastic_StVK::StressDiff(Mat3 F, Mat3 dF)
    {
        Mat3 I = Mat3(1.0);
		Mat3 E = 0.5 * (transpose(F) * F - I);
        Mat3 dE = 0.5 * (transpose(dF) * F + transpose(F) * dF);

        Mat3 dP = dF * ( 2 * m_mu * E + m_lambda * trace(E) * I )
                + F * ( 2 * m_mu * dE + m_lambda * trace(dE) * I );
        
        return dP;
    }


    /* Neohooken Model */

    Mat3 Elastic_neohookean::Piola(Mat3 F)
    {
        Mat3 cF = transpose(inverse(F));

        Mat3 P = m_mu * F 
                + ( m_lambda * log(determinant(F)) - m_mu) * cF;

        return P;
    }

    Mat3 Elastic_neohookean::StressDiff(Mat3 F, Mat3 dF)
    {
        Mat3 cF = transpose(inverse(F));

        Mat3 dP = m_mu * dF 
                + (m_mu - m_lambda * log(determinant(F))) * cF * transpose(dF) * cF
                + m_lambda * trace(transpose(cF)*dF) * cF;

        return dP;
    }

}
