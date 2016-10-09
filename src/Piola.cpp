#include "Types.h"
#include "Piola.h"

#include <glm/glm.hpp>
using namespace glm;

namespace BallonFEM
{
	float trace(Mat3 F)
	{
		return F[0][0] + F[1][1] + F[2][2];
	}

	Mat3 Piola_linear(Mat3 F, double mu, double lambda)
	{
		Mat3 I = Mat3(1.0);

		return Mat3(mu) * (F + transpose(F) - Mat3( 2 ) * I) + Mat3(lambda * trace(F - I)) * I;
	}

	Mat3 Piola_StVK(Mat3 F, double mu, double lambda)
	{
		Mat3 I = Mat3(1.0);
		Mat3 E = Mat3(0.5) * (transpose(F) * F - I);

		return F * ( Mat3(2 * mu) * E + Mat3(lambda * trace(E)) * I);
	}
}