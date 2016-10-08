#include "Types.h"
#include "Piola.h"

namespace BallonFEM
{

	Mat3 Piola_linear(Mat3 F, double mu, double lambda)
	{
		Mat3 I = Mat3::Identity();

		return mu * (F + F.transpose() - 2 * I) + lambda * (F - I).trace() * I;
	}

	Mat3 Piola_StVK(Mat3 F, double mu, double lambda)
	{
		Mat3 I = Mat3::Identity();
		Mat3 E = 0.5 * (F.transpose() * F - I);

		return F * (2 * mu * E + lambda * E.trace() * I);
	}
}