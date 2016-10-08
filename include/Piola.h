#ifndef _PIOLA_H_
#define _PIOLA_H_

#include "Types.h"

namespace BallonFEM
{
	Mat3 Piola_linear(Mat3 F, double mu, double lambda);
	Mat3 Piola_StVK(Mat3 F, double mu, double lambda);
	Mat3 Piola_corotate(Mat3 F, double mu, double lambda);
	Mat3 Piola_isotropic(Mat3 F, double mu, double lambda);
	Mat3 Piola_neohookean(Mat3 F, double mu, double lambda);
}
#endif //!_PIOLA_H_