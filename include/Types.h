#ifndef _TYPE_H_
#define _TYPE_H_

#include <vector>

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace BalloonFEM
{
    /* basic geometry element */
    class Vertex;
    class Face;
    class Piece;
    class Hindge;
    class Tetra;

    /* organized object */
    class Film;
	class Rigid;
    class Hole;
    class TetraMesh;

	/* some definition */
	typedef glm::dvec3	    Vec3;
	typedef glm::dmat3	    Mat3;
    typedef glm::dmat2      Mat2;
    typedef glm::dmat2x3    Mat3x2;
    typedef glm::ivec2      iVec2;
    typedef glm::ivec3	    iVec3;
    typedef glm::ivec4	    iVec4;
    typedef glm::dquat      Quat;

	typedef std::vector<Vec3> Vvec3;
	typedef std::vector<Quat> Vquat;

    /* Eigen vec & mat definition */
    typedef Eigen::VectorXd                 SpVec;
    typedef Eigen::SparseMatrix<double>     SpMat;
    typedef Eigen::Triplet<double>          T;

    /* iterator definition of geometry elements */
    typedef std::vector<Vertex>::iterator   VIter;
    typedef std::vector<Face>::iterator     FIter;
    typedef std::vector<Piece>::iterator    PIter;
    typedef std::vector<Hindge>::iterator   EIter;
    typedef std::vector<Tetra>::iterator    TIter;
    typedef std::vector<Film>::iterator     MIter;
    typedef std::vector<Rigid>::iterator    RIter;
    typedef std::vector<Hole>::iterator     HIter;

    typedef std::vector<Vertex>::const_iterator   VCIter;
    typedef std::vector<Face>::const_iterator     FCIter;
    typedef std::vector<Piece>::const_iterator    PCIter;
    typedef std::vector<Hindge>::const_iterator   ECIter;
    typedef std::vector<Tetra>::const_iterator    TCIter;
    typedef std::vector<Film>::const_iterator     MCIter;
    typedef std::vector<Rigid>::const_iterator    RCIter;
    typedef std::vector<Hole>::const_iterator     HCIter;
}

#endif // !_TYPE_H_
