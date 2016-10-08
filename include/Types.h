#ifndef _TYPE_H_
#define _TYPE_H_

#include <vector>
#include <Eigen/Dense>

namespace BallonFEM
{
    /* basic geometry element */
    class Vertex;
    class Face;
    class Tetra;

    /* organized object */
    class TetraMesh;

	typedef Eigen::Vector3d Vec3;
	typedef Eigen::Matrix3d Mat3;
    typedef Eigen::Vector3i iVec3;
    typedef Eigen::Vector4i iVec4;

	/* some definitions */
    typedef Eigen::Vector3d Vec3;
    typedef Eigen::Matrix3d Mat3;

    typedef std::vector<Vertex>::iterator   VIter;
    typedef std::vector<Face>::iterator     FIter;
    typedef std::vector<Tetra>::iterator    TIter;

    typedef std::vector<Vertex>::const_iterator   VCIter;
    typedef std::vector<Face>::const_iterator     FCIter;
    typedef std::vector<Tetra>::const_iterator    TCIter;
}

#endif // !_TYPE_H_
