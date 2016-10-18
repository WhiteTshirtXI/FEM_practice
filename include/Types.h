#ifndef _TYPE_H_
#define _TYPE_H_

#include <vector>
#include <glm/glm.hpp>

namespace BallonFEM
{
    /* basic geometry element */
    class Vertex;
    class Face;
    class Tetra;

    /* organized object */
	class Rigid;
    class TetraMesh;

	/* some definition */
	typedef glm::dvec3	Vec3;
	typedef glm::dmat3	Mat3;
    typedef glm::ivec3	iVec3;
    typedef glm::ivec4	iVec4;
    typedef glm::dquat  Quat;


    typedef std::vector<Vertex>::iterator   VIter;
    typedef std::vector<Face>::iterator     FIter;
    typedef std::vector<Tetra>::iterator    TIter;

    typedef std::vector<Vertex>::const_iterator   VCIter;
    typedef std::vector<Face>::const_iterator     FCIter;
    typedef std::vector<Tetra>::const_iterator    TCIter;
}

#endif // !_TYPE_H_
