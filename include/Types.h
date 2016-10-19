#ifndef _TYPE_H_
#define _TYPE_H_

#include <vector>
#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

namespace BallonFEM
{
    /* basic geometry element */
    class Vertex;
    class Face;
    class Tetra;

    /* organized object */
	class Rigid;
    class Hole;
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
    typedef std::vector<Rigid>::iterator    RIter;
    typedef std::vector<Hole>::iterator     HIter;

    typedef std::vector<Vertex>::const_iterator   VCIter;
    typedef std::vector<Face>::const_iterator     FCIter;
    typedef std::vector<Tetra>::const_iterator    TCIter;
    typedef std::vector<Rigid>::const_iterator    RCIter;
    typedef std::vector<Hole>::const_iterator     HCIter;
}

#endif // !_TYPE_H_
