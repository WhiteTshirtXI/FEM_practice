#ifndef _TETRAMESHIO_H_
#define _TETRAMESHIO_H_

/*  TetraMeshIO hadles input/output operations for TetraMesh objects.
 *
*/

#include <iosfwd>
#include <string>
#include <sstream>
#include <vector>
#include <set>

#include "Types.h"
#include "TetraMesh.h"

namespace BallonFEM{

    class TetraMeshData{
        public:
            TetraMeshData(){};
            std::vector<Vec3> vertices;
            std::vector<iVec4> tetrahedrons;
            std::vector<iVec3> surface;
            std::vector<size_t> fixed;
            std::vector<std::vector<size_t>> rigid;
            std::vector<std::vector<iVec3>> holes;

    };

    class TetraMeshIO
    {
        public:
            /* read a tetra mesh from a input stream */
            static int read( std::istream& in, TetraMesh& tetra);

            /* write a tetra mesh to a output stream*/
            static void write( std::ostream& out, const TetraMesh& tetra);

        protected:
            /* read tetra mesh data */
            static int readTetraData( std::istream& in, TetraMeshData& data);

            static int buildTetra( TetraMeshData& data, TetraMesh& tetra);

            /* when begin with v, read vertex position */
            static void readPosition( std::stringstream& ss, TetraMeshData& data );

            /* when begin with t, read tetra info */
            static void readTetra   ( std::stringstream& ss, TetraMeshData& data );

            /* when begin with f, read surface info */
            static void readSurf   ( std::stringstream& ss, TetraMeshData& data );

            /* when begin with x, read fixed vertices */
            static void readFixed   ( std::stringstream& ss, TetraMeshData& data );

            /* when begin with r, read rigid body by vertices */
            static void readRigid   ( std::stringstream& ss, TetraMeshData& data );

            /* when begin with r, read rigid body by vertices */
            static void readHole    ( std::stringstream& ss, TetraMeshData& data );
    };
}


#endif

