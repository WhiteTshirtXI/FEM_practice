#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include <glm/glm.hpp>

#include "Types.h"
#include "TetraMesh.h"
#include "TetraMeshIO.h"

using namespace std;

namespace BallonFEM
{
    
void Face::precomputation()
{
    Vec3 a = this->v[0]->m_pos - this->v[2]->m_pos;
    Vec3 b = this->v[1]->m_pos - this->v[2]->m_pos;

    this->m_normal = glm::normalize(glm::cross(a, b));
}

void Tetra::precomputation()
{
    Mat3 Dm;  /* coordinate vectors */
    Dm = Mat3(v[0]->m_cord, v[1]->m_cord, v[2]->m_cord);
	Dm -= Mat3(v[3]->m_cord, v[3]->m_cord, v[3]->m_cord);

    Bm = glm::inverse(Dm);

	/* volume of the tetrahedron */
    W = glm::determinant(Dm) / 6;
}

void TetraMesh::precomputation()
{
    /* make each tetrahedron correct direction and compute tetra propetry */
    printf("Precomputing properties of tetras.\n");
    for(TIter t = tetrahedrons.begin(); t != tetrahedrons.end(); t++)
    {
        t->precomputation();

        /* if t->W is negtive, then the arange of t-v is not correct */
        if ( t->W < 0 )
        {
            swap(t->v[0], t->v[1]);
            t->precomputation();
        }
        
        /* record vertex id in tetra */
        t->v_id = iVec4(t->v[0]->id, t->v[1]->id, t->v[2]->id, t->v[3]->id );
    }


    /* construct surface */

    printf("Constructing surfaces.\n");
    /* each face has two half faces and by counting we can get these faces */
    typedef set<int> S;
    set< S > existingFaces;
    map< S, int > faceCount;
    map< S, iVec3 > corespondFace;
    
    /* 4 faces of tetrahedron towards outside */
    iVec3 arange[] = { iVec3(0,1,2), iVec3(0,2,3), 
                       iVec3(0,3,1), iVec3(3,2,1),
                     };

    for(TIter t = tetrahedrons.begin(); t != tetrahedrons.end(); t++)
    {
        /* construct 4 faces of tetrahedron */
        for(int i = 0; i < 4; i++)
        {
            iVec3 halfFace = iVec3( t->v[ arange[i][0] ]->id, 
                                    t->v[ arange[i][1] ]->id, 
                                    t->v[ arange[i][2] ]->id
                                    );

            int triangle[] = {halfFace[0], halfFace[1], halfFace[2]}; 
            S face = S( triangle , triangle + 3);

            /* if the face already exists, add up 1; else create and record */
            if (existingFaces.find(face) != existingFaces.end())
            {
				faceCount[face] += 1;
            }
            else{
				existingFaces.insert(face);
				faceCount[face] = 1;
				corespondFace[face] = halfFace;
            }
        }
    }

    this->surface.clear();
    /* find single halffaces and add to tetrahedron's surface */
    for(map< S, int>::iterator m = faceCount.begin(); m != faceCount.end(); m++)
    {
        if (m->second == 1)
        {
            iVec3 &halfFace = corespondFace[m->first];
            Face nf;

            nf.v[0] = &( vertices[ halfFace[0] ] );
            nf.v[1] = &( vertices[ halfFace[1] ] );
            nf.v[2] = &( vertices[ halfFace[2] ] );

            this->surface.push_back( nf );
        }
        else if (m->second > 2)
        {
            iVec3 &halfFace = corespondFace[m->first];
            printf("face conflict on (%d, %d, %d)\n", halfFace[0], halfFace[1], halfFace[2] );
        }
    }

    printf("Total surface triangle %d \n", (int)surface.size());

	recomputeSurfaceNorm();
}

void TetraMesh::recomputeSurfaceNorm()
{
	printf("Recompute properties of surfaces.\n");
	for (FIter f = surface.begin(); f != surface.end(); f++)
	{
		f->precomputation();
	}
}

void TetraMesh::labelFixedId()
{
    fixed_ids.clear();
    for (VIter v = vertices.begin(); v != vertices.end(); v++)
    {
        if (v->m_fixed)
        {
            fixed_ids.push_back(v->id);
        }
    }
}

int TetraMesh::addRigidBody( const std::vector<size_t>& vertex_ids)
{
    Rigid r = Rigid(vertex_ids);

    Vec3 cord = Vec3(0);
    for(size_t i = 0; i < vertex_ids.size(); i++)
    {
        Vertex &v = vertices[vertex_ids[i]];
        if (v.rigid == NULL)
        {
            v.rigid = &r;
            cord += v.m_cord;
        }
        else{
            printf("Conflict with vertex %d.\n", v.id);
            return 1;
        }
    }

	r.m_cord = cord / (double)vertex_ids.size();
    r.m_pos = r.m_cord;

    rigid_bodies.push_back(r);

    return 0;
}

int TetraMesh::read( const string& filename )
{
    string inputFilename = filename;
    ifstream in( filename.c_str() );

    if( !in.is_open() )
    {
        cerr << "Error reading from mesh file " << filename << endl;
        return 1;
    }

    int rval;
    if( !( rval = TetraMeshIO::read( in, *this )))
    {
        this->precomputation();
    }
    return rval;
}

int TetraMesh::write( const string& filename )
{
    ofstream out( filename.c_str() );

    if( !out.is_open() )
    {
        cerr << "Error writing to mesh file " << filename << endl;
        return 1;
    }

    TetraMeshIO::write( out, *this );

    return 0;
}

}
