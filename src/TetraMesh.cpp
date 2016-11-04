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

namespace BalloonFEM
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
            printf("Invalid tetra vertices order! (%d %d %d %d)",
                    t->v_id.x, t->v_id.y, t->v_id.z, t->v_id.w);
            swap(t->v[0], t->v[1]);
            t->precomputation();
            /* record vertex id in tetra */
            t->v_id = iVec4(t->v[0]->id, t->v[1]->id, t->v[2]->id, t->v[3]->id );
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

	for (HIter h = holes.begin(); h != holes.end(); h++)
	{
		for (FIter f = h->holeface.begin(); f != h->holeface.end(); f++)
		{
			f->precomputation();
		}
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
    int r_id = rigid_bodies.size();
    for(size_t i = 0; i < vertex_ids.size(); i++)
    {
        Vertex &v = vertices[vertex_ids[i]];
        if (v.rigid == -1)
        {
            v.rigid = r_id;
            cord += v.m_cord;
        }
        else{
            printf("Conflict with vertex %d.\n", (int)v.id);
            return 1;
        }
    }

	r.m_cord = cord / (double)vertex_ids.size();
    r.m_pos = r.m_cord;

    rigid_bodies.push_back(r);

    return 0;
}

int TetraMesh::addHole( const std::vector<iVec3>& face_ids)
{
    Hole h;

    /* initialize hole */
    h.vertices.clear();
    h.holeface.clear();

    int h_id = holes.size();
    for(size_t i = 0; i < face_ids.size(); i++)
    {
        /* for each face, build face class and add vertices to hole vertices*/
        Face f;
        f.v_id = face_ids[i];
        for(size_t j = 0; j < 3; j++)
        {
            f.v[j] = &this->vertices[ f.v_id[j] ];
            if (f.v[j]->hole == -1)
            {
                f.v[j]->hole = h_id;
                h.vertices.push_back(f.v_id[j]);
            }
            else if (f.v[j]->hole != h_id){
                printf("Hole %d conflict with vertex %d.\n", h_id, (int)f.v_id[j]);
                return 1;
            }
        }

        h.holeface.push_back(f);
    }

    holes.push_back(h);
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

int TetraMesh::write( const string& filename, const string& additioninfo)
{
    ofstream out( filename.c_str() );

    if( !out.is_open() )
    {
        cerr << "Error writing to mesh file " << filename << endl;
        return 1;
    }

	out << additioninfo.c_str() << endl;
    TetraMeshIO::write( out, *this );

    return 0;
}

}
