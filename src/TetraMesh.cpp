#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>

#include <glm/glm.hpp>

#include "Types.h"
#include "TetraMesh.h"
#include "TetraMeshIO.h"

using namespace std;

namespace BalloonFEM
{
    
void Face::computeNorm()
{
    Vec3 a = this->v[0]->m_pos - this->v[2]->m_pos;
    Vec3 b = this->v[1]->m_pos - this->v[2]->m_pos;

    this->m_normal = glm::normalize(glm::cross(a, b));
}

void Piece::precomputation()
{
    Mat2 Dm(0);  /* coordinate vectors */
    Vec3 e1 = v[0]->m_cord - v[2]->m_cord;
    Vec3 e2 = v[1]->m_cord - v[2]->m_cord;

    double u = glm::length(e1);
    Dm[0][0] = u;
    Dm[1][0] = glm::dot(e2, e1) / u;
    e2 -= Dm[1][0] * e1 / u;
    Dm[1][1] = glm::length(e2);

    Bm = glm::inverse(Dm);

	/* surface area of face */
    W = glm::determinant(Dm) / 2;

    this->computeNorm();
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

void Film::computeHindges()
{
    /* struct used for storing piece_info and halfedge/edge, then sort them */
    struct hindge_tmp
    {
        hindge_tmp(iVec2 halfedge, iVec2 info)
        {
            piece_info = info;
            if (halfedge[0] > halfedge[1])
            {
                edge[0] = halfedge[1];
                edge[1] = halfedge[0];
            }
			else{
                edge = halfedge;
            }
        };
            
        bool operator<(const hindge_tmp& other)
        {
            if (this->edge[0] < other.edge[0])  return true;

            if (this->edge[0] > other.edge[0])  return false;
          
            return this->edge[1] < other.edge[1];
        };

        iVec2 edge;
        iVec2 piece_info;
    };
    
    std::vector<hindge_tmp> edge_list;
    edge_list.clear();
    edge_list.reserve( this->pieces.size() * 3 );

    for(size_t i = 0; i < this->pieces.size(); i++)
    {
        iVec3 v_id = pieces[i].v_id;
        for(size_t j = 0; j < 3; j++)
        {
            edge_list.push_back( 
                    hindge_tmp(
                        iVec2(v_id[(j+1)%3], v_id[(j+2)%3]),
                        iVec2(i, j)
                    ));
        }
    }

    std::sort(edge_list.begin(), edge_list.end());
    
    /* fill in hindges */
    this->hindges.clear();
    this->hindges.reserve(edge_list.size() / 2);
    size_t count = 0;
    while (count < edge_list.size() - 1)
    {
        if (glm::all(glm::equal(edge_list[count].edge, edge_list[count+1].edge)))
        {
            Hindge h;
            h.piece_info[0] = edge_list[count].piece_info;
            h.piece_info[1] = edge_list[count+1].piece_info;
            this->hindges.push_back(h);
            count += 2;
        }
        else
        {
            count++;
        }
    }

    /* fill in corresponding element of pieces */
    for(size_t i = 0; i < this->hindges.size(); i++)
    {
        Hindge &h = this->hindges[i];
        this->pieces[h.piece_info[0].x].hindge_id[h.piece_info[0].y] = i;
        this->pieces[h.piece_info[1].x].hindge_id[h.piece_info[1].y] = i;
    }

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

    printf("Precomputing properties of films.\n");
    for(MIter m = films.begin(); m != films.end(); m++)
    {
        for(PIter p = m->pieces.begin(); p != m->pieces.end(); p++)
        {
            p->precomputation();
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
		f->computeNorm();
	}

	for (HIter h = holes.begin(); h != holes.end(); h++)
	{
		for (FIter f = h->holeface.begin(); f != h->holeface.end(); f++)
		{
			f->computeNorm();
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

int TetraMesh::addFilm( const std::vector<iVec3>& piece_ids)
{
    Film f;
    f.pieces.clear();
    f.pieces.reserve( piece_ids.size() );
    f.hindges.clear();
    f.hindges.reserve( piece_ids.size() * 3 / 2 );

    for(size_t i = 0; i < piece_ids.size(); i++)
    {
        Piece p;
        /* for each face, build face class and add vertices to hole vertices*/
        p.v_id = piece_ids[i];
        for(size_t j = 0; j < 3; j++)
        {
            p.v[j] = &this->vertices[ p.v_id[j] ];
        }

        f.pieces.push_back(p);
    }

    /* compute hindges */
    f.computeHindges();

    films.push_back(f);
    return 0;
}

int TetraMesh::addRigidBody( const std::vector<size_t>& vertex_ids)
{
    Rigid r = Rigid(vertex_ids);

    Vec3 cord = Vec3(0);
    int r_id = rigids.size();
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

    rigids.push_back(r);

    return 0;
}

int TetraMesh::addHole( const std::vector<iVec3>& face_ids)
{
    Hole h;

    /* initialize hole */
    h.vertices.clear();
    h.vertices.reserve( face_ids.size() + 3);
    h.holeface.clear();
    h.holeface.reserve( face_ids.size() );

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
