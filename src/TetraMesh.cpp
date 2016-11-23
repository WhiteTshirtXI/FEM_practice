#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>

#include <glm/glm.hpp>
#include <glm/gtx/euler_angles.hpp>

#include "Types.h"
#include "TetraMesh.h"
#include "TetraMeshIO.h"

using namespace std;

namespace{
	/* return 1 when val >= 0, -1 when val < 0 */
	int sgn(double val){
		return ((double(0) <= val) - (val < double(0)));
	}

	double clap(double val){
		if (val < -1) return -1;
		if (val > 1) return 1;
		return val;
	}
}

namespace BalloonFEM
{
    
void Face::computeNorm()
{
    Vec3 a = this->v[0]->m_pos - this->v[2]->m_pos;
    Vec3 b = this->v[1]->m_pos - this->v[2]->m_pos;

    this->m_normal = glm::normalize(glm::cross(a, b));
}

Vec3 Face::color()
{
	return Vec3(1, 1, 1);
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

Vec3 Piece::color()
{
	const double min_value = 1;
	const double max_value = 4;
	const double resolution = 1e-4;

	Vec3 rgb(0, 0, 0);

	double value = aniso_sigma[1]/aniso_sigma[0];
	if (value < min_value)
		return Vec3(0,0,0);
	else if (value > max_value)
		return Vec3(1, 1, 1);

	float lamda = (value - min_value) / (max_value - min_value);
	if (resolution > 0.001)
		lamda = int(lamda / resolution) * resolution;

	float coef = lamda * 4.0f;

	if (coef <= 1.0f)
		rgb[0] = 0.0f, rgb[1] = coef, rgb[2] = 1.0f;
	else if (coef <= 2.0f)
		rgb[0] = 0.0f, rgb[1] = 1.0f, rgb[2] = 2.0f - coef;
	else if (coef <= 3.0f)
		rgb[0] = coef - 2.0f, rgb[1] = 1.0f, rgb[2] = 0.0f;
	else
		rgb[0] = 1.0f, rgb[1] = 4.0f - coef, rgb[2] = 0.0f;
	
	return rgb;
}

void Piece::computeStretch()
{
    Mat3x2 F = Mat3x2(v[0]->m_pos - v[2]->m_pos, v[1]->m_pos - v[2]->m_pos);
    F = F * Bm;

    Mat2 E = glm::transpose(F) * F;

    double tr = E[0][0] + E[1][1];
    double det = glm::determinant(E);
    double e = sqrt(tr*tr - 4 * det);
    if (e == 0)		/* if stretch is isotropic */
    {
		this->stretch_angle = 0;
        stretch = F;
    }
	else
    {
		/* take care if E[0][1] == 0, sgn should return 1, 
		 * this makes theta == pi case corrcet,
		 * or sigma_1 < sigma_2 
		 */
		this->stretch_angle = sgn(E[0][1]) * std::acos( clap((E[0][0] - E[1][1]) / e) ) / 2.0;
        Mat2 V = glm::orientate2(stretch_angle);
        stretch = F * V;
    }
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

		/* compute correspond angle */
		for (EIter h = m->hindges.begin(); h != m->hindges.end(); h++)
		{
			/* normal of piece_info[0] */
			Piece &p0 = m->pieces[h->piece_info[0].x];
			iVec3 &id0 = p0.v_id;
			Vec3 n0 = glm::cross(
				vertices[id0[0]].m_cord - vertices[id0[2]].m_cord,
				vertices[id0[1]].m_cord - vertices[id0[2]].m_cord
				);
			n0 /= glm::length(n0);

			/* normal of piece_info[1] */

			Piece &p1 = m->pieces[h->piece_info[1].x];
			iVec3 &id1 = p1.v_id;
			Vec3 n1 = glm::cross(
				vertices[id1[0]].m_cord - vertices[id1[2]].m_cord,
				vertices[id1[1]].m_cord - vertices[id1[2]].m_cord
				);
			n1 /= glm::length(n1);

			/* edge direction , is the positive direction of piece_info[0]*/
			int i = h->piece_info[0].y;
			int j = (i + 1) % 3, k = (i + 2) % 3;
			Vec3 e = vertices[id0[k]].m_cord - vertices[id0[j]].m_cord;
			e /= glm::length(e);

			/* signed theta is defined as positive when n0 n1 point away from each other */
			/* energy is 2*sin(x/2)^2 */
			double tmp = max(min(dot(n0, n1), 1.0), -1.0);
			h->theta = acos(tmp) * sgn(dot(e, cross(n0, n1)));
		}
    }

    printf("Total surface triangle %d \n", (int)surface.size());
	recomputeSurfaceNorm();

    /* compute mesh info */
    num_vertex = vertices.size();
    num_tetra = tetrahedrons.size();

    num_pieces = 0;
    num_hindges = 0;
    for(MIter f = films.begin(); f != films.end(); f++)
    {
        num_pieces += f->pieces.size();
        num_hindges += f->hindges.size();
    }

    num_holeface = 0;
    for(HIter h = holes.begin(); h != holes.end(); h++)
        num_holeface += h->holeface.size();
}

void TetraMesh::recomputeSurfaceNorm()
{
	printf("Recompute properties of surfaces.\n");
	for (FIter f = surface.begin(); f != surface.end(); f++)
	{
		f->computeNorm();
	}

	for (HIter h = holes.begin(); h != holes.end(); h++)
		for (FIter f = h->holeface.begin(); f != h->holeface.end(); f++)
		{
			f->computeNorm();
		}

    for (MIter f = films.begin(); f != films.end(); f++)
        for(PIter p = f->pieces.begin(); p != f->pieces.end(); p++)
        {
            p->computeNorm();
            p->computeStretch();
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

int TetraMesh::addFilm( const std::vector<iVec3>& piece_ids, const std::vector<PieceData>& piece_data)
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
        PieceData tmp = piece_data[i];
        p.h = tmp.h;
        p.aniso_angle = tmp.aniso_angle;
        p.aniso_sigma[0] = tmp.aniso_sigma_1;
        p.aniso_sigma[1] = tmp.aniso_sigma_2;

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
