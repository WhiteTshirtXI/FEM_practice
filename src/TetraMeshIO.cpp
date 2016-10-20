#include <iostream>

#include "Types.h"
#include "TetraMesh.h"
#include "TetraMeshIO.h"

using namespace std;

namespace BallonFEM{

int TetraMeshIO::read( istream& in, TetraMesh& tetra )
{
    TetraMeshData *data;
	data = new TetraMeshData();
    cout << "Reading from .obj file.\n";
    if ( readTetraData( in, *data ))
    {
        return 1;
    }

    cout << "Building tetra mesh.\n";
    if ( buildTetra( *data, tetra ))
    {
        return 1;
    }

    return 0;
}

int TetraMeshIO::readTetraData(istream& in, TetraMeshData& data)
{
    string line;
   
	data.vertices.clear();
	data.tetrahedrons.clear();
    data.surface.clear();
    data.fixed.clear();
    data.rigid.clear();
    data.holes.clear();

    while( getline( in, line ))
    {
        stringstream ss( line );
        string token;
   
        ss >> token;
   
        if( token == "v"  ) { readPosition( ss, data ); continue; } // vertex
        if( token == "t"  ) { readTetra( ss, data );    continue; } // tetra 
        if( token == "f"  ) { readSurf( ss, data );    continue; } // surface 
		if( token == "x"  ) { readFixed( ss, data );    continue; } // fixed vertices
        if( token == "r"  ) { readRigid( ss, data );    continue; } // rigid body
        if( token == "h"  ) { readHole(  ss, data );    continue; } // hole 
        if( token[0] == '#' ) continue; // comment
        if( token == "o" ) continue; // object name
        if( token == "g" ) continue; // group name
        if( token == "s" ) continue; // smoothing group
        if( token == "mtllib" ) continue; // material library
        if( token == "usemtl" ) continue; // material
        if( token == "" ) continue; // empty string

        cerr << "Error: does not appear to be a valid Wavefront OBJ file!" << endl;
        cerr << "(Offending line: " << line << ")" << endl;
        return 1;
    }

    return 0;
}

void TetraMeshIO::readPosition( stringstream& ss, TetraMeshData& data)
{
    double x, y, z;
    ss >> x >> y >> z;
    data.vertices.push_back( Vec3(x, y, z) );
}

void TetraMeshIO::readTetra( stringstream& ss, TetraMeshData& data)
{
    int x, y, z, w;
    ss >> x >> y >> z >> w;
	data.tetrahedrons.push_back(iVec4(x, y, z, w));
}

void TetraMeshIO::readSurf( stringstream& ss, TetraMeshData& data)
{
    int x, y, z;
    ss >> x >> y >> z;
	data.surface.push_back(iVec3(x, y, z));
}

void TetraMeshIO::readFixed( stringstream& ss, TetraMeshData& data)
{
    size_t f_id;
    while (ss >> f_id)
    {
        data.fixed.push_back(f_id);    
    }
}

void TetraMeshIO::readRigid( stringstream& ss, TetraMeshData& data)
{
    std::vector<size_t> rig;
    size_t f_id;
    while (ss >> f_id)
    {
        rig.push_back(f_id);    
    }

    data.rigid.push_back(rig);
}

void TetraMeshIO::readHole( stringstream& ss, TetraMeshData& data)
{
    size_t id, indexBias = 1;
    int x, y, z;
    ss >> id >> x >> y >> z;
    while (data.holes.size() < id)
    {
        std::vector<iVec3> hole;
        hole.clear();
        data.holes.push_back(hole);
    }

	data.holes[id - indexBias].push_back(iVec3(x, y, z));
}

int TetraMeshIO::buildTetra( TetraMeshData& data, TetraMesh& tetra)
{
    int indexBias = 1;
    /* pre-allocate */
    tetra.vertices.clear();
    tetra.tetrahedrons.clear();
    tetra.surface.clear();
    tetra.fixed_ids.clear();
    tetra.rigid_bodies.clear();
    tetra.holes.clear();

    tetra.vertices.reserve( data.vertices.size() );
    tetra.tetrahedrons.reserve( data.tetrahedrons.size() );
    tetra.surface.reserve( data.surface.size() );

    /* assign vertices */
    for(size_t i = 0; i < data.vertices.size(); i++)
    {
        Vertex nv( data.vertices[i] );
        nv.id = i;
        tetra.vertices.push_back( nv );
    }

    /* assign tetra */
    int v_size = tetra.vertices.size();
    for(size_t i = 0; i < data.tetrahedrons.size(); i++)
    {
        Tetra t;
        iVec4 &indeces = data.tetrahedrons[i];
        for(size_t j = 0; j < 4; j++)
        {
            /* stack overflow */
			if (indeces[j] - indexBias > v_size)
			{
				printf(" tetra index out of vertices range.\n");
				return 1;
			}
            
            t.v[j] = &tetra.vertices[ indeces[j] - indexBias ];
        }

        t.v_id = indeces - indexBias;

        tetra.tetrahedrons.push_back( t );
    }

    /* assign surface */
    for(size_t i = 0; i < data.surface.size(); i++)
    {
        Face f;
        iVec3 &indeces = data.surface[i];
        for(size_t j = 0; j < 3; j++)
        {
            /* stack overflow */
			if (indeces[j] - indexBias > v_size)
			{
				printf(" surface index out of vertices range.\n");
				return 1;
			}
            
            f.v[j] = &tetra.vertices[ indeces[j] - indexBias ];
        }

        f.v_id = indeces - indexBias;

        tetra.surface.push_back( f );
    }

    /* assign fixed vertices */
	tetra.fixed_ids.assign(data.fixed.begin(), data.fixed.end());
    for(size_t i = 0; i < data.fixed.size(); i++)
    {
        tetra.vertices[data.fixed[i] - indexBias].m_fixed = true;
		tetra.fixed_ids[i] -= indexBias;
    }
    

    /* assign rigid bodies */
    for(size_t i = 0; i < data.rigid.size(); i++)
    {
        std::vector<size_t> &rig = data.rigid[i];
        std::vector<size_t> tmp(rig.begin(), rig.end());

        for(size_t j = 0; j < tmp.size(); j++)
            tmp[j] -= indexBias;
        tetra.addRigidBody(tmp);
    }

    /* assign holes */
    for(size_t i = 0; i < data.holes.size(); i++)
    {
        std::vector<iVec3> &h = data.holes[i];
        std::vector<iVec3> tmp(h.begin(), h.end());

        for(size_t j = 0; j < tmp.size(); j++)
            tmp[j] -= indexBias;
        tetra.addHole(tmp);
    }

    return 0;
}

void TetraMeshIO::write( ostream& out, const TetraMesh& tetra)
{
    int indexBias = 1;
    for(VCIter v = tetra.vertices.begin(); v != tetra.vertices.end(); v++)
    {
		Vec3 pos = v->m_pos;
		out << "v " << pos[0] << " " << pos[1] << " " << pos[2] << endl; 
    }

    for(TCIter t = tetra.tetrahedrons.begin(); t !=tetra.tetrahedrons.end(); t++)
    {
        out << "t ";
        out << t->v[0]->id + indexBias << " ";
        out << t->v[1]->id + indexBias << " ";
        out << t->v[2]->id + indexBias << " ";
        out << t->v[3]->id + indexBias << endl;
    }
    
    for(FCIter f = tetra.surface.begin(); f !=tetra.surface.end(); f++)
    {
        out << "f ";
        out << f->v[0]->id + indexBias << " ";
        out << f->v[1]->id + indexBias << " ";
        out << f->v[2]->id + indexBias << endl;
    }

    out << "x";
    for(size_t i = 0; i < tetra.fixed_ids.size(); i++)
        out << " " << tetra.fixed_ids[i] + indexBias;
    out << endl;
    
    for(RCIter r = tetra.rigid_bodies.begin(); r != tetra.rigid_bodies.end(); r++)
    {
        out << "r";
        for(size_t i = 0; i < r->elements.size(); i++)
            out << r->elements[i] + indexBias;
        out << endl;
    }
 
    int count = indexBias;
     for(HCIter h = tetra.holes.begin(); h != tetra.holes.end(); h++)
	 {
		 for (FCIter f = h->holeface.begin(); f != h->holeface.end(); f++)
        {
            out << "h " << count << " ";
			out << f->v[0]->id + indexBias << " ";
			out << f->v[1]->id + indexBias << " ";
			out << f->v[2]->id + indexBias << endl;
        }
		count++;
    }
}

}
