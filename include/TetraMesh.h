#ifndef _TETRAMESH_H_
#define _TETRAMESH_H_

#include <vector>
#include <string>
#include <Eigen/Dense>

#include "Types.h"

namespace BalloonFEM{

    typedef Vertex* pVertex;

    class Vertex
    {
        public:
            Vertex();

            Vertex( Vec3 cord ): 
                m_cord( cord ), m_pos( cord ) {};

            Vertex( Vec3 cord, Vec3 pos ): 
                m_cord( cord ), m_pos( pos ) {};

            /* geometry property */
            Vec3 m_cord;                /* position in material coordinate */
            Vec3 m_pos;                 /* position in world space */
            Vec3 m_velocity = Vec3(0);  /* velocity in world space */
            Vec3 m_f_ext = Vec3(0);     /* exterio force */
            
            /* topology property */
            size_t id;                  /* index in tetramesh vertex vector */
            bool m_fixed = false;       /* if is fixed */
            int rigid = -1;             /* connected Rigid id, -1 if not */
            int hole = -1;              /* connected Hole id, -1 if not */
    };

    class Face
    {
        public:
            Face(){};
            
            /* compute normal */
            void computeNorm();
            
            /* geometry property */
            Vec3 m_normal;  

            /* topology property */
            pVertex v[3];   /* pointer to vertices */
            iVec3 v_id;     /* vertices ids */
    };
    
    class Peice : public Face
    {
        public:
            Peice():Face(){};

            /* compute parameterized coefficient */
            void precomputation();
 
            /* geometry properties */
            /* reverse matrix of parameter matrix [ u, v1; 0, v2] 
             * u = |r1 - r3|, 
             * v1 = (r2 - r3) * (r1 - r3) / u
             * v2 = |(r2 - r3) - (r1 - r3) * v/u|
             * */
            Mat2 Bm;    
            double W;   /* Area of face, |(r1 - r3) x (r2 - r3)|/2 */

            /* topology properties */
            iVec3 hindge_id = iVec3(-1);  /* index of hindge id correspond to vertex */
    };

    class Hindge
    {
        public:
            Hindge(){};

            /* geometry properties */

            /* topology properties */
            iVec2 peice_info[2];        /* information of peices, x stands for peice id, y for v local id*/
    };

    class Tetra
    {   
        public:
            Tetra(){};

            /* compute volume W and reverse matrix Bm */
            void precomputation();

            /* geometry property */
            Mat3 Bm;    /* reverse matrix of [r1-r4; r2-r4; r3-r4] */
            double W;   /* volume of tetrahedron */

            /* topology property */
            pVertex v[4];
            iVec4 v_id;
    };

    //////////////////////////////////////////////////// organized structures
    class Film
    {
        public:
            Film(){};
            
            /* compute hindges based on given peices, fill in hindges and correspond
             * part of peices elements*/
            void computeHindges();

            std::vector<Peice> peices;
            std::vector<Hindge> hindges;
    };

	class Rigid
	{
        public:
        Rigid(){};
        Rigid(const std::vector<size_t> &vertex_ids)
        {
            elements = vertex_ids;
        };

        /* topology property */
        std::vector<size_t> elements; 
        
        /* geometry property */
        Vec3 m_cord;
        Vec3 m_pos;
        Quat m_rot = Quat(1, 0, 0, 0);
        Vec3 m_omega = Vec3(0); 
	};

    class Hole
    {
        public:
            Hole(){};
            
            std::vector<size_t> vertices; /* related vertices */
            std::vector<Face> holeface;   /* related surfaces */
    };

    ///////////////////////////////////////////////////// total container
    class TetraMesh
    {
        public:
            TetraMesh(){};

            int read( const std::string& filename );
			int write(const std::string& filename, const std::string& additioninfo);

            /* perform after read to compute surface and tetra properties */
            void precomputation();

			/* recompute surface normal when position changed */
			void recomputeSurfaceNorm();

			/* re-exame vertices label to check if they are fixed */
			void labelFixedId();

            /* add hole information, by input iVec3 vector */
            int addFilm( const std::vector<iVec3>& peice_ids );

            /* add rigid body constrains, by input v_ids contained in vector */
            int addRigidBody( const std::vector<size_t>& vertex_ids );

            /* add hole information, by input iVec3 vector */
            int addHole( const std::vector<iVec3>& face_ids );

		    /* variables */
            std::vector<Vertex> vertices;       /* vectors of vertices */
            std::vector<Tetra> tetrahedrons;    /* tetra mesh */
            std::vector<Face>   surface;        /* boudary faces of tetramesh */
			std::vector<size_t> fixed_ids;      /* fixed ids */
            std::vector<Film>   films;          /* films */
            std::vector<Rigid>  rigids;   /* rigid parts */
            std::vector<Hole>   holes;          /* holes */
    };

}




#endif //_TETRAMESH_H_
