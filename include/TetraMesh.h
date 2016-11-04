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

            /* position in material coordinate */
            Vec3 m_cord;
            
            /* position in world space */
            Vec3 m_pos;

            /* velocity in world space */
            Vec3 m_velocity = Vec3(0);

            /* exterio force */
            Vec3 m_f_ext = Vec3(0);

            /* index in tetramesh vertex vector */
            size_t id;

            /* if is fixed */
            bool m_fixed = false;

            /* if is connected to a Rigid Body */
            int rigid = -1;

            /* if is connected to a hole */
            int hole = -1;
    };

    class Face
    {
        public:
            Face(){};
            
            /* compute normal */
            void precomputation();
            
            /* geometry property */
            Vec3 m_normal;

            /* topology property */
            pVertex v[3];
            iVec3 v_id;
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
            
            /* related vertices */
            std::vector<size_t> vertices;

            /* related surfaces */
            std::vector<Face> holeface;
    };

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

            /* add rigid body constrains, by input v_ids contained in vector */
            int addRigidBody( const std::vector<size_t>& vertex_ids );

            /* add hole information, by input iVec3 vector */
            int addHole( const std::vector<iVec3>& vertex_ids );

            /* vectors of vertices */
            std::vector<Vertex> vertices;

            /* boudary faces of tetramesh */
            std::vector<Face> surface;

            /* tetra mesh */
            std::vector<Tetra> tetrahedrons;

			/* fixed ids */
			std::vector<size_t> fixed_ids;

            /* rigid parts */
            std::vector<Rigid> rigid_bodies;

            /* holes */
            std::vector<Hole> holes;
    };

}




#endif //_TETRAMESH_H_
