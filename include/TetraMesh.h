#ifndef _TETRAMESH_H_
#define _TETRAMESH_H_

#include <vector>
#include <string>
#include <Eigen/Dense>

#include "Types.h"

namespace BallonFEM{

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
            Rigid* rigid = NULL;
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
        Rigid(std::vector<size_t> vertex_ids)
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

    class TetraMesh
    {
        public:
            TetraMesh(){};

            int read( const std::string& filename );
            int write( const std::string& filename );

            /* perform after read to compute surface and tetra properties */
            void precomputation();

			/* recompute surface normal when position changed */
			void recomputeSurfaceNorm();

			/* re-exame vertices label to check if they are fixed */
			void labelFixedId();

            /* add rigid body constrains */
            int addRigidBody( const std::vector<size_t>& vertex_ids );

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
    };

}




#endif //_TETRAMESH_H_
