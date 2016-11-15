#ifndef _OBJSTATE_H_
#define _OBJSTATE_H_

#include "Types.h"


namespace BalloonFEM
{
	/* ObjState fully represent the state of obj, and can be updated by
	* DeltaState */
	class ObjState
	{
	public:
		ObjState(){};
		ObjState(TetraMesh* tetra);
		ObjState(const ObjState& other);

		ObjState& operator=(const ObjState& other);

		virtual void input(TetraMesh* tetra);
		virtual void output();

		Vvec3 world_space_pos;

        /* recorded volumes of holes */
        std::vector<double> hole_volume;

        /* volume gradient of vertices belongs to hole */
        Vvec3 volume_gradient;

        SpVec thickness;	/* thickness of pieces */
		SpVec pressure;		/* pressure of each hole */

		/* project to real world space and record in world_space_pos
         * then calcualte volume of holes and volume gradients
         */
		void project();

		/* update from delata state represented by SpVec */
		virtual void update(SpVec& dpos);

		/* the freedom of state, including those fixed vertices */
		virtual size_t freedomDegree(){ return 3 * m_size + 6 * m_r_size; };

        /* calculate volume gradient difference based on dr */
        void volumeGradientDiff(Vvec3& dr, Vvec3& dg);

        /* calculate volume gradient difference Matrix */
        SpMat volumeGradientDiffMat();

        /* calculate hindge angle gradient mat */
        SpMat hindgeAngleGradient();

        /* projecti Mat, from unrestricted freedom to real world pos */
        SpMat projectMat(){return m_W;};

        /* the matrix specify those restricted vertices */
        SpMat restrictedMat(){return m_R;};

	protected:
		TetraMesh* m_tetra;

		
		size_t m_size;	/* number of vertices of tetramesh */
		size_t m_r_size; /* number of rigid bodys of tetramesh */
        size_t m_h_size; /* numer of holes of tetramesh */

		/* positions for those who is not fixed nor bounded to a rigid body
		* size equals to world_space_pos, meaningless nonzero value for
		* the fixed and bounded vertices.
		*/
		Vvec3 m_pos;

		Vvec3 m_r_pos;	/* positions for rigid body mass centers */
		Vquat m_r_rot;	/* rotation quaternion for rigid body */

		Quat omegaToQuat(Vec3 delta_phi); /* convert delta phi to quternion */

		void computeProjectMat();	 /* calculate projection mat */
		void computeRestrictedMat(); /* the matrix specify those restricted vertices */
        SpMat m_W;					 /* project Mat */
        SpMat m_R;					 /* restric mat */

        void holeVolume();		/* calculate hole volume */
        void volumeGradient();	 /* calculate hole volume and volume gradient */
	};
}

#endif //!_OBJSTATE_H_
