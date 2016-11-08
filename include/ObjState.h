#ifndef _OBJSTATE_H_
#define _OBJSTATE_H_

#include "Types.h"


namespace BalloonFEM
{
	/* ObjState fully represent the state of obj, and can be updated by
	* DeltaState */
	class ObjState;

	/* DeltaState is the infinitesimal displacement of ObjState, and can be
	* defined to perform linear algorithm
	* */
	class DeltaState;


	class ObjState
	{
	public:
		ObjState(){};
		ObjState(TetraMesh* tetra);
		ObjState(const ObjState& other);

		ObjState& operator=(const ObjState& other);

		void input(TetraMesh* tetra);
		void output();

		typedef std::vector<Vec3> Vvec3;
		Vvec3 world_space_pos;

        /* recorded volumes of holes */
        std::vector<double> hole_volume;

        /* volume gradient of vertices belongs to hole */
        Vvec3 volume_gradient;

		/* project to real world space and record in world_space_pos
         * then calcualte volume of holes and volume gradients
         */
		void project();

		/* update from DeltaState */
		void update(DeltaState& dpos);

        /* calculate volume gradient difference based on dr */
        void volumeGradientDiff(Vvec3& dr, Vvec3& dg);

        /* calculate volume gradient difference Matrix */
        SpMat volumeGradientDiffMat();

        /* calculate hindge angle gradient mat */
        SpMat hindgeAngleGradient();

		friend class DeltaState;

	private:
		TetraMesh* m_tetra;

		/* number of vertices of tetramesh */
		size_t m_size;

		/* number of rigid bodys of tetramesh */
		size_t m_r_size;

        /* numer of holes of tetramesh */
        size_t m_h_size;

		/* positions for those who is not fixed nor bounded to a rigid body
		* size equals to world_space_pos, meaningless nonzero value for
		* the fixed and bounded vertices.
		*/
		Vvec3 m_pos;

		/* positions for rigid body mass centers */
		Vvec3 m_r_pos;

		typedef std::vector<Quat> Vquat;
		/* rotation quaternion for rigid body */
		Vquat m_r_rot;

		/* convert delta phi to quternion */
		Quat omegaToQuat(Vec3 delta_phi);

        /* calculate hole volume */
        void holeVolume();

        /* calculate hole volume and volume gradient */
        void volumeGradient();


	};

	class DeltaState
	{
	public:
		DeltaState(){};
		DeltaState(ObjState& other);

		typedef std::vector<Vec3> Vvec3;
		Vvec3 world_space_pos;

		void clear();

		/* linear algebra algorithms */
		/* assign this DeltaState alpha * other */
		void assign(const double alpha, const DeltaState& other);
		/* add alpha * other to this DeltaState */
		void addup(const double alpha, const DeltaState& other);
		/* self-multiply alpha */
		void multiply(const double alpha);

		/* convert between obj state and real world spaces */
		/* project to real world space and record in world_space_pos */
		void project();

		/* project real world space data to constrained freedom state */
		void conterProject();

        /* the matrix enable inner freedom to project to real world */
        SpMat projectMat();

        /* the matrix specify those restricted vertices */
        SpMat restrictedMat();

        /* convert to SpVec, equal to world space vec multiply transpose projectMatrix() */
        SpVec toSpVec();

        /* read in SpVec and update the inner representation */
        void readSpVec(SpVec& stateVec);

		/* dot function */
		double dot(const DeltaState& other);

		friend class ObjState;

	private:
		TetraMesh* m_tetra;
		ObjState* m_state;

		/* number of vertices of tetramesh */
		size_t m_size;

		/* number of rigid bodys of tetramesh */
		size_t m_r_size;

		/* positions for those who is not fixed nor bounded
		* to a rigid body, size equals to world_space_pos,
		* but Vec3(0) value for the fixed and bounded vertices.
		*/
		Vvec3 dm_pos;

		/* positions for rigid body mass centers */
		Vvec3 dm_r_pos;

		/* rotation omega for rigid body */
		Vvec3 dm_r_rot;
	};
}

#endif //!_OBJSTATE_H_