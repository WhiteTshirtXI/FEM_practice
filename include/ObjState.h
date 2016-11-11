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

		/* update from delata state represented by SpVec */
		void update(SpVec& dpos);

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

        /* the freedom of state, including those fixed vertices */
        size_t freedomDegree()
        {
            return 3 * m_size + 6 * m_r_size;
        }

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

        /* project Mat */
        SpMat m_W;

        /* restric mat */
        SpMat m_R;

		/* convert delta phi to quternion */
		Quat omegaToQuat(Vec3 delta_phi);

        /* calculate hole volume */
        void holeVolume();

        /* calculate hole volume and volume gradient */
        void volumeGradient();

        /* calculate projection mat */
        void computeProjectMat();

        /* the matrix specify those restricted vertices */
        void computeRestrictedMat();
	};
}

#endif //!_OBJSTATE_H_
