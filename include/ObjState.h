
#include "Types.h"


namespace BallonFEM
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

		/* project to real world space and record in world_space_pos */
		void project();

		/* update from DeltaState */
		void update(DeltaState& dpos);

		friend class DeltaState;

	private:
		TetraMesh* m_tetra;

		/* number of vertices of tetramesh */
		size_t m_size;

		/* number of rigid bodys of tetramesh */
		size_t m_r_size;

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

	};

	class DeltaState
	{
	public:
		DeltaState(){};
		DeltaState(ObjState& other);

		typedef std::vector<Vec3> Vvec3;
		Vvec3 world_space_pos;

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