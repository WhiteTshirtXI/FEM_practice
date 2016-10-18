#ifndef _DYNAMIC_H_
#define _DYNAMIC_H_

#include <vector>

#include "Types.h"
#include "TetraMesh.h"
#include "ElasticModel.h"

namespace BallonFEM
{
    /* ObjState fully represent the state of obj, and can be updated by 
     * DeltaState */
    class ObjState;
    
    /* DeltaState is the infinitesimal displacement of ObjState, and can be
     * defined to perform linear algorithm
     * */
    class DeltaState;
    
    /* take tetramesh and build ObjState and then solve them */
    class Engine;

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

            size_t m_size;

            typedef std::vector<Quat> Vquat;

            Vvec3 m_pos;
	};

    class DeltaState
    {
        public:
            DeltaState(){};
            DeltaState(const ObjState& other);

            typedef std::vector<Vec3> Vvec3;
            Vvec3 world_space_pos;

            /* linear algebra algorithms */
            /* assign this DeltaState alpha * other */
            void assign(const double alpha, const DeltaState& other);
            /* add alpha * other to this DeltaState */
            void addup( const double alpha, const DeltaState& other);
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

            size_t m_size;

            typedef std::vector<Quat> Vquat;

            Vvec3 m_pos;

	};

	double deltaStateDot(const DeltaState& a, const DeltaState& b);

    class Engine
    {
        public:
            Engine(){};
            Engine(TetraMesh* tetra, ElasticModel* model);

            void setModel(ElasticModel* model){ m_model = model; };

			/* input data v_pos, v_velocity and f_ext */
			void inputData();

            /* output v_pos and v_velocity data to tetra mesh */
            void outputData();

            /* copy v_pos_next, to v_pos*/
            void stepToNext();

            /* solve v_pos and v_velocity for next timestep */
            void solveNextTimestep(double timestep);

            /* solve v_pos for quasistatic simulation or static 
             * output static position to v_pos_next 
             */
            void solveStaticPos();

            /* test the comuteElasticForces, deforme v_pos and output f_ela 
             * to v_velocity
             * */
            void forceTest();

        private:
            TetraMesh* m_tetra;
            ElasticModel* m_model;
			size_t m_size;

            typedef std::vector<Vec3> Vvec3;

            /* current position and current velocity */
            ObjState cur_state;
            Vvec3 f_ext;

            /* constrains */
            std::vector<size_t> fixed_id;
            
            /* position and velocity for next time step, need to be solved
             * by backward Euler method
             * */
            ObjState next_state;

            /* compute elastic force when given positions pos*/
            void computeElasticForces(ObjState &pos, Vvec3 &f_elas);

            /* compute delta elastic force when given position pos and disturbe dpos */
            void computeForceDifferentials(ObjState &pos, DeltaState &dpos, Vvec3 &df_elas);

    };

}

#endif // _DYNAMIC_H_
