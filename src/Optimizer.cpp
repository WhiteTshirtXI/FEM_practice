#include <iostream>
#include "Types.h"
#include "Viewer.h"
#include "Optimizer.h"

namespace
{
    using namespace BalloonFEM;

    void pushMat3x2(std::vector<T> &coeff, iVec3 v_id, int f_id, Mat3x2 H)
    {
        Mat3 F(H[0], H[1], -H[0]-H[1]);
        for (size_t i = 0; i < 3; i++)
            for(size_t j = 0; j < 3; j++)
            {
                coeff.push_back( T( 3 * v_id[i] + j, f_id, F[i][j]) );
            }
    }

    SpVec vvec3TospVec(Vvec3 &f)
    {
        SpVec u = SpVec::Zero(3 * f.size());
        for(size_t i = 0; i < f.size(); i++)
        {
            u( 3 * i     ) = f[i].x;
            u( 3 * i + 1 ) = f[i].y;
            u( 3 * i + 2 ) = f[i].z;
        }

        return u;
    }

}

extern View::Viewer *p_viewer;

namespace BalloonFEM
{
    ////////////////////////////// OptState Function ///////////////////////////
	size_t OptState::freedomDegree(){
		return 3 * m_size + 6 * m_r_size + m_tetra->num_pieces + m_tetra->holes.size();
	};

    void OptState::update(SpVec dpos)
    {
        size_t offset = 0;

		/* update vertex position */
        for (size_t i = 0; i < m_size; i++)
        {
            Vec3 &vpos = this->m_pos[i];
            vpos.x += dpos( 3 * i     + offset );
            vpos.y += dpos( 3 * i + 1 + offset );
            vpos.z += dpos( 3 * i + 2 + offset );
        }
        offset += 3 * m_size;

		/* update rigid body */
        for (size_t i = 0; i < m_r_size; i++)
        {
            Vec3 &rpos = this->m_r_pos[i];
            rpos.x += dpos( 6 * i     + offset );
            rpos.y += dpos( 6 * i + 1 + offset );
            rpos.z += dpos( 6 * i + 2 + offset );
            Quat &rrot = this->m_r_rot[i];
            Vec3 drot(0);
            drot.x = dpos( 6 * i + 3 + offset );
            drot.y = dpos( 6 * i + 4 + offset );
            drot.z = dpos( 6 * i + 5 + offset );
            rrot = omegaToQuat(drot) * rrot; 
        }
        offset += 6 * m_r_size;

        /* update thickness */
        thickness += dpos.segment(offset, m_tetra->num_pieces);
        offset += m_tetra->num_pieces;

        /* update pressure */
        pressure += dpos.tail(m_tetra->holes.size()); 

        this->project();
    }

    /////////////////////////////// Optimizer Function ////////////////////////

	Optimizer::Optimizer(
		TetraMesh* tetra,
		TetraMesh* target,
		ElasticModel* volume_model,
		AirModel* air_model,
		FilmModel* film_model,
		BendingModel* bend_model
		) :Engine(tetra, volume_model, air_model, film_model, bend_model)
	{
        /* copy from Engine, but change state to OptState */
        m_tetra = tetra;
        m_size = tetra->vertices.size();
        m_volume_model = volume_model;
        m_air_model = air_model;
		m_film_model = film_model;
        m_bend_model = bend_model;

		f_ext.assign(m_size, Vec3(0));

		cur_state = new OptState();
		next_state = new OptState();

		this->inputData();

		/* read in target */
		m_target = target;
		target_state = new OptState();
		target_state->input(m_target);
	};

	#define CONVERGE_ERROR_RATE 1e-4
    void Optimizer::solveOptimal()
    {
		/* initialize next_state */
		*next_state = *cur_state;
		next_state->project();
		SpVec f_sum;

		/* initialize temp variable for iterative implicit solving */
		/* f is total force on each vertex */
		/* K = - df/dr, here Force Diff Mat compute df/dr */
		SpMat K;
		computeForceAndGradient(*next_state, *target_state, f_sum, K);

		SpVec dstate = SpVec::Zero(next_state->freedomDegree());

		/* solver */
		Eigen::SimplicialLDLT<SpMat> solver;
		SpVec &b = f_sum;

		double err_f = f_sum.dot(f_sum);
		double err_begin = err_f;
		int count_iter = 0;		/* K dx = f iter count */
		/* while not converge f == 0, iterate */
		while ((err_f > CONVERGE_ERROR_RATE * err_begin) && (err_f > 1e-10))
		{
			/* debug use */
			count_iter++;
			printf("%d iter of K dv = f , err_felas = %.4e \n", count_iter, err_f);

			/* r0 = b - Ax0 */
			SpVec r = b - K * dstate;

			printf("building solver\n");
			solver.compute(K);
			if (solver.info() != Eigen::Success)
			{
				printf("decomposition failed!\n");
				printf("Number of non zeros: %d \n", K.nonZeros());
				return;
			}
			printf("solve delta_x \n");
			SpVec dstate = solver.solve(r);

			/* update v_pos_next and f_sum */
			((OptState*)next_state)->update(dstate);
			dstate.setZero();

			/* debug watch use*/
			next_state->output();
			p_viewer->refresh(1);
			std::cout << "pressure: " << next_state->pressure << std::endl;
			std::cout << "thickness: " << next_state->thickness << std::endl;

			/* update K and f*/
			computeForceAndGradient(*next_state, *target_state, f_sum, K);

			err_f = f_sum.dot(f_sum);
		}

		printf("f_sum error %f \n", err_f);
		printf("finish solving \n");
    }

    void Optimizer::computeForceAndGradient(ObjState &state, ObjState &target, SpVec &f, SpMat &A)
    {
        Vvec3 f_sum;
        f_sum.assign( m_size, Vec3(0.0));

        SpMat Tri( 3 * m_tetra->num_vertex, m_tetra->num_pieces );

        /* compute elastic forces by tetrahedrons */
        computeElasticForces(state, f_sum); 

		/* compute film forces by pieces */
		computeFilmForces(state, f_sum, Tri);
        
        /* compute forces by air pressure */
		computeAirForces(state, f_sum);

        for(size_t i = 0; i < m_size; i++)
            f_sum[i] += f_ext[i];

		/* compute forces diff by air pressure */
        SpMat K = computeAirDiffMat(state);

		/* compute film forces by pieces */
		K += computeFilmDiffMat(state);

		/* compute elastic forces diff by tetrahedrons */
		K += computeElasticDiffMat(state);

        /* compute bending force and gradient */
        K -= bendingForceAndGradient(state, f_sum);
        
        /* convert force to SpVec */
        SpVec f_real = vvec3TospVec( f_sum );
        SpVec x = vvec3TospVec( state.world_space_pos ) - vvec3TospVec(target.world_space_pos);
        SpVec h = state.thickness - target.thickness;
        SpVec press = vvec3TospVec( state.volume_gradient );

        /* tmp mat */
        size_t freedegree = state.freedomDegree();
        size_t kineticDegree = 3 * m_tetra->num_vertex + 6 * m_tetra->rigids.size();
        SpMat mat_a( 3 * m_tetra->num_vertex, freedegree);		/* target position displacement mat */
        SpMat mat_b( m_tetra->num_pieces, freedegree );			/* target thickness displacement mat */
        SpMat mat_c( 3 * m_tetra->num_vertex, freedegree);		/* total force intensity mat */
		SpMat mat_d(3 * m_tetra->num_vertex, freedegree);		/* restrict mat */

        SpMat I(m_tetra->num_pieces, m_tetra->num_pieces);
        I.setIdentity();

        mat_a.leftCols(kineticDegree) = state.projectMat();
        mat_b.middleCols(kineticDegree, m_tetra->num_pieces) = I;
        
        mat_c.leftCols(kineticDegree) = K * state.projectMat();
        mat_c.middleCols(kineticDegree, m_tetra->num_pieces) = Tri;
        mat_c.rightCols(1) = press.sparseView();
		mat_c = state.projectMat().transpose() * mat_c;

		mat_d.leftCols(kineticDegree) = state.restrictedMat();

        f = m_alpha * mat_a.transpose() * x 
            + m_beta * mat_b.transpose() * h 
			+ m_gamma * mat_c.transpose() * state.projectMat().transpose() * f_real;
		f = -f;

		A = m_alpha * mat_a.transpose() * mat_a
			+ m_beta * mat_b.transpose() * mat_b
			+ m_gamma * mat_c.transpose() * mat_c
			+ mat_d.transpose() * mat_d;

    }

	void Optimizer::computeFilmForces(ObjState &state, Vvec3 &f_sum, SpMat& Tri)
    {
        Vvec3 &pos = state.world_space_pos;

        Tri.setZero(); /* size should be (3*m_size, m_tetra->num_pieces) */
        std::vector<T> triangle_coeff;
        triangle_coeff.reserve( 9 * m_tetra->num_pieces );

		/* compute film elastic force */
        int piece_id = 0;
		for (MIter f = m_tetra->films.begin(); f != m_tetra->films.end(); f++)
		{
			for (PIter p = f->pieces.begin(); p != f->pieces.end(); p++)
			{
				iVec3 &id = p->v_id;
				Vec3 &v0 = pos[id[0]];
				Vec3 &v1 = pos[id[1]];
				Vec3 &v2 = pos[id[2]];

				/* calculate deformation in world space */
				Mat3x2 Ds = Mat3x2(v0 - v2, v1 - v2);

				/* calculate deformation gradient */
				Mat3x2 F = Ds * p->Bm;

				/* calculate Piola for this tetra */
				Mat3x2 P = m_film_model->Piola(F);

				/* calculate forces contributed from this tetra */
				Mat3x2 H = - p->W * P * transpose(p->Bm);

                pushMat3x2(triangle_coeff, id, piece_id, H);

                H *= state.thickness(piece_id);

				f_sum[id[0]] += H[0];
				f_sum[id[1]] += H[1];
				f_sum[id[2]] -= H[0] + H[1];

                piece_id ++;
			}
		}

        Tri.setFromTriplets(triangle_coeff.begin(), triangle_coeff.end());
    }

	void Optimizer::testFunc()
	{
		cur_state->thickness *= 2;
	}
}
