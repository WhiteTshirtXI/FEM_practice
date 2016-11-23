#include <iostream>

#include <glm/glm.hpp>
#include <glm/gtx/euler_angles.hpp>

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
		return 3 * m_size + 6 * m_r_size + 2 * m_tetra->num_pieces;
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
        // thickness += dpos.segment(offset, m_tetra->num_pieces);
        // offset += m_tetra->num_pieces;
    
        /* update aniso_sigma */
        aniso_sigma += dpos.segment(offset, 2 * m_tetra->num_pieces);
        offset += 2 * m_tetra->num_pieces;

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

        computeThicknessLap();
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
		double energy = computeGradientAndHessian(*next_state, *target_state, f_sum, K);

		SpVec dstate = SpVec::Zero(next_state->freedomDegree());

		/* solver */
		Eigen::SimplicialLDLT<SpMat> solver;
		SpVec &b = f_sum;

		double err_f = f_sum.dot(f_sum);
		double err_begin = err_f;
		int count_iter = 0;		/* K dx = f iter count */
		/* while not converge f == 0, iterate */
		while ((err_f > CONVERGE_ERROR_RATE * err_begin) && (err_f > 1e-10) && (count_iter < 5))
		{
			/* debug use */
			count_iter++;
			printf("%d iter of K dv = f , err_felas = %.4e, energy = %.4e \n", count_iter, err_f, energy);

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

			
			/* control step length */
			int cut_count = 0;
			double energy_next = computeGradientAndHessian(*next_state, *target_state, f_sum, K);
			while (energy < energy_next && cut_count < 10)
			{
				printf("cutting down dstate by half.\n");
				dstate /= 2.0;
				((OptState*)next_state)->update( - dstate);
				energy_next = computeGradientAndHessian(*next_state, *target_state, f_sum, K);
				cut_count++;
			}
			energy = energy_next;


			/* debug watch use*/
			next_state->output();
			p_viewer->refresh(1);
			std::cout << "pressure: " << next_state->pressure << std::endl;
			std::cout << "thickness: max " << next_state->thickness.maxCoeff() << std::endl;
			std::cout << "thickness: min " << next_state->thickness.minCoeff() << std::endl;
			std::cout << "sigma: max " << next_state->aniso_sigma.maxCoeff() << std::endl;
			std::cout << "sigma: min " << next_state->aniso_sigma.minCoeff() << std::endl;

			dstate.setZero();
		}

		printf("f_sum error %.4e energy %.4e\n", err_f, energy);
		printf("finish solving \n");
    }

	double Optimizer::computeGradientAndHessian(ObjState &state, ObjState &target, SpVec &f, SpMat &A)
    {
        Vvec3 f_sum;
        f_sum.assign( m_size, Vec3(0.0));

        SpMat Tri( 3 * m_tetra->num_vertex, m_tetra->num_pieces );
        SpMat Sigma( 3 * m_tetra->num_vertex, 2 * m_tetra->num_pieces );

        /* compute elastic forces by tetrahedrons */
        computeElasticForces(state, f_sum); 

		/* compute film forces by pieces */
		computeFilmForces(state, f_sum, Tri, Sigma);
        
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
        //K -= bendingForceAndGradient(state, f_sum);
        
        /* convert force to SpVec */
		SpVec f_freedeg = state.projectMat().transpose() * vvec3TospVec(f_sum);
        SpVec x = vvec3TospVec( state.world_space_pos ) - vvec3TospVec(target.world_space_pos);
		SpVec h = state.thickness;
		SpVec sig_delt = state.aniso_sigma.array() - 1.0;
		SpVec h_delt = m_L * h;

		/* used to balance influence of scale */
		double norm_coeff = m_tetra->num_vertex;
		norm_coeff *= norm_coeff;

		double energy = m_alpha / norm_coeff * x.squaredNorm()
			+ m_beta / norm_coeff * sig_delt.squaredNorm()
			+ m_gamma * norm_coeff * f_freedeg.squaredNorm();

		/* output errors */
		printf("pos_error = %.4e, sig_error = %.4e, f_error = %.4e, energy = %.4e \n", 
			x.squaredNorm(), sig_delt.squaredNorm(), f_freedeg.squaredNorm(), energy);

        /* tmp mat */
        size_t freedegree = state.freedomDegree();
        size_t kineticDegree = 3 * m_tetra->num_vertex + 6 * m_tetra->rigids.size();
        SpMat mat_a( 3 * m_tetra->num_vertex, freedegree);		/* target position displacement mat */
        SpMat mat_b( 2 * m_tetra->num_pieces, freedegree );			/* target thickness displacement mat */
        SpMat mat_c( 3 * m_tetra->num_vertex, freedegree);		/* total force intensity mat */
		SpMat mat_d(3 * m_tetra->num_vertex, freedegree);		/* restrict mat */

        SpMat I(2 * m_tetra->num_pieces, 2 * m_tetra->num_pieces);
        I.setIdentity();

        mat_a.leftCols(kineticDegree) = state.projectMat();
        mat_b.middleCols(kineticDegree, 2 * m_tetra->num_pieces) = I;
        
        mat_c.leftCols(kineticDegree) = K * state.projectMat();
        mat_c.rightCols(2 * m_tetra->num_pieces) = Sigma;
		mat_c = state.projectMat().transpose() * mat_c;

		mat_d.leftCols(kineticDegree) = state.restrictedMat();

		f = m_alpha/norm_coeff * mat_a.transpose() * x
			+ m_beta/norm_coeff * mat_b.transpose() * sig_delt
			+ m_gamma * norm_coeff * mat_c.transpose() * f_freedeg;

		f = -f;

		A = m_alpha/norm_coeff * mat_a.transpose() * mat_a
			+ m_beta/norm_coeff * mat_b.transpose() * mat_b
			+ m_gamma*norm_coeff * mat_c.transpose() * mat_c
			+ mat_d.transpose() * mat_d;

		return energy;
    }

    #define CONVERGE_ERROR_RATE 1e-4
	void Optimizer::solveOptimalGN()
	{
		/* initialize next_state */
		*next_state = *cur_state;
		next_state->project();
		size_t res_size = 3 * m_tetra->num_vertex + next_state->freedomDegree();
		SpVec res = SpVec::Zero(res_size);

		/* initialize temp variable for iterative implicit solving */
		/* res is the residual vector of energy */
		/* J =  d res / d x, the Jacobian of residual to parameters to be optimized */
		SpMat J;
		double energy = computeResidualAndJacobian(*next_state, *target_state, res, J);

		/* solver */
		Eigen::SimplicialLDLT<SpMat> solver;

		double energy_begin = energy;
		int count_iter = 0;		/* K dx = f iter count */
		/* while not converge f == 0, iterate */
		while ((energy > CONVERGE_ERROR_RATE * energy_begin) && (energy > 1e-10) && (count_iter < 5))
		{
			/* debug use */
			count_iter++;
			printf("%d iter of K dv = f , energy = %.4e \n", count_iter, energy);

			/* r0 = b - Ax0 */
			SpMat A = J.transpose() * J;

			printf("building solver\n");
			solver.compute(A);
			if (solver.info() != Eigen::Success)
			{
				printf("decomposition failed!\n");
				printf("Number of non zeros: %d \n", A.nonZeros());
				return;
			}
			printf("solve delta_x \n");

			SpVec dstate = - solver.solve(J.transpose() * res);

			/* update v_pos_next and f_sum */
			((OptState*)next_state)->update(dstate);

			/* update K and f*/
			double energy_next = computeResidualAndJacobian(*next_state, *target_state, res, J);

			/* control step length */
			int cut_count = 0;
			while (energy < energy_next && cut_count < 10)
			{
				printf("cutting down dstate by half.\n");
				dstate /= 2.0;
				((OptState*)next_state)->update(- dstate);
				next_state->output();
				p_viewer->refresh(1);
				energy_next = computeResidualAndJacobian(*next_state, *target_state, res, J);
				cut_count++;
			}

			/* debug watch use*/
			next_state->output();
			p_viewer->refresh(1);
			std::cout << "pressure: " << next_state->pressure << std::endl;
			std::cout << "thickness: max " << next_state->thickness.maxCoeff() << std::endl;
			std::cout << "thickness: min " << next_state->thickness.minCoeff() << std::endl;
			std::cout << "sigma: max " << next_state->aniso_sigma.maxCoeff() << std::endl;
			std::cout << "sigma: min " << next_state->aniso_sigma.minCoeff() << std::endl;

			energy = energy_next;
		}

		printf("energy %.4e\n", energy);
		printf("finish solving \n");
	}

	double Optimizer::computeResidualAndJacobian(ObjState &state, ObjState &target, SpVec &f, SpMat &A)
	{
		Vvec3 f_sum;
		f_sum.assign(m_size, Vec3(0.0));

		SpMat Tri(3 * m_tetra->num_vertex, m_tetra->num_pieces);
		SpMat Sigma(3 * m_tetra->num_vertex, 2 * m_tetra->num_pieces);
		
		computeElasticForces(state, f_sum); /* compute elastic forces by tetrahedrons */
		computeFilmForces(state, f_sum, Tri, Sigma); /* compute film forces by pieces */
		computeAirForces(state, f_sum); /* compute forces by air pressure */

		for (size_t i = 0; i < m_size; i++)
			f_sum[i] += f_ext[i];
		
		SpMat K = computeAirDiffMat(state); /* compute forces diff by air pressure */
		K += computeFilmDiffMat(state); /* compute film forces by pieces */
		K += computeElasticDiffMat(state); /* compute elastic forces diff by tetrahedrons */
		//K -= bendingForceAndGradient(state, f_sum); /* compute bending force and gradient */

		/* convert force to SpVec */
		SpVec f_freedeg = state.projectMat().transpose() * vvec3TospVec(f_sum);
		SpVec x = vvec3TospVec(state.world_space_pos) - vvec3TospVec(target.world_space_pos);
		SpVec h = state.thickness;
		SpVec sig_delt = state.aniso_sigma.array() - 1.0;
		SpVec h_delt = m_L * h;

		/* used to balance influence of scale */
	    size_t freedegree = state.freedomDegree();
		size_t kineticDegree = state.kineticDegree();

		double norm_coeff = m_tetra->num_vertex;
		double a = sqrt(m_alpha) / norm_coeff;
		double b = sqrt(m_beta) / norm_coeff;
		double c = sqrt(m_gamma) * norm_coeff;

        /* compute residual */
		f.head(3 * m_tetra->num_vertex) = a * x;
        f.segment(3 * m_tetra->num_vertex, 2 * m_tetra->num_pieces) = b * sig_delt;
	    f.tail(kineticDegree) =  c * f_freedeg;

		double energy = f.squaredNorm();

		/* output errors */
		printf("pos_error = %.4e, sig_error = %.4e, f_error = %.4e, energy = %.4e \n",
			x.squaredNorm(), sig_delt.squaredNorm(), f_freedeg.squaredNorm(), energy);

		/* tmp mat */
		SpMat mat_a(3 * m_tetra->num_vertex, freedegree);		/* target position displacement mat */
		SpMat mat_b(2 * m_tetra->num_pieces, freedegree);	    /* sigma regulazer mat */
		SpMat mat_c(3 * m_tetra->num_vertex, freedegree);		/* total force intensity mat */
		SpMat mat_d(3 * m_tetra->num_vertex, freedegree);		/* restriction mat */

		SpMat I(2 * m_tetra->num_pieces, 2 * m_tetra->num_pieces);
		I.setIdentity();

		mat_a.leftCols(kineticDegree) = state.projectMat();
		mat_b.middleCols(kineticDegree, 2 * m_tetra->num_pieces) = I;

		mat_c.leftCols(kineticDegree) = K * state.projectMat();
		mat_c.rightCols(2 * m_tetra->num_pieces) = Sigma;
		mat_c = state.projectMat().transpose() * mat_c;

        /* add restriction to force */
		mat_d.leftCols(kineticDegree) = state.restrictedMat();
		mat_c += mat_d;

		size_t res_size = mat_a.innerSize() + mat_b.innerSize() + mat_c.innerSize();
		SpMat Jt = SpMat(freedegree, res_size);			/* transpose of Jacobian matrix */

        Jt.leftCols( 3 * m_tetra->num_vertex ) = a * mat_a.transpose();
        Jt.middleCols( 3 * m_tetra->num_vertex, 2 * m_tetra->num_pieces) = b * mat_b.transpose();
        Jt.rightCols( kineticDegree ) = c * mat_c.transpose();

        A = Jt.transpose();

		return energy;
	}

	void Optimizer::computeFilmForces(ObjState &state, Vvec3 &f_sum, SpMat& Thk, SpMat& Sigma)
    {
        Vvec3 &pos = state.world_space_pos;

        Thk.setZero(); /* size should be (3*m_size, m_tetra->num_pieces) */
        std::vector<T> triangle_coeff;
        triangle_coeff.reserve( 9 * m_tetra->num_pieces );
        Sigma.setZero();
        std::vector<T> sigma_coeff;
        sigma_coeff.reserve( 9 * 2 * m_tetra->num_pieces);

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

                /* coordinate transform */
                Mat2 R = p->Bm * glm::orientate2(p->aniso_angle);

				/* calculate deformation gradient */
				Mat3x2 F = Ds * R;

                Vec2 sigma(state.aniso_sigma(2 * piece_id), state.aniso_sigma(2 * piece_id + 1));

				/* calculate Piola for this tetra */
				Mat3x2 P = m_film_model->Piola(F, sigma[0], sigma[1]);

				/* calculate forces contributed from this tetra */
				Mat3x2 H = - p->W * P * transpose(R);

                pushMat3x2(triangle_coeff, id, piece_id, H);

                H *= state.thickness(piece_id);

				f_sum[id[0]] += H[0];
				f_sum[id[1]] += H[1];
				f_sum[id[2]] -= H[0] + H[1];

                Mat3x2 S = - p->volume() * m_film_model->getMu() * F;
                Mat3x2 S0 = Mat3x2(S[0], Vec3(0)) * transpose(R) * 2.0 * sigma[0];
                Mat3x2 S1 = Mat3x2(Vec3(0), S[1]) * transpose(R) * 2.0 * sigma[1];
                pushMat3x2(sigma_coeff, id, 2 * piece_id    , S0);
                pushMat3x2(sigma_coeff, id, 2 * piece_id + 1, S1);

                piece_id ++;
			}
		}

        Thk.setFromTriplets(triangle_coeff.begin(), triangle_coeff.end());
        Sigma.setFromTriplets(sigma_coeff.begin(), sigma_coeff.end());
    }

	void Optimizer::computeThicknessLap()
	{
        std::vector<T> coeff;
        coeff.reserve(4 * m_tetra->num_hindges);
        SpVec diag = SpVec::Zero(m_tetra->num_pieces);
        
        int offset = 0;
        for(MIter f = m_tetra->films.begin(); f != m_tetra->films.end(); f++)
        {
            for(EIter h = f->hindges.begin(); h != f->hindges.end(); h++)
            {
                int i = h->piece_info[0].x + offset;
                int j = h->piece_info[1].x + offset;

                coeff.push_back(T(i, i, 1));
                coeff.push_back(T(i, j, -1));
                coeff.push_back(T(j, j, 1));
                coeff.push_back(T(j, i, -1));

                diag(i) += 1;
                diag(j) += 1;
            }
            offset += f->pieces.size();
        }

        m_L = SpMat(m_tetra->num_pieces, m_tetra->num_pieces);
        m_L.setFromTriplets(coeff.begin(), coeff.end());

        for(int i = 0; i < m_tetra->num_pieces; i++)
            diag(i) = 1.0 / diag(i);

        m_L = diag.asDiagonal() * m_L;
	}

	void Optimizer::testFunc()
	{
		cur_state->thickness *= 2;
	}
}
