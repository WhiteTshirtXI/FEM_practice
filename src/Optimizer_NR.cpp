#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Types.h"
#include "Viewer.h"
#include "Optimizer.h"

extern View::Viewer *p_viewer;

namespace{

	using namespace BalloonFEM;

	SpVec vvec3TospVec(Vvec3 &f)
	{
		SpVec u = SpVec::Zero(3 * f.size());
		for (size_t i = 0; i < f.size(); i++)
		{
			u(3 * i) = f[i].x;
			u(3 * i + 1) = f[i].y;
			u(3 * i + 2) = f[i].z;
		}

		return u;
	}
}

namespace BalloonFEM{

	#define CONVERGE_ERROR_RATE 1e-4
    void Optimizer::solveOptimal_NR()
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
		f_sum.assign(m_size, Vec3(0.0));

		SpMat Tri(3 * m_tetra->num_vertex, m_tetra->num_pieces);
		SpMat Sigma(3 * m_tetra->num_vertex, 2 * m_tetra->num_pieces);

		computeElasticForces(state, f_sum);          /* compute elastic forces by tetrahedrons */
		computeFilmForces(state, f_sum, Tri, Sigma); /* compute film forces by pieces */
		computeAirForces(state, f_sum);              /* compute forces by air pressure */

		for (size_t i = 0; i < m_size; i++) f_sum[i] += f_ext[i];   /* add external force */

		SpMat K = computeAirDiffMat(state);     /* compute forces diff by air pressure */
		K += computeFilmDiffMat(state);         /* compute film forces by pieces */
		K += computeElasticDiffMat(state);      /* compute elastic forces diff by tetrahedrons */
		//K -= bendingForceAndGradient(state, f_sum); /* compute bending force and gradient */
        
        /* convert force to SpVec */
		SpVec f_kinetic = state.projectMat().transpose() * vvec3TospVec(f_sum);
        SpVec dx = vvec3TospVec( state.world_space_pos ) - vvec3TospVec(target.world_space_pos);
		SpVec h = state.thickness;
		SpVec dsig = state.aniso_sigma.array() - 1.0;
		SpVec dh = m_L * h;

        size_t opt_para_deg = state.freedomDegree();
        size_t kinetic_deg = state.kineticDegree(); 
		
        /* used to balance influence of scale */
		double norm_coeff = m_tetra->num_vertex;
		norm_coeff *= norm_coeff;
		double a = m_alpha / norm_coeff;
		double b = m_beta / norm_coeff;
		double c = m_gamma * norm_coeff;

		double energy = a * dx.squaredNorm()
			+ b * dsig.squaredNorm()
			+ c * f_kinetic.squaredNorm();

		/* output errors */
		printf("pos_error = %.4e, sig_error = %.4e, f_error = %.4e, energy = %.4e \n", 
			dx.squaredNorm(), dsig.squaredNorm(), f_kinetic.squaredNorm(), energy);

        /* tmp mat */
      
        SpMat mat_a( dx.size(), opt_para_deg);		/* target position displacement mat */
        SpMat mat_b( dsig.size(), opt_para_deg );	/* target thickness displacement mat */
        SpMat mat_c( 3 * m_tetra->num_vertex, opt_para_deg);/* total force intensity mat */
		SpMat mat_d(kinetic_deg, opt_para_deg);		/* restrict mat */

    	SpMat I(dsig.size(), dsig.size());
		I.setIdentity();

        mat_a.leftCols(kinetic_deg) = state.projectMat();
        mat_b.middleCols(kinetic_deg, dsig.size()) = I;
        
        mat_c.leftCols(kinetic_deg) = K * state.projectMat();
        mat_c.rightCols(dsig.size()) = Sigma;
		mat_c = state.projectMat().transpose() * mat_c;

		mat_d.leftCols(kinetic_deg) = state.restrictedMat().transpose() * state.restrictedMat();

		f = a * mat_a.transpose() * dx
			+ b * mat_b.transpose() * dsig
			+ c * mat_c.transpose() * f_kinetic;

		f = -f;

		A = a * mat_a.transpose() * mat_a
			+ b * mat_b.transpose() * mat_b
			+ c * mat_c.transpose() * mat_c
			+ mat_d.transpose() * mat_d;

		return energy;
    }
}
