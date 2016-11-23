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
	void Optimizer::solveOptimal_GN()
	{
		/* solver */
		Eigen::SimplicialLDLT<SpMat> solver;

		/* initialize next_state */
		*next_state = *cur_state;
		next_state->project();

		/* res is the residual vector of energy */
		/* J =  d res / d x, the Jacobian of residual to parameters to be optimized */
		size_t res_size = 3 * m_tetra->num_vertex + next_state->freedomDegree();
		SpVec res = SpVec::Zero(res_size);
		SpMat J;
		double energy = computeResidualAndJacobian(*next_state, *target_state, res, J);

		double energy_begin = energy;
		int count_iter = 0;		/* \beta = - inv(JtJ) * Jt * res iter count */
		while ((energy > CONVERGE_ERROR_RATE * energy_begin) && (energy > 1e-10) && (count_iter < 5))
		{
			/* debug use */
			count_iter++;
			printf("%d iter of \beta = - inv(JtJ) * Jt * res , energy = %.4e \n", count_iter, energy);

			/* inv(JtJ) */
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

			/* update next state */
			((OptState*)next_state)->update(dstate);

			/* update J and res */
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
			energy = energy_next;

			/* debug watch use*/
			next_state->output();
			p_viewer->refresh(1);
			std::cout << "pressure: " << next_state->pressure << std::endl;
			std::cout << "thickness: max " << next_state->thickness.maxCoeff() << std::endl;
			std::cout << "thickness: min " << next_state->thickness.minCoeff() << std::endl;
			std::cout << "sigma: max " << next_state->aniso_sigma.maxCoeff() << std::endl;
			std::cout << "sigma: min " << next_state->aniso_sigma.minCoeff() << std::endl;

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
		SpVec dx = vvec3TospVec(state.world_space_pos) - vvec3TospVec(target.world_space_pos);
		SpVec h = state.thickness;
		SpVec dsig = state.aniso_sigma.array() - 1.0;
		SpVec dh = m_L * h;
        
        size_t opt_para_deg = state.freedomDegree();
		size_t kinetic_deg = state.kineticDegree();

		/* used to balance influence of scale */
		double norm_coeff = m_tetra->num_vertex;
		double a = sqrt(m_alpha) / norm_coeff;
		double b = sqrt(m_beta) / norm_coeff;
		double c = sqrt(m_gamma) * norm_coeff;

        /* compute residual */
		f.head(dx.size()) = a * dx;
        f.segment(dx.size(), dsig.size()) = b * dsig;
	    f.tail(f_kinetic.size()) =  c * f_kinetic;

		double energy = f.squaredNorm();

		/* output errors */
		printf("pos_error = %.4e, sig_error = %.4e, f_error = %.4e, energy = %.4e \n",
			dx.squaredNorm(), dsig.squaredNorm(), f_kinetic.squaredNorm(), energy);

		/* tmp mat */
		SpMat mat_a(dx.size(), opt_para_deg);		/* target position displacement mat */
		SpMat mat_b(dsig.size(), opt_para_deg);	    /* sigma regulazer mat */
		SpMat mat_c(3 * m_tetra->num_vertex, opt_para_deg);		/* total force intensity mat */
		SpMat mat_d(kinetic_deg, opt_para_deg);		/* restriction mat */

		SpMat I(dsig.size(), dsig.size());
		I.setIdentity();

		mat_a.leftCols(kinetic_deg) = state.projectMat();
		mat_b.middleCols(kinetic_deg, dsig.size()) = I;

		mat_c.leftCols(kinetic_deg) = K * state.projectMat();
		mat_c.rightCols(dsig.size()) = Sigma;
		mat_c = state.projectMat().transpose() * mat_c;

        /* add restriction to force */
		mat_d.leftCols(kinetic_deg) = state.restrictedMat().transpose() * state.restrictedMat();
		mat_c += mat_d;

		size_t res_size = dx.size() + dsig.size() + f_kinetic.size();
		SpMat Jt = SpMat(opt_para_deg, res_size);			/* transpose of Jacobian matrix */

        Jt.leftCols( dx.size() ) = a * mat_a.transpose();
        Jt.middleCols( dx.size() , dsig.size() ) = b * mat_b.transpose();
        Jt.rightCols( kinetic_deg ) = c * mat_c.transpose();

        A = Jt.transpose();

		return energy;
	}
}
