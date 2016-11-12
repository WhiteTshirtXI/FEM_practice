#include "Types.h"
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

namespace BalloonFEM
{
    void Optimizer::solveOptimal()
    {
        this->solveStaticPos();
    }

    void Optimizer::computeForceAndGradient(ObjState &state, SpVec &f, SpMat &A)
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
        SpVec x = vvec3TospVec( state.world_space_pos );
        SpVec &h = state.thickness;
        SpVec press = vvec3TospVec( state.volume_gradient );

        /* tmp mat */
        size_t freedegree = state.freedomDegree() + m_tetra->num_pieces + m_tetra->holes.size();
        SpMat mat_a( 3 * m_tetra->num_vertex, freedegree);
        SpMat mat_b( m_tetra->num_pieces, freedegree );
        SpMat mat_c( 3 * m_tetra->num_vertex, freedegree);
        SpMat I(m_tetra->num_pieces, m_tetra->num_pieces);
        I.setIdentity();

        mat_a.leftCols(state.freedomDegree()) = state.projectMat();
        mat_b.middleCols(state.freedomDegree(), m_tetra->num_pieces) = I;
        
        mat_c.leftCols(state.freedomDegree()) = K * state.projectMat();
        mat_c.middleCols(state.freedomDegree(), m_tetra->num_pieces) = Tri;
        mat_c.rightCols(1) = press;

        f = m_alpha * mat_a.transpose() * x 
            + m_beta * mat_b.transpose() * h 
            + m_gamma * mat_c.transpose() * f_real;

        A = m_alpha * mat_a.transpose() * mat_a 
            + m_beta * mat_b.transpose() * mat_b
            + m_gamma * mat_c.transpose() * mat_c;

        A.topLeftCorner(state.freedomDegree(), state.freedomDegree()) += state.restrictedMat();
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
}
