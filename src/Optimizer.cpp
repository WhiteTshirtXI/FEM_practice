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
