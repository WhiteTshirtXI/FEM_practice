#include <iostream>
#include <sstream>

//#include "../include/Dependence.h"
#include "../include/Types.h"
#include "../include/TetraMesh.h"
#include "../include/TetraMeshIO.h"
#include "../include/Dynamic.h"

#include <glm/glm.hpp>
using namespace glm;

#include "gtest/gtest.h"


using namespace BalloonFEM;

class BendingTest : public::testing::Test
{

    protected:
        virtual void SetUp()
        {
            /* build a squar of length sqrt(2) on xoy plane */
            std::stringstream in;
            in << "v 0 -1 0\n";   /* A */
            in << "v 1  0 0\n";   /* B */
            in << "v 0  1 0\n";   /* C */
            in << "v -1 0 0\n";   /* D */
            in << "fl 1 1 2 3\n"; /* ABC */
            in << "fl 1 1 3 4\n"; /* ACD */

            TetraMeshIO::read(in, m_tetra);
			m_tetra.precomputation();

            m_engine = new Engine(&m_tetra);
        }

        TetraMesh m_tetra;
        Engine* m_engine;
};

TEST_F(BendingTest, TestSetUp){
    EXPECT_EQ(4, m_tetra.num_vertex);
	EXPECT_EQ(2, m_tetra.num_pieces);
	EXPECT_EQ(1, m_tetra.num_hindges);
}


TEST_F(BendingTest, TestForce){

	const double PI = 3.1415926;
	const int DIV = 100;
	Vvec3 force;
	SpMat K;

	for (int i = 0; i < DIV; i++)
	{
		/* bend the B(1, 0, 0) to B'(cos(), 0, sin()) */
		m_tetra.vertices[1].m_pos = Vec3(cos(PI*i / DIV), 0, sin(PI*i / DIV));
		m_engine->inputData();

		K = m_engine->forceTest(force);

		EXPECT_GT(1e-6, length(force[3] - Vec3(0, 0, -sin(PI*i / DIV)))) << force[3].z << " " << -sin(PI*i / DIV) << std::endl;
	}

 //   /* bending force on A, C should point at (-1, 0, 1) */
	//Vec3 d = normalize(Vec3(-1, 0, 1));
	//EXPECT_GT(1e-4, length(normalize(force[0]) - d));
	//EXPECT_GT(1e-4, length(normalize(force[2]) - d));

 //   /* bending force on B should parallel to (1, 0, 0) */
	//d = Vec3(1, 0, 0);
	//EXPECT_GT(1e-4, length(normalize(force[1]) - d));

 //   /* bending force on D should parallel to (0, 0, -1) */
	//d = Vec3(0, 0, -1);
	//EXPECT_GT(1e-4, length(normalize(force[3]) - d));
}

TEST_F(BendingTest, TestStiffness){
	/* bend the B(1, 0, 0) to B'(-1/sqrt(2), 0, 1) */
	m_tetra.vertices[1].m_pos = Vec3(-sqrt(2)/2, 0, sqrt(2)/2);
	m_engine->inputData();

	Vvec3 force;
	SpMat K = m_engine->forceTest(force);

	/* exame the force */
    SpVec dp = SpVec::Zero(m_tetra.num_vertex * 3);
    dp( 3 * 1 + 0 ) = -1;       /* move B'(0, 0, 1) towards -x */

    SpVec df = K * dp;

	std::cout << df << std::endl;

    /* exame */
    /* D should be more compressed */
	EXPECT_EQ(0, df(3 * 3));
	EXPECT_EQ(0, df(3 * 3 + 1));
	EXPECT_GT(0, df(3 * 3 + 2)) << df << std::endl ;

	/* A, C should be symmetry by xoz plane */
	EXPECT_EQ(df(3 * 0), df(3 * 2));
	EXPECT_EQ(df(3 * 0 + 1), -df(3 * 2 + 1));
	EXPECT_EQ(df(3 * 0 + 2), df(3 * 2 + 2));
}
