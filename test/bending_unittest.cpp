#include <iostream>
#include <sstream>

//#include "../include/Dependence.h"
#include "../include/Types.h"
#include "../include/TetraMesh.h"
#include "../include/TetraMeshIO.h"
#include "../include/Dynamic.h"

#include <glm/glm.hpp>

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


//TEST_F(BendingTest, TestForce){
//    /* bend the B(1, 0, 0) to B'(0, 0, 1) */
//    m_tetra.vertices[1].m_pos = Vec3(0, 0, 1);
//    m_engine.inputData();
//
//    Vvec3 force;
//    SpMat K = m_engine.forceTest(force);
//
//    /* bending force on A, C should be zero */
//    EXPECT_GT( 1e-4, glm::length(force[0]));
//    EXPECT_GT( 1e-4, glm::length(force[2]));
//
//    /* bending force on B should parallel to (1, 0, 0) */
//	EXPECT_EQ(1, glm::dot(glm::normalize(force[1]), Vec3(1, 0, 0)));
//
//    /* bending force on D should parallel to (0, 0, -1) */
//	EXPECT_EQ(1, glm::dot(glm::normalize(force[3]), Vec3(0, 0, -1)));
//}

