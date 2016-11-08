#include "windows.h"
#include <cstdio>

//#pragma comment(lib, "glfw3.lib")
//#pragma comment(lib, "gtest.lib")


#include "gtest/gtest.h"

#include "../include/Viewer.h"

View::Viewer *p_viewer = NULL;

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
	RUN_ALL_TESTS();
	std::printf("Press Any Key To Exit.\n");
	std::getchar();
    return 0;
}
