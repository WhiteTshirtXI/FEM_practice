cmake_minimum_required (VERSION 2.6)
project (tetra)

# original sources
include_directories( ${PROJECT_SOURCE_DIR}/include )
file( GLOB SRC_LIST  "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_executable(tetra ${SRC_LIST})

# required third party 
# OpenGL
find_package(OpenGL REQUIRED)
if (OPENGL_FOUND)
    include_directories(${OPENGL_INCLUDE_DIR})
    target_link_libraries(tetra ${OPENGL_gl_LIBRARY})
    target_link_libraries(tetra ${OPENGL_glu_LIBRARY})
    message(STATUS "OpenGL include dir: " ${OPENGL_INCLUDE_DIR})
endif (OPENGL_FOUND)

# GLEW
find_package(GLEW REQUIRED)
if (GLEW_FOUND)
    include_directories($(GLEW_INCLUCDE_DIRS))
    message(STATUS "GLEW include dir: " ${GLEW_INCLUCDE_DIRS})
endif (GLEW_FOUND)

# GLFW
find_package(glfw3 3.2 REQUIRED)
if (GLFW_FOUND)
    include_directories($(GLFW_INCLUDE_DIRS))
    message(STATUS "GLFW include dir: " ${GLFW_INCLUDE_DIRS})
endif (GLEW_FOUND)

