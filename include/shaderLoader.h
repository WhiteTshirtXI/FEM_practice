#ifndef _SHADERLOADER_H_
#define _SHADERLOADER_H_

#include <fstream>
#include <string>
#include <vector>
#include <GL/glew.h>

using namespace std;

GLuint shaderLoader(const char* vertexShader, const char* fragmentShader)
{
    GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
    GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

    /* load source code */
    ifstream vertexFile(vertexShader, std::ios::in);
    string vertexCont((istreambuf_iterator<char>(vertexFile)), istreambuf_iterator<char>());

    ifstream fragmentFile(fragmentShader, std::ios::in);
    string fragmentCont((istreambuf_iterator<char>(fragmentFile)), istreambuf_iterator<char>());

    GLint Result = GL_FALSE;
    int InfoLogLength;
    
    /* compile vertex shader */
    printf("Compiling vertex shader : %s\n", vertexShader);
    char const * vertexContPointer = vertexCont.c_str();
    glShaderSource(VertexShaderID, 1, &vertexContPointer, NULL);
    glCompileShader(VertexShaderID);

    /* check vertex shader */
    glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    vector<char> VertexShaderErrorMessage(InfoLogLength);
    glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
    fprintf(stdout, "%s\n", &VertexShaderErrorMessage[0]);
    
    /* compile fragment shader */
    printf("Compiling fragment shader: %s\n", fragmentShader);
    char const * fragmentContPointer = fragmentCont.c_str();
    glShaderSource(FragmentShaderID, 1, &fragmentContPointer, NULL);
    glCompileShader(FragmentShaderID);

    /* check fragment shader*/
    glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    vector<char> FragmentShaderErrorMessage(InfoLogLength);
    glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
    fprintf(stdout, "%s\n", &FragmentShaderErrorMessage[0]);

    /* link program */
    fprintf(stdout, "Linking program\n");
    GLuint ProgramID = glCreateProgram();
    glAttachShader(ProgramID, VertexShaderID);
    glAttachShader(ProgramID, FragmentShaderID);
    glLinkProgram(ProgramID);

    /* check program*/
    glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
    glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    vector<char> ProgramErrorMessage(InfoLogLength);
    glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
    fprintf(stdout, "%s\n", &ProgramErrorMessage[0]);
    
    glDeleteShader(VertexShaderID);
    glDeleteShader(FragmentShaderID);

    return ProgramID;
}

#endif
