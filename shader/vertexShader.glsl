#version 430 core
layout (location = 0) in vec3 vPos;
layout (location = 1) in vec3 vNorm;

out vec3 Position_worldspace;
out vec3 Normal_worldspace;
out vec3 CameraDirection_worldspace;
out vec3 EyeDirection_worldspace;

uniform mat4 MVP;
uniform mat4 M;
uniform vec3 CameraPosition_worldspace;
uniform vec3 EyePosition_modelspace;

void main(){
    /* vertex location after projection */
    gl_Position = MVP * vec4(vPos, 1.0);

    /* worldspace position and normal */
    Position_worldspace = (M * vec4(vPos, 1.0)).xyz;
    Normal_worldspace = (M * vec4(vNorm, 1.0)).xyz;

    /* worldspace vec3 from vertex to camera */
    CameraDirection_worldspace = CameraPosition_worldspace - Position_worldspace;

    /* worldspace vec3 from vertex to eye */
    vec3 EyePosition_worldspace = (M * vec4(EyePosition_modelspace, 1.0)).xyz;
    EyeDirection_worldspace = EyePosition_worldspace - Position_worldspace;
}

