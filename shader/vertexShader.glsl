#version 430 core
layout (location = 0) in vec3 vPos;
layout (location = 1) in vec3 vNorm;
layout (location = 2) in vec3 vColor;

out vec3 Position_worldspace;
out vec3 Normal_worldspace;
out vec3 CameraDirection_worldspace;
out vec3 SurfaceColor;

uniform mat4 MVP;
uniform mat4 M;
uniform vec3 CameraPosition_worldspace;

void main(){
    /* vertex location after projection */
    gl_Position = MVP * vec4(vPos, 1.0);

    /* worldspace position and normal */
    Position_worldspace = (M * vec4(vPos, 1.0)).xyz;
    Normal_worldspace = (M * vec4(vNorm, 1.0)).xyz;

    /* worldspace vec3 from vertex to camera */
    CameraDirection_worldspace = CameraPosition_worldspace - Position_worldspace;

    /* vertex color */
    SurfaceColor = vColor;
}

