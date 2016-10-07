#version 430 core

in vec3 Position_worldspace;
in vec3 Normal_worldspace;
in vec3 CameraDirection_worldspace;
in vec3 EyeDirection_worldspace;

//uniform Light{
//    vec4 lightColor;
//    vec3 lightPos;
//    vec4 globalAmb;
//}

out vec4 color;

void main()
{
    vec3 l = normalize(CameraDirection_worldspace);
    vec3 n = normalize(Normal_worldspace);
    vec3 e = normalize(EyeDirection_worldspace);
    float cosAlpha = clamp(dot(n, l), 0, 1);
    float cosBeta = dot(n, e);

    vec3 amb = vec3(0.1);
    vec3 light = vec3(0.9);

    vec3 sight = amb + light * cosAlpha;
    color = vec4(sight, 0.4);

    if (abs(cosBeta)<7e-3){
        color = vec4(1.0, 0.0, 0.0, 1.0);
    }
}
