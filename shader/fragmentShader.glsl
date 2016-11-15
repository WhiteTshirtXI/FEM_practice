#version 430 core

in vec3 Position_worldspace;
in vec3 Normal_worldspace;
in vec3 CameraDirection_worldspace;
in vec3 SurfaceColor;

out vec4 color;

void main()
{
    vec3 l = normalize(CameraDirection_worldspace);
    vec3 n = normalize(Normal_worldspace);
    float cosAlpha = clamp(dot(n, l), 0, 1);

    vec3 surf = SurfaceColor * 0.3;
    vec3 light = vec3(0.7);

    vec3 sight = surf + light * cosAlpha;
    color = vec4(sight, 0.4);
}
