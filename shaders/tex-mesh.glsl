#version 430

#ifdef VERTEX_SHADER
in vec3 a_position;

uniform mat4 gl_ModelViewMatrix;
uniform mat4 gl_ProjectionMatrix;

out vec3 worldPos;

void main(void)
{
    mat4 MVP    = gl_ProjectionMatrix * gl_ModelViewMatrix;
    gl_Position = MVP * (vec4(a_position, 1.0));
    worldPos    = a_position;
} 
#endif

#ifdef FRAGMENT_SHADER
uniform vec2 u_worldMin;
uniform vec2 u_worldSize;
uniform sampler2D u_texture;

in vec3 worldPos;

out vec4 fragment;

void main()
{
    vec2 uv = (worldPos.xy - u_worldMin)/u_worldSize;
    vec4 c = texture(u_texture, uv);
    fragment = vec4(c.rgb, 1.0);
}
#endif
