#version 430 core

#ifdef VERTEX_SHADER
void main(void)
{
    vec4 vertices[4] = vec4[4](vec4(-1.0, -1.0, 1.0, 1.0),
                               vec4( 1.0, -1.0, 1.0, 1.0),
                               vec4(-1.0,  1.0, 1.0, 1.0),
                               vec4( 1.0,  1.0, 1.0, 1.0));
    vec4 pos = vertices[gl_VertexID];

    gl_Position = pos;
}
#endif

#ifdef FRAGMENT_SHADER
layout (depth_greater) out float gl_FragDepth;

layout(location = 0) uniform vec3 CamPos;
layout(location = 1) uniform vec3 CamLookAt;
layout(location = 2) uniform vec3 CamUp;
layout(location = 3) uniform vec2 iResolution;


out vec4 color;

vec3 BuildRd()
{
	vec3 ro = CamPos;
	vec3 ta = CamLookAt;
	vec3 camUp  = CamUp;

	// Calculate orthonormal camera reference system
	vec3 camDir   = normalize(ta-ro); // direction for center ray
	vec3 camRight = normalize(cross(camDir,camUp));

	vec2 coord =-1.0+2.0*gl_FragCoord.xy/iResolution.xy;
	coord.x *= iResolution.x/iResolution.y;

	// Get direction for this pixel
	vec3 rd = normalize(camDir + (coord.x*camRight + coord.y*camUp)) ;

	return rd;
}

// Compute sky color 
// d  Ray direction
vec3 SkyShadeBlue(in vec3 d)
{
  	// light direction
	vec3 lig = normalize(vec3( 0.3,0.5, 0.6));
	float sun = 0.5*(dot(lig,d)+1.0);
	vec3 color = vec3(0.35,0.45,0.75)*(0.75-0.25*d.z);
	color += vec3(0.65,0.6,0.55)*pow( sun, 18.0 );
	return color;
}

void main(void)
{
	vec3 rd = BuildRd();

	color = vec4(SkyShadeBlue(rd),1.0);

	gl_FragDepth = 1.0;
}
#endif
