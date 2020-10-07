#version 400

#ifdef VERTEX_SHADER

// model attributes
in vec3 a_position;
in vec3 a_normal;

// instance attributes
in vec3  i_translation;
in float i_height;

// uniforms
uniform mat4 gl_ModelViewMatrix;
uniform mat4 gl_ProjectionMatrix;
uniform vec2 u_cellSize;
uniform ivec2 u_gridCells;

// output
out vec3  worldPos;
out vec2  cellCenter;
out vec3  eyePos;
out vec3  worldNormal;
out float cubeHeight;
out ivec2 cellCoords;

void main(void)
{
    mat4 MVP    = gl_ProjectionMatrix * gl_ModelViewMatrix;
    
    worldPos    = a_position*vec3(u_cellSize, i_height) + i_translation;
    eyePos      = vec3(gl_ModelViewMatrix*vec4(worldPos, 1.0));
    worldNormal = a_normal;
    cubeHeight  = i_height;
    
    cellCoords = ivec2(gl_InstanceID/u_gridCells.y, gl_InstanceID%u_gridCells.y);
    cellCenter = (vec2(cellCoords) + 0.5) * u_cellSize;
    
    gl_Position = MVP * vec4(worldPos, 1.0); 
} 
#endif

#ifdef FRAGMENT_SHADER

#define PI 3.1415926538

// input from vertex shader
in vec3  worldPos;
in vec2  cellCenter;
in vec3  eyePos;
in vec3  worldNormal;
in float cubeHeight;
in ivec2 cellCoords;

// uniforms
uniform vec4  u_materialAmbient;
uniform vec4  u_materialDiffuse;
uniform vec4  u_materialSpecular;
uniform float u_materialShininess;
uniform vec3  u_lightPos; 

uniform sampler2D u_texture;
uniform sampler2D u_arrowTexture;

uniform int  u_useTexture; // 0: ambient material color, 1,2: texture, 3: custom shader defined
uniform bool u_shadeCube;
uniform vec2 u_worldSize;
uniform vec2 u_cellSize;

// raw data for custom shader texturing
uniform sampler2D u_datatexBed;
uniform sampler2D u_datatexIce;
uniform vec3 u_bedRange; // min, max, range
uniform vec3 u_iceRange;

// edit cursor
uniform bool  u_drawCursor;
uniform vec3  u_cursorPoint;
uniform float u_cursorRadius;

// output
out vec4 fragColor;

vec4 phongModel(vec4 c_amb, vec4 c_diff, vec4 c_spec) {
    
    // Returned color set to ambient
    vec3 c = c_amb.rgb;
  
    // Modified diffuse lighting
    vec3 L = normalize(u_lightPos - worldPos);
    vec3 N = normalize(worldNormal);
    float d = 0.5*(1.0 + dot(N, L));
    c += c_diff.rgb * vec3(d*d);

    // Specular 
    vec3  R = reflect(L, N);
    float l = 0.5*(1.0 + dot(R, eyePos));
    float s = pow(l*l, u_materialShininess);
    c += c_spec.rgb * vec3(s);
  
    return vec4(c, 1.0);
}


// All components are in the range [0…1], including hue.
vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

vec2 rotate(vec2 v, float a) {
    float s = sin(a);
    float c = cos(a);
    mat2 m = mat2(c, -s, s, c);
    return m * v;
}


void main()
{
    vec4 colorAmbient;
    
    // uniform color
    if (u_useTexture == 0) {
        colorAmbient = u_materialAmbient;
    }
    
    // shading using a texture with blocky appearance
    else if (u_useTexture == 1) {
        vec2 texCoords = cellCenter.xy/u_worldSize.xy;
        colorAmbient = texture(u_texture, texCoords);
    }
    
    // shading using a texture with interpolation
    else if (u_useTexture == 2) {
        vec2 texCoords = worldPos.xy/u_worldSize.xy;        
        colorAmbient = texture(u_texture, texCoords);
    }
    
    // custom shading with access to bed and ice layers, for debugging purposes
    else if (u_useTexture == 3) {
        
        vec2 texCoords = cellCenter.xy/u_worldSize.xy;
        vec2 texelSize = u_cellSize/u_worldSize.xy;
        
        float h = u_iceRange.x + u_iceRange.z*texture(u_datatexIce, texCoords).r;
        float b = u_bedRange.x + u_bedRange.z*texture(u_datatexBed, texCoords).r;
        float s = h + b;
        float hxm  = u_iceRange.x + u_iceRange.z*texture(u_datatexIce, texCoords + vec2(-1, 0)*texelSize).r;
        float hxp  = u_iceRange.x + u_iceRange.z*texture(u_datatexIce, texCoords + vec2( 1, 0)*texelSize).r;
        float hym  = u_iceRange.x + u_iceRange.z*texture(u_datatexIce, texCoords + vec2( 0 -1)*texelSize).r;
        float hyp  = u_iceRange.x + u_iceRange.z*texture(u_datatexIce, texCoords + vec2( 0, 1)*texelSize).r;
        float bxm  = u_bedRange.x + u_bedRange.z*texture(u_datatexBed, texCoords + vec2(-1, 0)*texelSize).r;
        float bxp  = u_bedRange.x + u_bedRange.z*texture(u_datatexBed, texCoords + vec2( 1, 0)*texelSize).r;
        float bym  = u_bedRange.x + u_bedRange.z*texture(u_datatexBed, texCoords + vec2( 0 -1)*texelSize).r;
        float byp  = u_bedRange.x + u_bedRange.z*texture(u_datatexBed, texCoords + vec2( 0, 1)*texelSize).r;
        float sxm  = bxm + hxm;
        float sxp  = bxp + hxp;
        float sym  = bym + hym;
        float syp  = byp + hyp;
        
        float grad_s_x = 0.5*(sxp - sxm)/u_cellSize.x;
        float grad_s_y = 0.5*(syp - sym)/u_cellSize.y;
        float grad_s  = sqrt(grad_s_x*grad_s_x + grad_s_y*grad_s_y);
        vec3 surfNormal = vec3(-grad_s_x, -grad_s_y, 1.0);
        float slopeAngle = 0.5*PI - acos(length(vec3(grad_s_x, grad_s_y, 0))/length(surfNormal));
        vec2 flowDir = -normalize(vec2(grad_s_x, grad_s_y));
        float flowDirAngle = atan(flowDir.y, flowDir.x);
        
        vec3 c;
        float ct = PI/8.0;
        if (slopeAngle <= ct)
            c = mix(vec3( 97.0,130.0,234.0), vec3(221.0,220.0,219.0), slopeAngle/ct);
        else
            c = mix(vec3(221.0,220.0,219.0), vec3(220.0, 94.0, 75.0), clamp(slopeAngle/ct - 1.0, 0.0, 1.0));
        colorAmbient = vec4(c/255.0, 1.0);
        
        // draw arrow
        vec2 arrowCoords = 1.6*(fract(worldPos.xy/u_cellSize.xy) - 0.5);
        arrowCoords = rotate(arrowCoords, flowDirAngle);
        vec4 arrowTex = texture(u_arrowTexture, 0.5 + arrowCoords);
        float aa = pow(arrowTex.a, 5);
        colorAmbient = aa*arrowTex + (1 - aa)*colorAmbient;
    }
    
    vec3 color;
    if (u_shadeCube)
        color = phongModel(colorAmbient, u_materialDiffuse, u_materialSpecular).rgb;
    else
        color = colorAmbient.rgb;
    
    if (u_drawCursor) {
        float t = length(worldPos.xy - u_cursorPoint.xy)/u_cursorRadius;
        t = 0.5*(1.0 - smoothstep(0.5, 1.0, t));
        color = color*(1 - t) + vec3(0.4,0.8,0.0)*t;
    }
    
    float alpha = min(cubeHeight/1.0, 1.0);
    
    fragColor = vec4(color, alpha);
}
#endif
