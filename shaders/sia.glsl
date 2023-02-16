#version 430 core
#extension GL_ARB_compute_shader : enable
#extension GL_ARB_shader_storage_buffer_object : enable

# ifdef COMPUTE_SHADER

uniform int GridSizeX;
uniform int GridSizeY;
uniform double dx;
uniform double dy; 
uniform double dt;

uniform double betaAccum;
uniform double betaAblate;
uniform double factorUdeform;
uniform double factorUslide;


layout(std430, binding=1) buffer Bedrock {
    double bedrock[];
};

layout(std430, binding=2) buffer InIce {
    double inIce[];
};

layout(std430, binding=3) buffer MapELA {
    double mapELA[];
};

layout(std430, binding=4) buffer MapBeta {
    double mapBeta[];
};

layout(std430, binding=5) buffer OutIce {
    double outIce[];
};

layout(std430, binding=6) buffer Diffusion {
    double diffusion[];
};


layout(local_size_x = WORK_GROUP_SIZE_X,  local_size_y = WORK_GROUP_SIZE_Y, local_size_z = 1) in;



// avoid divisions by zero
const double eps = 1e-6;

// Glenn law exponent
const double n = 3;

// These parameters are defined in section 7.1
const double A = 7.57e-17;//1e-16;
const double g = 9.81;
const double rho = 910.0;
const double rg = rho * g;
const double Gamma = 2.0 * A * rg*rg*rg / (n + 2);

const double Gamma_d = 7.26e-5;
const double Gamma_s = 3.27;



int GetOffset(int i, int j)
{
    return i * GridSizeY + j;
}

double H(int i, int j)
{
    if (i < 0) return 0;
    if (j < 0) return 0;
    if (i >= GridSizeX) return 0;
    if (j >= GridSizeY) return 0;
    return inIce[GetOffset(i, j)];
}

double B(int i, int j)
{
    if (i < 0) return 0;
    if (j < 0) return 0;
    if (i >= GridSizeX) return 0;
    if (j >= GridSizeY) return 0;
    return bedrock[GetOffset(i, j)];
}

double ELA(int i, int j)
{
    if (i < 0) return 0;
    if (j < 0) return 0;
    if (i >= GridSizeX) return 0;
    if (j >= GridSizeY) return 0;
    return mapELA[GetOffset(i, j)];
}

double BetaFactor(int i, int j)
{
    if (i < 0) return 0;
    if (j < 0) return 0;
    if (i >= GridSizeX) return 0;
    if (j >= GridSizeY) return 0;
    return mapBeta[GetOffset(i, j)];
}

double mdot(double z, double ela, double k)
{
    if (z >= ela) return betaAccum  * (z - ela) * k;
    else          return betaAblate * (z - ela);
}

double phi(double r)
{
    const double b = 2;
    return max(0, max(min(b*r, 1), min(r, b)));
}

double diffusivity(double grad_s, double h_p, double h_m, double s_p, double s_m)
{
    double D_p   = Gamma * h_p*h_p*h_p*h_p*h_p * grad_s;
    double D_m   = Gamma * h_m*h_m*h_m*h_m*h_m * grad_s;
    double D_min = min(D_p, D_m);
    double D_max = max(D_p, D_m);
    if (s_p <= s_m && h_m <= h_p) return D_min;
    if (s_p <= s_m && h_m >  h_p) return D_max;
    if (s_p >  s_m && h_m <= h_p) return D_max;
    if (s_p >  s_m && h_m >  h_p) return D_min;
    return 0;
}

double diffusivityWithSliding(double grad_s, double h_p, double h_m, double s_p, double s_m)
{
    double D_p   = h_p*h_p*h_p * (factorUdeform*Gamma_d*h_p*h_p + factorUslide*Gamma_s) * grad_s;
    double D_m   = h_m*h_m*h_m * (factorUdeform*Gamma_d*h_m*h_m + factorUslide*Gamma_s) * grad_s;
    double D_min = min(D_p, D_m);
    double D_max = max(D_p, D_m);
    if (s_p <= s_m && h_m <= h_p) return D_min;
    if (s_p <= s_m && h_m >  h_p) return D_max;
    if (s_p >  s_m && h_m <= h_p) return D_max;
    if (s_p >  s_m && h_m >  h_p) return D_min;
    return 0;
}

void main()
{
    int i = int(gl_GlobalInvocationID.x);
    int j = int(gl_GlobalInvocationID.y);

    if (i < 0) return;
    if (j < 0) return;
    if (i >= GridSizeX) return;
    if (j >= GridSizeY) return;
    
    double ela  = ELA(i, j);
    double beta = BetaFactor(i, j);

# ifdef SIA_DIR_AXES
    
    double dx2 = dx*dx;
    double dy2 = dy*dy;

    // access ice height positions in buffer
    double h     = H(i,j);
    double h_ip  = H(i+1,j);
    double h_ipp = H(i+2,j);
    double h_im  = H(i-1,j);
    double h_imm = H(i-2,j);
    double h_jp  = H(i,j+1);
    double h_jpp = H(i,j+2);
    double h_jm  = H(i,j-1);
    double h_jmm = H(i,j-2);

    // access bedrock positions and compute surface
    double z    = B(i,j);
    double s    = h    + z;
    double s_ip = h_ip + B(i+1,j);
    double s_im = h_im + B(i-1,j);
    double s_jp = h_jp + B(i,j+1);
    double s_jm = h_jm + B(i,j-1);

# endif 
# ifdef SIA_DIR_DIAGONALS

    double dx2 = dx*dx + dy*dy;
    double dy2 = dx*dx + dy*dy;

    // access ice height positions in buffer
    double h     = H(i,j);
    double h_ip  = H(i+1,j+1);
    double h_ipp = H(i+2,j+2);
    double h_im  = H(i-1,j-1);
    double h_imm = H(i-2,j-2);
    double h_jp  = H(i-1,j+1);
    double h_jpp = H(i-2,j+2);
    double h_jm  = H(i+1,j-1);
    double h_jmm = H(i+2,j-2);
    
    // access bedrock positions and compute surface
    double z    = B(i,j);
    double s    = h    + z;
    double s_ip = h_ip + B(i+1,j+1);
    double s_im = h_im + B(i-1,j-1);
    double s_jp = h_jp + B(i-1,j+1);
    double s_jm = h_jm + B(i+1,j-1);
    
# endif
    
    // compute downstream to upstream ice thickness ratio
    double r_i  = (h    - h_im) /(h_ip  - h    + eps);
    double r_ip = (h_ip - h)    /(h_ipp - h_ip + eps);
    double r_im = (h_im - h_imm)/(h     - h_im + eps);
    
    double r_j  = (h    - h_jm) /(h_jp  - h    + eps);
    double r_jp = (h_jp - h)    /(h_jpp - h_jp + eps);
    double r_jm = (h_jm - h_jmm)/(h     - h_jm + eps);
    
    
    // ice thickness at cell boundary (staggered grid)
    // up = +1/2, dn = -1/2
    double h_iup_m = h    + 0.5 * phi(r_i)  * (h_ip  - h);
    double h_iup_p = h_ip - 0.5 * phi(r_ip) * (h_ipp - h_ip);
    double h_idn_m = h_im + 0.5 * phi(r_im) * (h     - h_im);
    double h_idn_p = h    - 0.5 * phi(r_i)  * (h_ip  - h);
    
    double h_jup_m = h    + 0.5 * phi(r_j)  * (h_jp  - h);
    double h_jup_p = h_jp - 0.5 * phi(r_jp) * (h_jpp - h_jp);
    double h_jdn_m = h_jm + 0.5 * phi(r_jm) * (h     - h_jm);
    double h_jdn_p = h    - 0.5 * phi(r_j)  * (h_jp  - h);
    
    
    // slope gradients
    double grad_s_iup = (s_ip - s)*(s_ip - s)/dy2;
    double grad_s_idn = (s - s_im)*(s - s_im)/dy2;
    double grad_s_jup = (s_jp - s)*(s_jp - s)/dx2;
    double grad_s_jdn = (s - s_jm)*(s - s_jm)/dx2;
    
    
    // diffusivities at the 4 cell boundaries
    double D_iup = diffusivityWithSliding(grad_s_iup, h_iup_p, h_iup_m, s_ip, s);
    double D_idn = diffusivityWithSliding(grad_s_idn, h_idn_p, h_idn_m, s, s_im);
    double D_jup = diffusivityWithSliding(grad_s_jup, h_jup_p, h_jup_m, s_jp, s);
    double D_jdn = diffusivityWithSliding(grad_s_jdn, h_jdn_p, h_jdn_m, s, s_jm);
    
    
    // flux q divergence
    double div_q_i = (D_iup*(s_ip - s) - D_idn*(s - s_im))/dy2;
    double div_q_j = (D_jup*(s_jp - s) - D_jdn*(s - s_jm))/dx2;
    double div_q   = div_q_i + div_q_j;
    
    
    // mass balance
    double m = mdot(s, ela, beta);
    
    
    // explicit time stepping
    double h_res = h + dt*(m + div_q);
    h_res = max(h_res, 0);
    
    // max diffusion value, needed for adaptive timestep in next iteration
    double maxD = max(max(abs(D_iup), abs(D_idn)), max(abs(D_jup), abs(D_jdn)));
    // ice delta, used for checking glacier stability
    double dIce = h_res - h;
    
    
    // output
    int idx = GetOffset(i, j);
    outIce[idx] = h_res;
    diffusion[idx] = maxD;
    
}

#endif