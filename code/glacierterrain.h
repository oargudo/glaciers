#ifndef GLACIERTERRAIN_H
#define GLACIERTERRAIN_H

#include "glew.h"
#include "core.h"
#include "scalarfield2.h"
#include <QtCore/QElapsedTimer>

class GlacierTerrain
{

public:
    GlacierTerrain();
    GlacierTerrain(int nx, int ny);
    GlacierTerrain(const ScalarField2& bedrock);
    GlacierTerrain(const ScalarField2& bedrock, const ScalarField2& ice);
    ~GlacierTerrain();

    void initSimulationGL();
    void configSimulation(const ScalarField2& mapELA, const ScalarField2& mapAccum, double accumRate, double ablateRate, double factorUdeform, double factorUslide);
    void resetSimulation();
    double runSimulationSteps(unsigned int numSimulSteps, double maxSimulTime = -1);

    void setIceMap(const ScalarField2& icemap);
    ScalarField2 remapIceSurface(const ScalarField2& hiresBed) const;


    // ScalarFields
    ScalarField2& GetBedrock() { return bedrock; }
    ScalarField2& GetIce() { return ice; }
    ScalarField2& GetELA() { return ela; }
    ScalarField2& GetAccumRate() { return beta; }
    ScalarField2 GetHeightfield() const { return bedrock + ice; }
    const ScalarField2& GetBedrock() const { return bedrock; }
    const ScalarField2& GetIce() const { return ice; }
    const ScalarField2& GetELA() const {return ela; }
    const ScalarField2& GetAccumRate() const { return beta; }
    const ScalarField2& GetDiffusivity() const { return diffusivity; }

    // Getters
    double Bedrock(const Vector2& p) const { return bedrock.value(p); }
    double Bedrock(int i, int j) const { return bedrock.at(i,j); }
    double Ice(const Vector2& p) const { return ice.value(p); }
    double Ice(int i, int j) const { return ice.at(i, j); }
    double Height(const Vector2& p) const { return bedrock.value(p) + ice.value(p); }
    double Height(int i, int j) const { return bedrock.at(i,j) + ice.at(i,j); }

    // Physical magnitudes
    double terrainArea() { return terrainSize[0]*terrainSize[1]; }
    double cellWidth() const { return terrainSize[0]/nx; }
    double cellHeight() const { return terrainSize[1]/ny; }
    double cellArea() const { return cellWidth()*cellHeight(); }
    bool isEmpty() const { return (nx <= 0) || (ny <= 0); }
    int numCellsX() const { return nx; }
    int numCellsY() const { return ny; }
	int vertexIndex(int i, int j) { return i + nx * j; }
    Box2 getDomain() const { return bedrock.getDomain(); }
    Box3 getBoundingBox() const {
        Box2 dom = getDomain();
        double hmin, hmax;
        bedrock.getRange(hmin, hmax);
        return Box3(Vector3(dom.getMin()[0], dom.getMin()[1], hmin),
                    Vector3(dom.getMax()[0], dom.getMax()[1], hmax));
    }

    double iceVolume() const;
    double elaValue(int, int) const;
    bool   aboveELA(int i, int j) const;
    bool   accumArea(int i, int j) const;
    double iceThickness(int, int) const;
    Vector2 iceGradient(int, int) const;
    double iceStress(int, int) const;
    double iceSpeed(int, int) const;
    double iceSpeedDeform(int, int) const;
    double iceSpeedSlide(int, int) const;
    Vector2 iceVelocity(int, int) const;
    double diffusionTerm(int, int) const;
    double diffusionTermRaw(int, int) const;
    double yearlyIce(bool withIce) const;
    double glaciatedArea() const;


    // Modifiers
    void SmoothRock(int = 1);
    void SmoothIce(int = 1);

    // Query
    bool Intersect(const Ray& ray, double& t, Vector3& q) const;
    QString GetStats() const;

    // GL Buffers
    GLuint bufferBedrockId() const { return glbufferBedrock; };
    GLuint bufferIceInId() const { return glbufferIceIn; };
    GLuint bufferIceOutId() const { return glbufferIceOut; };
    GLuint bufferDiffusionId() const { return glbufferD; };

    void updateBedrockGPU();
    void updateIceGPU();
    void updateELAGPU();
    void updateAccumRateGPU();


protected:

    int nx, ny;
    Vector2 terrainSize;
    ScalarField2 bedrock;   // Bedrock
    ScalarField2 ice;       // Ice cover
    ScalarField2 ela;       // ELA map
    ScalarField2 beta;      // Precipitation modifier

    // Simulation params and constants
    double simAccumRate, simAblateRate;
    double simFactorUdeform, simFactorUslide;
    static constexpr double g = 9.81;           // gravity
    static constexpr double rho = 910.0;        // ice density
    static constexpr double Gamma_d = 7.26e-5;  // from [Headley et al. 2012]
    static constexpr double Gamma_s = 3.27;     // from [Headley et al. 2012]

    // Simulation compute shader
    static const unsigned int WORK_GROUP_SIZE_X = 32;
    static const unsigned int WORK_GROUP_SIZE_Y = 32;
    GLuint shaderSIA;
    GLuint shaderSIAaxis, shaderSIAdiag;
    GLuint glbufferBedrock;
    GLuint glbufferIceIn, glbufferIceOut;
    GLuint glbufferELA, glbufferBeta;
    GLuint glbufferD;

    // buffers
    bool initializedBuffers = false;
    bool doPassthrough = true;
    ScalarField2  diffusivity;
    unsigned int  bufferElems;

    // simulation functions
    void stepSIA(double dt);
    double getAdaptiveTimestep();

    // internal funcs
    void updateELA(double avgELA);

};


inline double GlacierTerrain::iceVolume() const
{
    const ScalarField2& field = GetIce();
    double v = 0;
    double a = cellArea();
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            v += field.at(i, j);
        }
    }
    return a * v;
}

inline double GlacierTerrain::elaValue(int i, int j) const
{
    return ela.at(i, j);
}

inline bool GlacierTerrain::aboveELA(int i, int j) const
{
    return (bedrock.at(i, j) + ice.at(i,j)) > ela.at(i,j);
}

inline bool GlacierTerrain::accumArea(int i, int j) const
{
    return aboveELA(i, j);
}

inline double GlacierTerrain::iceThickness(int i, int j) const
{
    return ice.at(i, j);
}

inline Vector2 GlacierTerrain::iceGradient(int i, int j) const
{
    double dip = i < nx - 1 ? Height(i + 1, j) : Height(i, j);
    double dim = i > 0 ? Height(i - 1, j) : Height(i, j);
    double dx = i > 0 && i < nx - 1 ? 2 * cellWidth() : cellWidth();
    double djp = j < ny - 1 ? Height(i, j + 1) : Height(i, j);
    double djm = j > 0 ? Height(i, j - 1) : Height(i, j);
    double dy = j > 0 && j < ny - 1 ? 2 * cellHeight() : cellHeight();
    return Vector2((dip - dim)/dx, (djp - djm)/dy);
}

inline double GlacierTerrain::iceStress(int i, int j) const
{
    return rho * g * iceThickness(i, j) * Norm(iceGradient(i, j));
}

inline double GlacierTerrain::iceSpeed(int i, int j) const
{
    return iceSpeedDeform(i, j) + iceSpeedSlide(i, j);
}

inline double GlacierTerrain::iceSpeedDeform(int i, int j) const
{
    double h = iceThickness(i, j);
    double s = Norm(iceGradient(i, j));
    return simFactorUdeform * Gamma_d * h*h*h*h * s*s;
}

inline double GlacierTerrain::iceSpeedSlide(int i, int j) const
{
    double h = iceThickness(i, j);
    double s = Norm(iceGradient(i, j));
    return simFactorUslide * Gamma_s * h*h * s*s;
}

inline Vector2 GlacierTerrain::iceVelocity(int i, int j) const
{
    return iceSpeed(i, j) * iceGradient(i, j)/Norm(iceGradient(i,j));
}

inline double GlacierTerrain::diffusionTerm(int i, int j) const
{
    return diffusivity.at(i, j);
}

inline double GlacierTerrain::diffusionTermRaw(int i, int j) const
{
    return iceThickness(i,j) * iceSpeed(i, j);
}

inline double GlacierTerrain::yearlyIce(bool withIce) const {
    double iceAccum = 0;
    for (int i = 0; i < nx*ny; i++) {
        iceAccum += std::max(0.0, simAccumRate * beta.at(i) * (bedrock.at(i) + withIce*ice.at(i) - ela.at(i)));
    }
    return iceAccum * cellArea();
}

inline double GlacierTerrain::glaciatedArea() const {
    int iceCells = 0;
    for (int i = 0; i < nx*ny; i++) {
        if (ice.at(i) > 0) iceCells++;
    }
    return iceCells * cellArea();
}

#endif
