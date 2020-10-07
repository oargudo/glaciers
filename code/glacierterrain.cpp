#include "glacierterrain.h"
#include "shader-utils.h"
#include <chrono>
#include <iostream>


GlacierTerrain::GlacierTerrain()
{
    initializedBuffers = false;
	nx = ny = 0;
}

GlacierTerrain::GlacierTerrain(int nx, int ny) : bedrock(Box2(Vector2(0), Vector2(1)), nx, ny),
                                                 ice(Box2(Vector2(0), Vector2(1)), nx, ny),
                                                 ela(Box2(Vector2(0), Vector2(1)), nx, ny),
                                                 beta(Box2(Vector2(0), Vector2(1)), nx, ny),
                                                 diffusivity(Box2(Vector2(0), Vector2(1)), nx, ny),
												 nx(nx), ny(ny)
{
	initializedBuffers = false;
    terrainSize = Vector2(1, 1);
}

GlacierTerrain::GlacierTerrain(const ScalarField2& b) : bedrock(b), nx(b.getSizeX()), ny(b.getSizeY()),
                                                        ice(b.getDomain(), b.getSizeX(), b.getSizeY()),
                                                        ela(b.getDomain(), b.getSizeX(), b.getSizeY()),
                                                        beta(b.getDomain(), b.getSizeX(), b.getSizeY(), 1.0),
                                                        diffusivity(b.getDomain(), b.getSizeX(), b.getSizeY())

{
    initializedBuffers = false;
    terrainSize = Vector2(b.getDomain().width(), b.getDomain().height());
}

GlacierTerrain::GlacierTerrain(const ScalarField2& b, const ScalarField2& i) : bedrock(b), ice(i), nx(b.getSizeX()), ny(b.getSizeY()),
                                                                               ela(b.getDomain(), b.getSizeX(), b.getSizeY()),
                                                                               beta(b.getDomain(), b.getSizeX(), b.getSizeY(), 1.0),
                                                                               diffusivity(b.getDomain(), b.getSizeX(), b.getSizeY())
{
    initializedBuffers = false;
    terrainSize = Vector2(b.getDomain().width(), b.getDomain().height());
}



GlacierTerrain::~GlacierTerrain() {
    if (initializedBuffers) {
        glDeleteBuffers(1, &glbufferBedrock);
        glDeleteBuffers(1, &glbufferIceIn);
        glDeleteBuffers(1, &glbufferIceOut);
        glDeleteBuffers(1, &glbufferD);
        glDeleteBuffers(1, &glbufferBeta);
        glDeleteBuffers(1, &glbufferELA);
    }
}


void GlacierTerrain::initSimulationGL()
{
    // load shader
    std::string definitions = "";
    definitions += "#define WORK_GROUP_SIZE_X " + std::to_string(WORK_GROUP_SIZE_X) + "\n";
    definitions += "#define WORK_GROUP_SIZE_Y " + std::to_string(WORK_GROUP_SIZE_Y) + "\n";
    shaderSIAaxis = read_program("./Shaders/sia.glsl", (definitions + "#define SIA_DIR_AXES\n").c_str());
    shaderSIAdiag = read_program("./Shaders/sia.glsl", (definitions + "#define SIA_DIR_DIAGONALS\n").c_str());
    std::cout << "Compute shader loaded!" << std::endl;

    if (initializedBuffers) return;

    // create buffers
    bufferElems = ice.getSizeX() * ice.getSizeY();

    glGenBuffers(1, &glbufferBedrock);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferBedrock);
    glBufferData(GL_SHADER_STORAGE_BUFFER, bufferElems * sizeof(double), (const void*)(&bedrock[0]), GL_STATIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    glGenBuffers(1, &glbufferIceIn);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferIceIn);
    glBufferData(GL_SHADER_STORAGE_BUFFER, bufferElems * sizeof(double), (const void*)(&ice[0]), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    glGenBuffers(1, &glbufferELA);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferELA);
    glBufferData(GL_SHADER_STORAGE_BUFFER, bufferElems * sizeof(double), (const void*)(&ela[0]), GL_STATIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    glGenBuffers(1, &glbufferBeta);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferBeta);
    glBufferData(GL_SHADER_STORAGE_BUFFER, bufferElems * sizeof(double), (const void*)(&beta[0]), GL_STATIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    glGenBuffers(1, &glbufferIceOut);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferIceOut);
    glBufferData(GL_SHADER_STORAGE_BUFFER, bufferElems * sizeof(double), NULL, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    glGenBuffers(1, &glbufferD);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferD);
    glBufferData(GL_SHADER_STORAGE_BUFFER, bufferElems * sizeof(double), NULL, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    // run one simulation step with dt=0, to initialize diffusivity buffer
    initializedBuffers = true;
    shaderSIA = shaderSIAaxis;
    doPassthrough = true;
}


void GlacierTerrain::resetSimulation()
{
    ice.fill(0);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferIceIn);
    glBufferData(GL_SHADER_STORAGE_BUFFER, bufferElems * sizeof(double), (const void*)(&ice[0]), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    doPassthrough = true;
}


void GlacierTerrain::setIceMap(const ScalarField2& icemap)
{
    for (int i = 0; i < icemap.getSizeX(); i++) {
        for (int j = 0; j < icemap.getSizeY(); j++) {
            ice(i, j) = icemap.at(i, j);
        }
    }
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferIceIn);
    glBufferData(GL_SHADER_STORAGE_BUFFER, bufferElems * sizeof(double), (const void*)(&ice[0]), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    doPassthrough = true;
}


void GlacierTerrain::updateBedrockGPU()
{
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferBedrock);
    glBufferData(GL_SHADER_STORAGE_BUFFER, bufferElems * sizeof(double), (const void*)(&bedrock[0]), GL_STATIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    doPassthrough = true;
}

void GlacierTerrain::updateIceGPU()
{
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferIceIn);
    glBufferData(GL_SHADER_STORAGE_BUFFER, bufferElems * sizeof(double), (const void*)(&ice[0]), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    doPassthrough = true;
}

void GlacierTerrain::updateELAGPU()
{
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferELA);
    glBufferData(GL_SHADER_STORAGE_BUFFER, bufferElems * sizeof(double), (const void*)(&ela[0]), GL_STATIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}

void GlacierTerrain::updateAccumRateGPU()
{
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferBeta);
    glBufferData(GL_SHADER_STORAGE_BUFFER, bufferElems * sizeof(double), (const void*)(&beta[0]), GL_STATIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}


void GlacierTerrain::configSimulation(const ScalarField2& mapELA, const ScalarField2& mapAccum,
    double accumRate, double ablateRate, double factorUdeform, double factorUslide)
{
    ela = mapELA;
    beta = mapAccum;
    simAccumRate = accumRate;
    simAblateRate = ablateRate;
    simFactorUdeform = factorUdeform;
    simFactorUslide = factorUslide;

    std::cerr << "Config simulation " << std::endl;
    std::cerr << "  - accum rate: " << accumRate << std::endl;
    std::cerr << "  - ablate rate: " << ablateRate << std::endl;
    std::cerr << "  - Udeform * " << factorUdeform << std::endl;
    std::cerr << "  - Uslide * " << factorUslide << std::endl;

    glGenBuffers(1, &glbufferELA);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferELA);
    glBufferData(GL_SHADER_STORAGE_BUFFER, bufferElems * sizeof(double), (const void*)(&ela[0]), GL_STATIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    glGenBuffers(1, &glbufferBeta);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferBeta);
    glBufferData(GL_SHADER_STORAGE_BUFFER, bufferElems * sizeof(double), (const void*)(&beta[0]), GL_STATIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}


double GlacierTerrain::runSimulationSteps(unsigned int numSimulSteps, double maxSimulTime)
{
    const double minTimeStep = 1e-8;
    const double maxTimeStep = 1.0;

    // run n steps and count the ellapsed time
    double deltaTime = 0;
    for (unsigned int i = 0; i < numSimulSteps && deltaTime < maxSimulTime; i++) {
        shaderSIA = i % 2 == 0 ? shaderSIAaxis : shaderSIAdiag;

        double dt = (doPassthrough ? minTimeStep : 1) * std::max(minTimeStep, std::min(maxTimeStep, getAdaptiveTimestep()));
        dt = std::max(0.0, std::min(maxSimulTime - deltaTime, dt));
        stepSIA(dt);

        deltaTime += dt;
        std::swap(glbufferIceIn, glbufferIceOut); // double buffer for ice in/out
        doPassthrough = false;
    }

    // copy the ice buffer from GPU to cpu
    glUseProgram(shaderSIA);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferIceIn); // note we have the most recent buffer here due to swap!
    void* ptr = glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0, bufferElems * sizeof(double), GL_MAP_READ_BIT);
    memcpy(&ice[0], ptr, bufferElems * sizeof(double));
    glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    glUseProgram(0);

    return deltaTime;
}


void GlacierTerrain::stepSIA(double dt)
{
    glUseProgram(shaderSIA);

    glProgramUniform1i(shaderSIA, glGetUniformLocation(shaderSIA, "GridSizeX"), ny); // attention X,Y
    glProgramUniform1i(shaderSIA, glGetUniformLocation(shaderSIA, "GridSizeY"), nx);
    glProgramUniform1d(shaderSIA, glGetUniformLocation(shaderSIA, "dx"), cellHeight());
    glProgramUniform1d(shaderSIA, glGetUniformLocation(shaderSIA, "dy"), cellWidth());
    glProgramUniform1d(shaderSIA, glGetUniformLocation(shaderSIA, "dt"), dt);
    glProgramUniform1d(shaderSIA, glGetUniformLocation(shaderSIA, "betaAccum"), simAccumRate);
    glProgramUniform1d(shaderSIA, glGetUniformLocation(shaderSIA, "betaAblate"), simAblateRate);
    glProgramUniform1d(shaderSIA, glGetUniformLocation(shaderSIA, "factorUdeform"), simFactorUdeform);
    glProgramUniform1d(shaderSIA, glGetUniformLocation(shaderSIA, "factorUslide"), simFactorUslide);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, glbufferBedrock);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, glbufferIceIn);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, glbufferELA);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, glbufferBeta);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, glbufferIceOut);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, glbufferD);

    glDispatchCompute((ny / WORK_GROUP_SIZE_X) + 1, (nx / WORK_GROUP_SIZE_Y) + 1, 1); // attention X,Y
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, 0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, 0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, 0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, 0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, 0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, 0);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, glbufferD);
    void* ptr = glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0, bufferElems * sizeof(double), GL_MAP_READ_BIT);
    memcpy(&diffusivity[0], ptr, bufferElems * sizeof(double));
    glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    glUseProgram(0);
}

double GlacierTerrain::getAdaptiveTimestep()
{
    const double c_stab = 0.125;
    double dx = cellWidth();
    double dy = cellHeight();

    double dmax = 0;
    for (unsigned int i = 0; i < bufferElems; i++) {
        dmax = std::max(dmax, diffusivity[i]);
    }

    return c_stab * std::min(dx*dx, dy*dy)/std::max(dmax, 1.0);
}

ScalarField2 GlacierTerrain::remapIceSurface(const ScalarField2& hiresBed) const
{
    int rx = hiresBed.getSizeX();
    int ry = hiresBed.getSizeY();

    // upsample ice surface (set to 0 where there is no ice)
    ScalarField2 surface(bedrock.getDomain(), nx, ny, 0);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            if (ice.at(i, j) > 0) {
                surface(i, j) = ice.at(i,j) + bedrock.at(i, j);
            }
        }
    }
    ScalarField2 surfaceUp = surface.setResolution(rx, ry);

    // iceUp = surfaceUp - hiresBed
    ScalarField2 iceUp(bedrock.getDomain(), rx, ry, 0);
    for (int i = 0; i < rx; i++) {
        for (int j = 0; j < ry; j++) {
            iceUp(i, j) = std::max(0.0, surfaceUp.at(i, j) - hiresBed.at(i, j));
        }
    }

    // cover holes
    for (int k = 0; k < 2*(rx/nx)*(ry/ny); k++) {
        ScalarField2 surfaceUp(hiresBed);
        surfaceUp += iceUp;

        for (int i = 1; i < rx-1; i++) {
            for (int j = 1; j < ry-1; j++) {

                bool hole = true;
                int numNeighbors = 0;
                double surfaceHeight = 0;

                for (int di = -1; hole && di <= 1; di++) {
                    for (int dj = -1; hole && dj <= 1; dj++) {
                        if (di == 0 && dj == 0) continue;

                        hole = hole && (surfaceUp.at(i + di, j + dj) > surfaceUp(i, j));
                        if (iceUp.at(i + di, j + dj) > 0) {
                            surfaceHeight += surfaceUp.at(i + di, j + dj);
                            numNeighbors++;
                        }
                    }
                }

                if (hole) {
                    iceUp(i, j) = std::max(0.0, surfaceHeight/numNeighbors - hiresBed.at(i,j));
                }
                else if (iceUp.at(i,j) <= 0 && numNeighbors >= 3) {
                    if (iceUp.at(i - 1, j) > 0 && iceUp.at(i + 1, j) <= 0 && surfaceUp.at(i + 1, j) > surfaceUp(i - 1, j))
                        iceUp(i, j) = std::max(0.0, surfaceUp.at(i - 1, j) - hiresBed.at(i, j));
                    if (iceUp.at(i + 1, j) > 0 && iceUp.at(i - 1, j) <= 0 && surfaceUp.at(i - 1, j) > surfaceUp(i + 1, j))
                        iceUp(i, j) = std::max(0.0, surfaceUp.at(i + 1, j) - hiresBed.at(i, j));

                    if (iceUp.at(i, j - 1) > 0 && iceUp.at(i, j + 1) <= 0 && surfaceUp.at(i, j + 1) > surfaceUp(i, j - 1))
                        iceUp(i, j) = std::max(0.0, surfaceUp.at(i, j - 1) - hiresBed.at(i, j));
                    if (iceUp.at(i, j + 1) > 0 && iceUp.at(i, j - 1) <= 0 && surfaceUp.at(i, j - 1) > surfaceUp(i, j + 1))
                        iceUp(i, j) = std::max(0.0, surfaceUp.at(i, j + 1) - hiresBed.at(i, j));
                }
            }
        }
    };

    return iceUp;
}

bool GlacierTerrain::Intersect(const Ray &ray, double &t, Vector3 &q) const
{
    // Compute bounding box
    Box2 box = getDomain();
    double ta, tb;

    // Check the intersection with the bounding box
    Vector3 r0 = ray(0);
    Vector3 r1 = ray(1);
    if (!box.Intersect(Vector2(r0[0], r0[1]), Vector2(r1[0], r1[1]), ta, tb))
        return false;

    t = ta + 0.0001;
    if (ta < 0.0)
    {
        t = 0.0;
    }

    // Ray marching
    while (t < tb)
    {
        // Point along the ray
        Vector3 p = ray(t);
        double h = Height(Vector2(p[0], p[1]));
        if (h > p[2]) {
            q = Vector3(p[0], p[1], h);
            return true;
        }
        else {
            t += 1.0;
        }
    }
    return false;
}
