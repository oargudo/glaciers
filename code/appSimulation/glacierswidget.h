#ifndef GLACIERSWIDGET_H
#define GLACIERSWIDGET_H

#include "glew.h"
#include <QtGui/QKeyEvent>
#include <QtGui/QMouseEvent>
#include <QtWidgets/QOpenGLWidget>
#include "scalarfield2.h"
#include "glacierterrain.h"


class GlaciersWidget : public QOpenGLWidget
{
    Q_OBJECT

public:

    enum class GlacierTextureType {
        NONE,
        ELA,
        THICKNESS,
        GRADIENT,
        STRESS,
        DIFFUSIVITY,
        DIFFUSIVITY_RAW,
        SPEED_DEFORM,
        SPEED_SLIP,
        SPEED_DOMINANT_TYPE,
        VELOCITY,
        ICE_DIFFERENCE,
        FEAT_SHADER
    };


    GlaciersWidget(QWidget* = nullptr);
    ~GlaciersWidget();

    void setGlacierTerrain(GlacierTerrain* glacier);

    void SetCamera(const Camera& cam) { camera = cam; };
    Camera GetCamera() { return camera; }

    void setRenderType(bool voxels, bool bedrock, bool ice);
	void setTextureType(GlacierTextureType texType) { activeTexType = texType; updateTexture(); }
	void updateBedrockTexture();
    void updateTexture();
    const QImage& getTexture() { return texImg; };
    void setTexture(const QImage& img);

    void setSteadyCondition(double d) { dIceSteadyCondition = d; }
	double getLastdIce() { return simuldIce; }
    void setMinimumSimulationYears(double y) { minSimulYears = y; }

    void runSimulation(bool pauseStationary); // remember to config parameters before this call, directly on the glacier heightfield
    void pauseSimulation();
    void resetSimulation();
	
	QImage shadedBedrock();


public slots:
    virtual void mousePressEvent(QMouseEvent*);
    virtual void mouseMoveEvent(QMouseEvent*);
	virtual void mouseReleaseEvent(QMouseEvent*);
    virtual void wheelEvent(QWheelEvent* event);

    void simulationTask();

    void updateGeometry();
    void reloadShader();
	void updateBrushRadius(double d) { brushRadius = d; };

signals:
    void _signalEditSceneLeft(const Vector3&);
    void _signalEditSceneRight(const Vector3&);

protected:

    struct CubeInstanceData {
        float x, y, z;
        float height;
    };

    virtual void initializeGL();
    virtual void resizeGL(int, int);
    virtual void paintGL();

    void initializeCubeVAO();

    void fillCubeInstances(const ScalarField2& fieldBottom, const ScalarField2& fieldTop, std::vector<CubeInstanceData>& cubeInstances);

    void drawVoxels();


protected:

    // Glacier simulator
    GlacierTerrain* glacierTerrain = nullptr;
    Box3 terrainBBox;
    int nx = 0, ny = 0;
    ScalarField2 icePrev;
    const unsigned int SIMUL_STEPS_BATCH = 16;
    QTimer* simulQtimer = nullptr;

    // OpenGL render
    GLuint shaderCubes = 0;
    GLuint vaoCube = 0;
    GLuint vboInstanceData = 0;
    GLuint texBed = 0;
    GLuint texId = 0;
    GLuint texArrow = 0;
    GLuint skyboxShader = 0;
    GLuint skyboxVAO = 0;
    std::vector<CubeInstanceData> cubeInstancesBedrock;
    std::vector<CubeInstanceData> cubeInstancesIce;

    bool renderAsVoxels, showBedrock, showIce;
    ScalarField2 texValues;
    QImage texImg;
    GlacierTextureType activeTexType;

	// Custom shader data
	std::vector<float> datatexBedVals;
	std::vector<float> datatexIceVals;
	GLuint datatexBedrock = 0;
	GLuint datatexIce = 0;
	double bedMin, bedMax;
	double iceMin, iceMax;


    // Camera
    Camera camera;
    bool MoveAt = false;
    int x0 = 0, y0 = 0;
    Vector3 currentAt = Vector3(0);
    Vector3 toAt = Vector3(0);
    int stepAt = 0;

    // Timer
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    int nbframes;
    int fps = 0;

    // Simul
    double maxSimulYears = 1e5;
    double simulYear = 0;
    double simulComputeTime = 0;
    double simulPerf = 0;
    double simuldIce = 0;
    double simulCurrMaxIce = 0;
    int simulSteps = 0;
    bool steadyState = false;
    bool simulationRunning = true;
    bool pauseWhenSteady = true;
	double dIceSteadyCondition = 0.001;
    double minSimulYears = 0;
    double lastIceVol = 0;
    double yearlyIce = 0;

    // Editor
	bool validAnchor;
    Vector3 anchor;
	double brushRadius = 750.0;
};

#endif // GLACIERSWIDGET_H
