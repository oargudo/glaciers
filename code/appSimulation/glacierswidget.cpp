#include "glacierswidget.h"
#include <QDateTime>
#include <QTimer>
#include <QPainter>
#include <QMatrix4x4>
#include <QVector3D>
#include <iostream>
#include "shader-utils.h"


GlaciersWidget::GlaciersWidget(QWidget* parent) : QOpenGLWidget(parent)
{
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);
    start = std::chrono::high_resolution_clock::now();
    nbframes = 0;
    shaderCubes = 0;
}

GlaciersWidget::~GlaciersWidget()
{
}


void GlaciersWidget::setGlacierTerrain(GlacierTerrain* glacier)
{
	makeCurrent();

    glacierTerrain = glacier;
    nx = glacier->numCellsX();
    ny = glacier->numCellsY();

    glacierTerrain->initSimulationGL();

    Box2 domain = glacierTerrain->GetBedrock().getDomain();
    terrainBBox = glacierTerrain->getBoundingBox();

    cubeInstancesBedrock.resize(nx * ny);
    cubeInstancesIce.resize(nx * ny);

    fillCubeInstances(ScalarField2(domain, nx, ny, 0),
                      glacierTerrain->GetBedrock(),
                      cubeInstancesBedrock);
    fillCubeInstances(glacierTerrain->GetBedrock(),
                      glacierTerrain->GetIce(),
                      cubeInstancesIce);

    icePrev = glacierTerrain->GetIce();
    texValues = ScalarField2(domain, nx, ny, 0);
    texImg = QImage(nx, ny, QImage::Format::Format_RGBA8888);

	datatexBedVals.resize(nx * ny);
	datatexIceVals.resize(nx * ny);

	updateBedrockTexture();
    updateTexture();
}


void GlaciersWidget::initializeGL()
{
    glewExperimental = true;
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        std::cout << "Error : " << glewGetErrorString(err) << std::endl;
        exit(-1);
    }

    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    shaderCubes = read_program("./shaders/instanced-cube.glsl");
    skyboxShader = read_program("./shaders/skybox.glsl");

    glGenVertexArrays(1, &skyboxVAO);
    initializeCubeVAO();

    unsigned char foo = 255;
    glGenTextures(1, &texId);
    glBindTexture(GL_TEXTURE_2D, texId);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, (void*)(&foo));
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenTextures(1, &texBed);
    glBindTexture(GL_TEXTURE_2D, texBed);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, (void*)(&foo));
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);
	
    QImage arrowImg;
    arrowImg.load("./data/arrow.png");
    glGenTextures(1, &texArrow);
    glBindTexture(GL_TEXTURE_2D, texArrow);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, arrowImg.width(), arrowImg.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, arrowImg.bits());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glGenerateMipmap(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);

	float ffoo = 1.0f;
	glGenTextures(1, &datatexBedrock);
	glBindTexture(GL_TEXTURE_2D, datatexBedrock);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 1, 1, 0, GL_RED, GL_FLOAT, (void*)(&ffoo));
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, 0);

	glGenTextures(1, &datatexIce);
	glBindTexture(GL_TEXTURE_2D, datatexIce);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 1, 1, 0, GL_RED, GL_FLOAT, (void*)(&ffoo));
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, 0);
}

void GlaciersWidget::initializeCubeVAO()
{
    const int faceVerts[36] = { 0, 4, 5, 0, 5, 1,   // plane ymin
                                0, 1, 3, 0, 3, 2,   // plane xmin
                                2, 3, 7, 2, 7, 6,   // plane ymax
                                4, 6, 7, 4, 7, 5,   // plane xmax
                                0, 2, 6, 0, 6, 4,   // plane zmin
                                1, 5, 7, 1, 7, 3 }; // plane zmax
    const float vertNormal[18] = {
        0, -1, 0,
        -1, 0, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, -1,
        0, 0, 1
    };

    std::vector<float> bufferV(36 * 3, 0);
    std::vector<float> bufferN(36 * 3, 0);
    for (int k = 0; k < 36; k++) {
        bufferV[k*3    ] = faceVerts[k] & 0x04 ? 1.f : 0.f;
        bufferV[k*3 + 1] = faceVerts[k] & 0x02 ? 1.f : 0.f;
        bufferV[k*3 + 2] = faceVerts[k] & 0x01 ? 1.f : 0.f;
        bufferN[k*3    ] = vertNormal[int(k / 6) * 3];
        bufferN[k*3 + 1] = vertNormal[int(k / 6) * 3 + 1];
        bufferN[k*3 + 2] = vertNormal[int(k / 6) * 3 + 2];
    }

    glUseProgram(shaderCubes);
    GLuint attribVertexLoc = glGetAttribLocation(shaderCubes, "a_position");
    GLuint attribNormalLoc = glGetAttribLocation(shaderCubes, "a_normal");
    GLuint instanceTranslation = glGetAttribLocation(shaderCubes, "i_translation");
    GLuint instanceHeight = glGetAttribLocation(shaderCubes, "i_height");

    glGenVertexArrays(1, &vaoCube);
    glBindVertexArray(vaoCube);

    GLuint vboVertex;
    glGenBuffers(1, &vboVertex);
    glBindBuffer(GL_ARRAY_BUFFER, vboVertex);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*bufferV.size(), &bufferV[0], GL_STATIC_DRAW);
    glVertexAttribPointer(attribVertexLoc, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(attribVertexLoc);

    GLuint vboNormal;
    glGenBuffers(1, &vboNormal);
    glBindBuffer(GL_ARRAY_BUFFER, vboNormal);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*bufferN.size(), &bufferN[0], GL_STATIC_DRAW);
    glVertexAttribPointer(attribNormalLoc, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(attribNormalLoc);

    glGenBuffers(1, &vboInstanceData);
    glBindBuffer(GL_ARRAY_BUFFER, vboInstanceData);
    glVertexAttribPointer(instanceTranslation, 3, GL_FLOAT, GL_FALSE, sizeof(CubeInstanceData), 0);
    glVertexAttribPointer(instanceHeight, 1, GL_FLOAT, GL_FALSE, sizeof(CubeInstanceData), (GLvoid*)(3*sizeof(float)));
    glEnableVertexAttribArray(instanceTranslation);
    glEnableVertexAttribArray(instanceHeight);
    glVertexAttribDivisor(instanceTranslation, 1);
    glVertexAttribDivisor(instanceHeight, 1);

    glBindVertexArray(0);
    glUseProgram(0);
}


void GlaciersWidget::updateGeometry()
{
    fillCubeInstances(ScalarField2(glacierTerrain->GetBedrock().getDomain(), nx, ny, 0),
                      glacierTerrain->GetBedrock(),
                      cubeInstancesBedrock);
    fillCubeInstances(glacierTerrain->GetBedrock(),
                      glacierTerrain->GetIce(),
                      cubeInstancesIce);
    updateTexture();
}

void GlaciersWidget::reloadShader()
{
    if (shaderCubes) glDeleteProgram(shaderCubes);
    shaderCubes = read_program("./shaders/instanced-cube.glsl");
}

void GlaciersWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, (GLint)w, (GLint)h);
}

void GlaciersWidget::paintGL()
{
    double camDist = Norm(camera.getEye() - terrainBBox.center());
    double tradius = terrainBBox.radius();
    if (camDist < terrainBBox.radius())
        camera.setPlanes(100, 2*tradius);
    else
        camera.setPlanes(camDist - tradius, 2 * tradius + camDist);

    // Clear
    glClearColor(0.62f, 0.74f, 0.85f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Sky
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glUseProgram(skyboxShader);
    glBindVertexArray(skyboxVAO);
    glUniform3f(0, camera.getEye()[0], camera.getEye()[1], camera.getEye()[2]);
    glUniform3f(1, camera.getAt()[0], camera.getAt()[1], camera.getAt()[2]);
    glUniform3f(2, camera.getUp()[0], camera.getUp()[1], camera.getUp()[2]);
    glUniform2f(3, width(), height());
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    // Terrain
    if (glacierTerrain) {
        drawVoxels();
    }

    // FPS Computation
    nbframes++;
    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count();
    double seconds = double(microseconds) / 1000000.0;
    if (seconds >= 1.0)
    {
        fps = nbframes / seconds;
        nbframes = 0;
        start = std::chrono::high_resolution_clock::now();
    }

    // Draw text box overlay
    // We need to unbind VAO and program for now because it causes problem with commands below
    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    QPainter painter;
    painter.begin(this);

    painter.setRenderHint(QPainter::Antialiasing);
    QPen penLineGrey(QColor(50, 50, 50));
    QPen penLineWhite(QColor(250, 250, 250));

    const int bX = 10;
    const int bY = 10;
    const int sizeX = 210;
    const int sizeY = 130;

    // Background
    painter.setPen(penLineGrey);
    painter.fillRect(QRect(bX, bY, sizeX, sizeY), QColor(0, 0, 255, 25));
    painter.drawRect(bX, bY, sizeX, sizeY);

    // Text
    QFont f1, f2;
    f1.setBold(true);
    f1.setPointSize(18);
    f2.setPointSize(10);
    painter.setFont(f1);
    painter.setPen(penLineWhite);
    painter.setFont(f2);
    painter.drawText(10 + 5, bY + 10 +  10, "Year:\t" + QString::number(simulYear, 'f', 1) + (steadyState ? " (steady)" : ""));
	painter.drawText(10 + 5, bY + 10 +  30, "Time:\t" + QString::number(simulComputeTime, 'f', 1) + " s");
    painter.drawText(10 + 5, bY + 10 +  50, "Perf:\t" + QString::number(simulPerf, 'f', 2) + " years/s");
    painter.drawText(10 + 5, bY + 10 +  70, "Avg dIce:\t" + QString::number(simuldIce, 'f', 4) + " m/y");
    painter.drawText(10 + 5, bY + 10 +  90, "Ice Max:\t" + QString::number(simulCurrMaxIce, 'f', 1) + " m");
    painter.drawText(10 + 5, bY + 10 + 110, "Ice Vol:\t" + QString::number(lastIceVol/1e9, 'f', 1) + " km3");
    painter.end();

    // Reset GL depth test
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Schedule next draw
    update();
}


void GlaciersWidget::drawVoxels()
{
    float colorBedrock[4] = { 0.5f, 0.5f, 0.5f, 1.0f };
    float colorIce[4] = { 0.6f, 0.8f, 0.8f, 1.0f };

    glUseProgram(shaderCubes);
    glBindVertexArray(vaoCube);

    // camera
    QMatrix4x4 matPerspective;
    matPerspective.perspective(Math::RadianToDegree(camera.getAngleOfViewV(width(), height())),
                    (GLdouble)width() / (GLdouble)height(),
                    camera.getNearPlane(), camera.getFarPlane());
    glUniformMatrix4fv(glGetUniformLocation(shaderCubes, "ProjectionMatrix"), 1, GL_FALSE, matPerspective.data());

    QMatrix4x4 matView;
    matView.lookAt(QVector3D(camera.getEye()[0], camera.getEye()[1], camera.getEye()[2]),
                   QVector3D(camera.getAt()[0],  camera.getAt()[1],  camera.getAt()[2]),
                   QVector3D(camera.getUp()[0],  camera.getUp()[1],  camera.getUp()[2]));
    glUniformMatrix4fv(glGetUniformLocation(shaderCubes, "ModelViewMatrix"), 1, GL_FALSE, matView.data());

    // common uniforms
    glUniform2f(glGetUniformLocation(shaderCubes, "u_worldSize"), terrainBBox.width(), terrainBBox.height());
    glUniform2f(glGetUniformLocation(shaderCubes, "u_cellSize"),  glacierTerrain->cellWidth(), glacierTerrain->cellHeight());
    glUniform2i(glGetUniformLocation(shaderCubes, "u_gridCells"), nx, ny);
    glUniform3f(glGetUniformLocation(shaderCubes, "u_lightPos"),  terrainBBox.center()[0],
                                                                  terrainBBox.center()[1],
                                                                  100000);
    glUniform4f(glGetUniformLocation(shaderCubes, "u_materialSpecular"), 0.f, 0.f, 0.f, 0.0f);
    glUniform1f(glGetUniformLocation(shaderCubes, "u_shininess"), 20.0f);
	glUniform1i(glGetUniformLocation(shaderCubes, "u_shadeCube"), 1);

	glUniform1i(glGetUniformLocation(shaderCubes, "u_drawCursor"), validAnchor);
	glUniform3f(glGetUniformLocation(shaderCubes, "u_cursorPoint"), anchor[0], anchor[1], anchor[2]);
	glUniform1f(glGetUniformLocation(shaderCubes, "u_cursorRadius"), brushRadius);

    // draw bedrock
    if (showBedrock) {
        glUniform4fv(glGetUniformLocation(shaderCubes, "u_materialAmbient"), 1, colorBedrock);
		glUniform4f(glGetUniformLocation(shaderCubes, "u_materialDiffuse"), 0.08f, 0.08f, 0.08f, 1.0f);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texBed);
        glUniform1i(glGetUniformLocation(shaderCubes, "u_texture"), 0);
        glUniform1i(glGetUniformLocation(shaderCubes, "u_useTexture"), 2);

        glBindBuffer(GL_ARRAY_BUFFER, vboInstanceData);
        glBufferData(GL_ARRAY_BUFFER, sizeof(CubeInstanceData) * cubeInstancesBedrock.size(), &cubeInstancesBedrock[0], GL_DYNAMIC_DRAW);
        glDrawArraysInstanced(GL_TRIANGLES, 0, 36, GLsizei(cubeInstancesBedrock.size()));
    }

    // draw ice
    if (showIce) {
        glUniform4fv(glGetUniformLocation(shaderCubes, "u_materialAmbient"), 1, colorIce);
		glUniform4f(glGetUniformLocation(shaderCubes, "u_materialDiffuse"), 0.15f, 0.15f, 0.15f, 1.0f);
		
        // compute directly the texture using the shader
        if (activeTexType == GlacierTextureType::FEAT_SHADER) {
            glUniform1i(glGetUniformLocation(shaderCubes, "u_useTexture"), 3);

			glUniform3f(glGetUniformLocation(shaderCubes, "u_bedRange"), float(bedMin), float(bedMax), float(bedMax - bedMin));
			glUniform3f(glGetUniformLocation(shaderCubes, "u_iceRange"), float(iceMin), float(iceMax), float(iceMax - iceMin));

            glActiveTexture(GL_TEXTURE1);
            glBindTexture(GL_TEXTURE_2D, texArrow);
            glUniform1i(glGetUniformLocation(shaderCubes, "u_arrowTexture"), 1);

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, datatexBedrock);
			glUniform1i(glGetUniformLocation(shaderCubes, "u_datatexBed"), 2);

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, datatexIce);
			glUniform1i(glGetUniformLocation(shaderCubes, "u_datatexIce"), 3);
        }
		else {
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, texId);
			glUniform1i(glGetUniformLocation(shaderCubes, "u_texture"), 0);
			glUniform1i(glGetUniformLocation(shaderCubes, "u_useTexture"), 1);
		}

        glBindBuffer(GL_ARRAY_BUFFER, vboInstanceData);
        glBufferData(GL_ARRAY_BUFFER, sizeof(CubeInstanceData) * cubeInstancesIce.size(), &cubeInstancesIce[0], GL_DYNAMIC_DRAW);
        glDrawArraysInstanced(GL_TRIANGLES, 0, 36, GLsizei(cubeInstancesIce.size()));
    }

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}


void GlaciersWidget::fillCubeInstances(const ScalarField2& fieldBase, const ScalarField2& fieldHeight, std::vector<CubeInstanceData>& cubeInstances)
{
    double dx = glacierTerrain->cellWidth();
    double dy = glacierTerrain->cellHeight();
    Vector3 pminBox = glacierTerrain->getBoundingBox().getMin();

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            cubeInstances[i * ny + j].x = pminBox[0] + i*dx;
            cubeInstances[i * ny + j].y = pminBox[1] + j*dy;
            cubeInstances[i * ny + j].z = fieldBase.at(i, j);
            cubeInstances[i * ny + j].height = fieldHeight.at(i, j);
        }
    }
}

void GlaciersWidget::runSimulation(bool pauseSteady)
{
    simulationRunning = true;
    steadyState = false;
    pauseWhenSteady = pauseSteady;
    yearlyIce = glacierTerrain->yearlyIce(false);

    if (!simulQtimer) {
        simulQtimer = new QTimer(this);
        connect(simulQtimer, SIGNAL(timeout()), this, SLOT(simulationTask()));
    }
    simulQtimer->start(0);
}

void GlaciersWidget::pauseSimulation()
{
    simulationRunning = false;
}

void GlaciersWidget::resetSimulation()
{
    simulYear = 0;
    steadyState = false;
    simulComputeTime = 0;
    simulSteps = 0;
    yearlyIce = glacierTerrain->yearlyIce(false);

    fillCubeInstances(glacierTerrain->GetBedrock(), glacierTerrain->GetIce(), cubeInstancesIce);
}


void GlaciersWidget::simulationTask()
{
    if (simulationRunning && (!steadyState || !pauseWhenSteady)) {

        yearlyIce = glacierTerrain->yearlyIce(true);
        icePrev = glacierTerrain->GetIce();

        // simulation
        auto start = std::chrono::high_resolution_clock::now();

        double dt = glacierTerrain->runSimulationSteps(SIMUL_STEPS_BATCH, 1e3);
        simulYear += dt;
        simulSteps += SIMUL_STEPS_BATCH;
        minSimulYears -= dt;

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = finish - start;
        simulComputeTime += duration.count();

        double currIceVol = glacierTerrain->iceVolume();
        double dIce = currIceVol - lastIceVol;
        double dIce_dt = dIce / dt;
        double dIce_dt_avg = dIce_dt / glacierTerrain->terrainArea();
		        
        simulPerf = dt / duration.count();
        double minIce;
        glacierTerrain->GetIce().getRange(minIce, simulCurrMaxIce);

		/*steadyState = minSimulYears < 0 && (std::abs(dIceRate) < 0.01 * yearlyIce * steadyRatio
										 || std::abs(dIceDensityRate) < 0.001*steadyRatio);*/
		steadyState = minSimulYears < 0 && std::abs(dIce_dt_avg) < dIceSteadyCondition;
		lastIceVol = currIceVol;
        simuldIce = dIce_dt_avg;
		
        // graphics
        fillCubeInstances(glacierTerrain->GetBedrock(), glacierTerrain->GetIce(), cubeInstancesIce);
        updateTexture();

        if (!steadyState) {
            simulQtimer->start(0);
        }
        update();
    }
}

QImage GlaciersWidget::shadedBedrock()
{
	QImage shading(nx, ny, QImage::Format_RGBA8888);

	const ScalarField2& hf = glacierTerrain->GetBedrock();
	double a, b;
	hf.getRange(a, b);

	Vector3 lightDir = Normalized(Vector3(-1.0, -1.0, 2.5));
	ColorPalette reliefPalette = ColorPalette::Relief();
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			// Height shading: color is a linear interpolation of height colors
			double t = 300.0 * Math::LinearStep(hf.at(i, j), -b, 1.5*b);
			Vector3 cz = reliefPalette.getColor(t);

			double light = dot(hf.normal(i, j), lightDir);
			light = 0.5 * (1.0 + light);
			double cs = 0.9 * light;

			// Cosine like
			cs *= cs;

			// Normal shading: color is a combination of cool and cold colors according to the orientation
			Vector3 ambient = 0.25 * Vector3(1.0, 1.0, 1.0);
			Vector3 c1 = 0.25 * cs * cz;
			Vector3 c2 = 0.5 * cs * Vector3(1.0, 1.0, 1.0);
			Vector3 c = ambient + c1 + c2;
			c = c - Vector3(0.05);

			shading.setPixelColor(i, j, c.toQColor().rgb());
		}
	}

	return shading;
}

void GlaciersWidget::updateBedrockTexture()
{
	QImage shading = shadedBedrock();
	glBindTexture(GL_TEXTURE_2D, texBed);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, nx, ny, 0, GL_RGBA, GL_UNSIGNED_BYTE, shading.bits());
	glGenerateMipmap(GL_TEXTURE_2D);
}

void GlaciersWidget::updateTexture()
{
    double maxAbsVal = 0;

    switch (activeTexType) {

        case GlacierTextureType::NONE:
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
					texImg.setPixelColor(i, j, QColor(178, 188, 188, 255));
                }
            }
            break;

        case GlacierTextureType::ELA:
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    if (glacierTerrain->aboveELA(i, j))
                        texImg.setPixelColor(i, j, QColor(200, 200, 200, 255));
                    else
                        texImg.setPixelColor(i, j, QColor(178, 188, 188, 255));
                }
            }
            break;

        case GlacierTextureType::THICKNESS:
            texImg = glacierTerrain->GetIce().CreateImage(ColorPalette::CoolWarm()).convertToFormat(QImage::Format_RGBA8888);
            break;

        case GlacierTextureType::GRADIENT:
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    texValues(i, j) = Norm(glacierTerrain->iceGradient(i, j));
                }
            }
            texImg = texValues.CreateImage(ColorPalette::CoolWarm()).convertToFormat(QImage::Format_RGBA8888);
            break;

        case GlacierTextureType::STRESS:
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    texValues(i, j) = glacierTerrain->iceStress(i, j);
                }
            }
            texImg = texValues.CreateImage(ColorPalette::CoolWarm()).convertToFormat(QImage::Format_RGBA8888);
            break;

        case GlacierTextureType::DIFFUSIVITY:
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    texValues(i, j) = glacierTerrain->diffusionTerm(i, j);
                }
            }
            texImg = texValues.CreateImage(ColorPalette::CoolWarm()).convertToFormat(QImage::Format_RGBA8888);
            break;

        case GlacierTextureType::DIFFUSIVITY_RAW:
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    texValues(i, j) = glacierTerrain->diffusionTermRaw(i, j);
                }
            }
            texImg = texValues.CreateImage(ColorPalette::CoolWarm()).convertToFormat(QImage::Format_RGBA8888);
            break;

        case GlacierTextureType::SPEED_DEFORM:
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    texValues(i, j) = glacierTerrain->iceSpeedDeform(i, j);
                }
            }
            texImg = texValues.CreateImage(ColorPalette::CoolWarm()).convertToFormat(QImage::Format_RGBA8888);
            break;

        case GlacierTextureType::SPEED_SLIP:
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    texValues(i, j) = glacierTerrain->iceSpeedSlide(i, j);
                }
            }
            texImg = texValues.CreateImage(ColorPalette::CoolWarm()).convertToFormat(QImage::Format_RGBA8888);
            break;

        case GlacierTextureType::SPEED_DOMINANT_TYPE:
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    if (glacierTerrain->iceThickness(i, j) <= 0) {
                        texImg.setPixelColor(i, j, QColor(50, 50, 50, 255));
                        continue;
                    }
                    double ud = glacierTerrain->iceSpeedDeform(i, j);
                    double us = glacierTerrain->iceSpeedSlide(i, j);
                    if (ud > us)      texImg.setPixelColor(i, j, QColor(50, 50, 200, 255));
                    else if (ud < us) texImg.setPixelColor(i, j, QColor(200, 50, 50, 255));
                    else	          texImg.setPixelColor(i, j, QColor(200, 200, 200, 255));
                }
            }
            break;

        case GlacierTextureType::VELOCITY:
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    Vector2 v = glacierTerrain->iceGradient(i, j) / Norm(glacierTerrain->iceGradient(i, j));
                    texImg.setPixelColor(i, j, QColor((unsigned char)(v[0] * 128 + 128), (unsigned char)(v[1] * 128 + 128), 128, 255));
                }
            }
            break;

        case GlacierTextureType::ICE_DIFFERENCE:
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    texValues(i, j) = -(glacierTerrain->Ice(i, j) - icePrev.at(i, j));
                    maxAbsVal = std::max(std::abs(texValues.at(i,j)), maxAbsVal);
                }
            }
            maxAbsVal = std::max(maxAbsVal, 1e-6);
            texImg = texValues.CreateImage(-maxAbsVal, maxAbsVal, ColorPalette::CoolWarm()).convertToFormat(QImage::Format_RGBA8888);
            break;

        case GlacierTextureType::FEAT_SHADER:
        {
			glacierTerrain->GetBedrock().getRange(bedMin, bedMax);
			glacierTerrain->GetIce().getRange(iceMin, iceMax);

			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					datatexBedVals[glacierTerrain->vertexIndex(i, j)] = float((glacierTerrain->Bedrock(i, j) - bedMin) / (bedMax - bedMin));
					datatexIceVals[glacierTerrain->vertexIndex(i, j)] = float((glacierTerrain->iceThickness(i, j) - iceMin) / (iceMax - iceMin));
				}
			}

			glBindTexture(GL_TEXTURE_2D, datatexBedrock);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, nx, ny, 0, GL_RED, GL_FLOAT, (void*)(&datatexBedVals[0]));
			glGenerateMipmap(GL_TEXTURE_2D);

			glBindTexture(GL_TEXTURE_2D, datatexIce);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, nx, ny, 0, GL_RED, GL_FLOAT, (void*)(&datatexIceVals[0]));
			glGenerateMipmap(GL_TEXTURE_2D);
			
            return;
        }
        break;
        default:
            break;
    }

    setTexture(texImg);
}

void GlaciersWidget::setTexture(const QImage& img)
{
    glBindTexture(GL_TEXTURE_2D, texId);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, nx, ny, 0, GL_RGBA, GL_UNSIGNED_BYTE, img.bits());
    glGenerateMipmap(GL_TEXTURE_2D);
}

void GlaciersWidget::setRenderType(bool voxels, bool bedrock, bool ice)
{
    renderAsVoxels = voxels;
    showBedrock = bedrock;
    showIce = ice;
}


void GlaciersWidget::mousePressEvent(QMouseEvent * e)
{
    x0 = e->globalX();
    y0 = e->globalY();

    if (e->modifiers() & Qt::ControlModifier) {
		Ray ray = camera.pixelToRay(e->pos().x(), e->pos().y() - 1, width(), height());
		double t;
		validAnchor = glacierTerrain->Intersect(ray, t, anchor);

		if (validAnchor) {
			if (e->buttons() & Qt::LeftButton) {
				emit _signalEditSceneLeft(anchor);
			}
			if (e->buttons() & Qt::RightButton) {
				emit _signalEditSceneRight(anchor);
			}
		}
    }
	else {
		validAnchor = false;
	}

    update();
}

void GlaciersWidget::mouseMoveEvent(QMouseEvent * e)
{
    int x = e->globalX();
    int y = e->globalY();

    if (e->modifiers() & Qt::ControlModifier) {
		const QPoint& pixel = e->pos();
		Ray ray = camera.pixelToRay(pixel.x(), pixel.y() - 1, width(), height());
		double t;
		validAnchor = glacierTerrain->Intersect(ray, t, anchor);

		if (validAnchor) {
			if (e->buttons() & Qt::LeftButton) {
				emit _signalEditSceneLeft(anchor);
			}
			if (e->buttons() & Qt::RightButton) {
				emit _signalEditSceneRight(anchor);
			}
		}
    }
	else {
		if (e->buttons() & Qt::LeftButton) {
			camera.leftRightRound((x0 - x) * 0.01);
			camera.upDownRound((y0 - y) * 0.005);
		}
		else if (e->buttons() & Qt::RightButton) {
			Vector3 previousAt = camera.getAt();
			Vector3 direction = camera.getAt() - camera.getEye();
			double currentDist = 0.002f * Norm(direction);
			camera.backForth((y - y0) * currentDist);
			camera.setAt(previousAt);
		}
		else if (e->buttons() & Qt::MidButton) {
			camera.leftRightPlane((x - x0) * 1.0);
			camera.upDownPlane((y - y0) * 1.0);
		}

		x0 = e->globalX();
		y0 = e->globalY();

		validAnchor = false;
	}

    update();
}

void GlaciersWidget::mouseReleaseEvent(QMouseEvent*)
{
	validAnchor = false;
	updateBedrockTexture();
	update();
}

void GlaciersWidget::wheelEvent(QWheelEvent* e)
{
    int numDegrees = e->angleDelta().y() / 8;
    int numSteps = numDegrees / 15;
	
    if (!(e->modifiers()&Qt::ShiftModifier)) {
        Vector3 direction = camera.getAt() - camera.getEye();
        double currentDist = Norm(direction);
        direction = Normalized(direction);
        double speed = 0.1*currentDist;

        camera.setEye(camera.getEye() + direction * speed * numSteps);
    }

    update();
}
