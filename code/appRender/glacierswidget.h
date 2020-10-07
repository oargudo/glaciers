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
    GlaciersWidget(QWidget* = nullptr);
    ~GlaciersWidget();

    void SetCamera(const Camera& cam) { camera = cam; };
    Camera GetCamera() { return camera; }

    void setHeightfield(const ScalarField2& hf, bool centerCamera = true);
    void setTexture(const QImage& img);

public slots:
    virtual void mousePressEvent(QMouseEvent*);
    virtual void mouseMoveEvent(QMouseEvent*);
    virtual void mouseReleaseEvent(QMouseEvent*);
    virtual void wheelEvent(QWheelEvent* event);

signals:
    void _signalEditSceneLeft(const Vector3&);
    void _signalEditSceneRight(const Vector3&);

protected:
    virtual void initializeGL();
    virtual void resizeGL(int, int);
    virtual void paintGL();

protected:

    // OpenGL render
    Box3 terrainBBox;
    GLuint shader = 0;
    GLuint meshVAO = 0;
    GLuint bufferVerts = 0, bufferIndices = 0;
    GLuint numTriangles = 0;
    GLuint texId = 0;
    GLuint skyboxShader = 0;
    GLuint skyboxVAO = 0;

    // Camera
    Camera camera;
    bool MoveAt = false;
    int x0 = 0, y0 = 0;
    Vector3 currentAt = Vector3(0);
    Vector3 toAt = Vector3(0);
    int stepAt = 0;

};

#endif // GLACIERSWIDGET_H
