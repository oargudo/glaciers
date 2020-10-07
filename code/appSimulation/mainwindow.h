#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "glacierswidget.h"
#include "glacierterrain.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

public slots:

    void presetScene1();
    void presetScene2();
    void presetScene3();
    void presetScene4();
    void presetScene5();
    void presetScene6();
    void presetScene7();
    void presetScene8();
    void presetScene9();

	void selectDEM();
    void loadDEM();
    void loadELA();
    void loadPrecipitation();
    void loadInitialIce();

    void Render(bool resetCamera = true);
    void RenderUpdate();

    void editingSceneLeft(const Vector3& p);
    void editingSceneRight(const Vector3& p);
	void updateBrushRadius(double d);

    void runGlacierSimulation();
    void pauseGlacierSimulation();
    void resetGlacierSimulation();
    void upscaleGlacier();
	void updateSteadyCondition(double);

    void saveScene();

private:
    void createActions();

    void configSimulation();

	void loadBedrocks(const QString& filename, const std::vector<int>& downscales, double terrainKM, const std::vector<double>& elevRange = {});
    void loadMapsELA(const QString& filename);
    void loadMapsPrecipitation(const QString& filename);
    void loadMapsIce(const QString& filename);

    ScalarField2 getVariableELA(double avgELA, double devELA);
	void updatePreviewELA(const ScalarField2& ela);
	void updatePreviewAccum(const ScalarField2& accum);

private:
    Ui::MainWindow *ui;

    GlaciersWidget *glaciersWidget;
    GlacierTerrain  glacierTerrain;

    // multiresolution bedrocks and control maps
	QString demFilename = "";
    std::vector<ScalarField2> bedrocksMultires;
    std::vector<ScalarField2> factorElaDeviation;
    std::vector<ScalarField2> factorBetaMaps;
    std::vector<ScalarField2> initialIceMaps;
    int multiresIndex;

    // brushes
    double brushRadius = 750;
    const double rockStrength = 20;
    const double iceStrength = 10;
    const double elaStrength = 20;
    const double precStrenth = 0.1;

};
#endif // MAINWINDOW_H
