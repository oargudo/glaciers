#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include <QDir>
#include <QFileDialog>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    glaciersWidget = new GlaciersWidget();
	QGridLayout* GLlayout = new QGridLayout;
	GLlayout->addWidget(glaciersWidget, 0, 0);
	GLlayout->setContentsMargins(0, 0, 0, 0);
	ui->glWidget->setLayout(GLlayout);
    glaciersWidget->SetCamera(Camera(Vector3(-10.0, -10.0, 10.0), Vector3(0.0, 0.0, 0.0)));

    createActions();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::createActions()
{
    connect(glaciersWidget, SIGNAL(_signalEditSceneLeft(const Vector3&)), this, SLOT(editingSceneLeft(const Vector3&)));
    connect(glaciersWidget, SIGNAL(_signalEditSceneRight(const Vector3&)), this, SLOT(editingSceneRight(const Vector3&)));

    connect(ui->btnLoadScene1, SIGNAL(clicked()), this, SLOT(presetScene1()));
    connect(ui->btnLoadScene2, SIGNAL(clicked()), this, SLOT(presetScene2()));
    connect(ui->btnLoadScene3, SIGNAL(clicked()), this, SLOT(presetScene3()));
    connect(ui->btnLoadScene4, SIGNAL(clicked()), this, SLOT(presetScene4()));
    connect(ui->btnLoadScene5, SIGNAL(clicked()), this, SLOT(presetScene5()));
    connect(ui->btnLoadScene6, SIGNAL(clicked()), this, SLOT(presetScene6()));
    connect(ui->btnLoadScene7, SIGNAL(clicked()), this, SLOT(presetScene7()));
    connect(ui->btnLoadScene8, SIGNAL(clicked()), this, SLOT(presetScene8()));
    connect(ui->btnLoadScene9, SIGNAL(clicked()), this, SLOT(presetScene9()));

	connect(ui->btnSelectDEM, SIGNAL(clicked()), this, SLOT(selectDEM()));
	connect(ui->btnLoadDEM, SIGNAL(clicked()), this, SLOT(loadDEM()));
    connect(ui->btnLoadELA, SIGNAL(clicked()), this, SLOT(loadELA()));
    connect(ui->btnLoadPrec, SIGNAL(clicked()), this, SLOT(loadPrecipitation()));
    connect(ui->btnLoadInitialIce, SIGNAL(clicked()), this, SLOT(loadInitialIce()));

    connect(ui->btnGlacierSimRun, SIGNAL(clicked()), this, SLOT(runGlacierSimulation()));
    connect(ui->btnGlacierSimPause, SIGNAL(clicked()), this, SLOT(pauseGlacierSimulation()));
    connect(ui->btnGlacierSimReset, SIGNAL(clicked()), this, SLOT(resetGlacierSimulation()));
    connect(ui->btnGlacierUpsample, SIGNAL(clicked()), this, SLOT(upscaleGlacier()));
	connect(ui->sb_dIce, SIGNAL(valueChanged(double)), this, SLOT(updateSteadyCondition(double)));

	connect(ui->sb_brushRadius, SIGNAL(valueChanged(double)), this, SLOT(updateBrushRadius(double)));

    connect(ui->btnSaveScene, SIGNAL(clicked()), this, SLOT(saveScene()));

    connect(ui->radioButton_Voxels, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));
    connect(ui->radioButton_Mesh, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));

	connect(ui->checkBox_RenderBedrock, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));
	connect(ui->checkBox_RenderIce, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));

    connect(ui->radioButton_ShadingIceNeutral, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));
    connect(ui->radioButton_ShadingIceELA, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));
    connect(ui->radioButton_ShadingIceThick, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));
    connect(ui->radioButton_ShadingIceGradient, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));
    connect(ui->radioButton_ShadingIceStress, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));
    connect(ui->radioButton_ShadingIceDiffusivity, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));
    connect(ui->radioButton_ShadingIceDiffusivityRaw, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));
    connect(ui->radioButton_ShadingIceUDeform, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));
    connect(ui->radioButton_ShadingIceUSlip, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));
    connect(ui->radioButton_ShadingIceUDominant, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));
    connect(ui->radioButton_ShadingIceDifference, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));
    connect(ui->radioButton_ShadingFeatShader, SIGNAL(toggled(bool)), this, SLOT(RenderUpdate()));

	connect(ui->btnReloadShader, SIGNAL(clicked()), glaciersWidget, SLOT(reloadShader()));
}

void MainWindow::presetScene1()
{
    loadBedrocks("data/terrains/aiguestortes_20m.png", { 1, 2, 3, 6 }, 30);
    loadMapsELA("data/terrains/aiguestortes_sun.png");
    //loadMapsPrecipitation("data/terrains/aiguestortes_precipitation.png");
    ui->checkUseELAMap->setChecked(true);
    ui->checkUsePrecMap->setChecked(false);
    ui->sb_ELA->setValue(2500);
    ui->sb_ELAdev->setValue(200);
    ui->sb_AccumRate->setValue(0.002);
	ui->sb_AblateRate->setValue(0.001);
}

void MainWindow::presetScene2()
{
    loadBedrocks(QString("data/terrains/ecrins_30m.png"), { 1, 2, 4, 8 }, 84);
    loadMapsELA("data/terrains/ecrins_sun.png");
    factorBetaMaps.clear();
    ui->checkUseELAMap->setChecked(true);
    ui->checkUsePrecMap->setChecked(false);
    ui->sb_ELA->setValue(2600.0);
    ui->sb_ELAdev->setValue(200.0);
    ui->sb_AccumRate->setValue(0.001);
	ui->sb_AblateRate->setValue(0.001);
}

void MainWindow::presetScene3()
{
    loadBedrocks(QString("data/terrains/canyon_10m.png"), { 1, 2, 4, 6, 12 }, 30);
    loadMapsELA("data/terrains/canyon_sun.png");
    loadMapsPrecipitation("data/terrains/canyon_precipitation.png");
    ui->checkUseELAMap->setChecked(true);
    ui->checkUsePrecMap->setChecked(true);
    ui->sb_ELA->setValue(2000.0);
    ui->sb_ELAdev->setValue(300.0);
    ui->sb_AccumRate->setValue(0.002);
	ui->sb_AblateRate->setValue(0.001);
}

void MainWindow::presetScene4()
{
    loadBedrocks("data/terrains/mtperdu_20m.png", { 1, 2, 3, 6 }, 30);
    loadMapsELA("data/terrains/mtperdu_sun.png");
    factorBetaMaps.clear();
    ui->checkUseELAMap->setChecked(true);
    ui->checkUsePrecMap->setChecked(false);
    ui->sb_ELA->setValue(2500);
    ui->sb_ELAdev->setValue(300);
    ui->sb_AccumRate->setValue(0.002);
	ui->sb_AblateRate->setValue(0.001);
}

void MainWindow::presetScene5()
{
    loadBedrocks("data/terrains/sthelens_10m.png", { 2, 4, 8 }, 20);
    loadMapsELA("data/terrains/sthelens_sun.png");
    factorBetaMaps.clear();
    ui->checkUseELAMap->setChecked(true);
    ui->checkUsePrecMap->setChecked(false);
    ui->sb_ELA->setValue(2000);
    ui->sb_ELAdev->setValue(200);
    ui->sb_AccumRate->setValue(0.005);
	ui->sb_AblateRate->setValue(0.001);
}

void MainWindow::presetScene6()
{
    loadBedrocks(QString("data/terrains/ruapehu_25m.png"), { 1, 2, 3, 6 }, 30);
    loadMapsELA("data/terrains/ruapehu_sun.png");
    factorBetaMaps.clear();
    ui->checkUseELAMap->setChecked(true);
    ui->checkUsePrecMap->setChecked(false);
    ui->sb_ELA->setValue(2000);
    ui->sb_ELAdev->setValue(300);
    ui->sb_AccumRate->setValue(0.003);
	ui->sb_AblateRate->setValue(0.001);
}

void MainWindow::presetScene7()
{
    loadBedrocks(QString("data/terrains/smokies_10m.png"), { 1, 2, 4, 6, 12 }, 30);
    loadMapsELA("data/terrains/smokies_sun.png");
    factorBetaMaps.clear();
    ui->checkUseELAMap->setChecked(true);
    ui->checkUsePrecMap->setChecked(false);
    ui->sb_ELA->setValue(1500);
    ui->sb_ELAdev->setValue(0);
    ui->sb_AccumRate->setValue(0.001);
	ui->sb_AblateRate->setValue(0.001);
}

void MainWindow::presetScene8()
{
    loadBedrocks(QString("data/terrains/kungsleden_25m.png"), { 1, 3, 6, 12 }, 60);
    loadMapsELA("data/terrains/kungsleden_sun.png");
    factorBetaMaps.clear();
    ui->checkUseELAMap->setChecked(true);
    ui->checkUsePrecMap->setChecked(false);
    ui->sb_ELA->setValue(1400);
    ui->sb_ELAdev->setValue(0);
    ui->sb_AccumRate->setValue(0.005);
	ui->sb_AblateRate->setValue(0.001);
}

void MainWindow::presetScene9()
{
    loadBedrocks(QString("data/terrains/chartreuse_edit_20m.png"), { 1, 2, 5, 10 }, 50);
    loadMapsELA("data/terrains/chartreuse_sun.png");
    factorBetaMaps.clear();
    ui->checkUseELAMap->setChecked(true);
    ui->checkUsePrecMap->setChecked(false);
    ui->sb_ELA->setValue(1400);
    ui->sb_ELAdev->setValue(100);
    ui->sb_AccumRate->setValue(0.003);
	ui->sb_AblateRate->setValue(0.001);
}

void MainWindow::selectDEM()
{
	QString filename = QFileDialog::getOpenFileName(nullptr,
		"Choose a filename to load", "", "png file (*.png)");
	if (!filename.isNull()) {
		demFilename = filename;
		ui->btnSelectDEM->setText(QFile(filename).fileName());
	}
}

void MainWindow::loadDEM()
{
	if (demFilename.length() < 1) return;

	double demKM = ui->sb_DEM_km->value();
	double hmin = ui->sb_DEM_hmin->value();
	double hmax = ui->sb_DEM_hmax->value();
	std::vector<int> scales;
	QStringList ds = ui->le_DEM_scales->text().split(",");
	for (QString ss : ds) {
		bool ok = false;
		int s = ss.toInt(&ok);
		if (!ok) s = 1;
		scales.push_back(s);
	}

	loadBedrocks(demFilename, scales, demKM, { hmin, hmax });
}

void MainWindow::loadELA()
{
    QString filename = QFileDialog::getOpenFileName(nullptr,
        "Choose a filename to load", "", "png file (*.png)");
    if (!filename.isNull()) {
        loadMapsELA(filename);
    }
}

void MainWindow::loadPrecipitation()
{
    QString filename = QFileDialog::getOpenFileName(nullptr,
        "Choose a filename to load", "", "png file (*.png)");
    if (!filename.isNull()) {
        loadMapsPrecipitation(filename);
    }
}

void MainWindow::loadInitialIce()
{
    QString filename = QFileDialog::getOpenFileName(nullptr,
        "Choose a filename to load", "", "png file (*.png)");
    if (!filename.isNull()) {
        loadMapsIce(filename);
    }
}

void MainWindow::Render(bool resetCamera)
{
    if (glacierTerrain.isEmpty())
        return;

    if (resetCamera) {
        Box3 tbox = glacierTerrain.getBoundingBox();
        Camera cc = Camera::View(tbox);
        glaciersWidget->SetCamera(cc);
    }

    // Geometry
    bool b = ui->checkBox_RenderBedrock->isChecked();
    bool i = ui->checkBox_RenderIce->isChecked();
    if (ui->radioButton_Voxels->isChecked()) {
        glaciersWidget->setRenderType(true, b, i);
    }
    else {
        glaciersWidget->setRenderType(false, b, i);
    }

    // Texture
    if (ui->radioButton_ShadingIceNeutral->isChecked())        glaciersWidget->setTextureType(GlaciersWidget::GlacierTextureType::NONE);
    if (ui->radioButton_ShadingIceELA->isChecked())            glaciersWidget->setTextureType(GlaciersWidget::GlacierTextureType::ELA);
    if (ui->radioButton_ShadingIceThick->isChecked())          glaciersWidget->setTextureType(GlaciersWidget::GlacierTextureType::THICKNESS);
    if (ui->radioButton_ShadingIceGradient->isChecked())       glaciersWidget->setTextureType(GlaciersWidget::GlacierTextureType::GRADIENT);
    if (ui->radioButton_ShadingIceStress->isChecked())         glaciersWidget->setTextureType(GlaciersWidget::GlacierTextureType::STRESS);
    if (ui->radioButton_ShadingIceDiffusivity->isChecked())    glaciersWidget->setTextureType(GlaciersWidget::GlacierTextureType::DIFFUSIVITY);
    if (ui->radioButton_ShadingIceDiffusivityRaw->isChecked()) glaciersWidget->setTextureType(GlaciersWidget::GlacierTextureType::DIFFUSIVITY_RAW);
    if (ui->radioButton_ShadingIceUDeform->isChecked())        glaciersWidget->setTextureType(GlaciersWidget::GlacierTextureType::SPEED_DEFORM);
    if (ui->radioButton_ShadingIceUSlip->isChecked())          glaciersWidget->setTextureType(GlaciersWidget::GlacierTextureType::SPEED_SLIP);
    if (ui->radioButton_ShadingIceUDominant->isChecked())      glaciersWidget->setTextureType(GlaciersWidget::GlacierTextureType::SPEED_DOMINANT_TYPE);
    if (ui->radioButton_ShadingIceDifference->isChecked())     glaciersWidget->setTextureType(GlaciersWidget::GlacierTextureType::ICE_DIFFERENCE);
    if (ui->radioButton_ShadingFeatShader->isChecked())        glaciersWidget->setTextureType(GlaciersWidget::GlacierTextureType::FEAT_SHADER);

    glaciersWidget->update();
}

void MainWindow::RenderUpdate()
{
    Render(false);
}

void MainWindow::editingSceneLeft(const Vector3& anchor)
{
    Vector2 ctr = Vector2(anchor[0], anchor[1]);
    if (ui->radioButton_brushBed->isChecked()) {
        glacierTerrain.GetBedrock().addGaussian(ctr, brushRadius, rockStrength);
        glacierTerrain.updateBedrockGPU();
        glaciersWidget->updateGeometry();
    }
    else if (ui->radioButton_brushIce->isChecked()) {
        glacierTerrain.GetIce().addGaussian(ctr, brushRadius, iceStrength);
        glacierTerrain.updateIceGPU();
        glaciersWidget->updateGeometry();
    }
    else if (ui->radioButton_brushELA->isChecked()) {
        glacierTerrain.GetELA().addGaussian(ctr, brushRadius, elaStrength);
        glacierTerrain.updateELAGPU();
		updatePreviewELA(glacierTerrain.GetELA());
    }
    else if (ui->radioButton_brushPrec->isChecked()) {
        glacierTerrain.GetAccumRate().addGaussian(ctr, brushRadius, precStrenth);
        glacierTerrain.updateAccumRateGPU();
		updatePreviewAccum(glacierTerrain.GetAccumRate());
    }

    Render(false);
}

void MainWindow::editingSceneRight(const Vector3& anchor)
{
	Vector2 ctr = Vector2(anchor[0], anchor[1]);
    if (ui->radioButton_brushBed->isChecked()) {
        glacierTerrain.GetBedrock().addGaussian(ctr, brushRadius, -rockStrength);
        glacierTerrain.updateBedrockGPU();
        glaciersWidget->updateGeometry();
    }
    else if (ui->radioButton_brushIce->isChecked()) {
        glacierTerrain.GetIce().addGaussian(ctr, brushRadius, -iceStrength);
        glacierTerrain.updateIceGPU();
        glaciersWidget->updateGeometry();
    }
    else if (ui->radioButton_brushELA->isChecked()) {
        glacierTerrain.GetELA().addGaussian(ctr, brushRadius, -elaStrength);
        glacierTerrain.updateELAGPU();
		updatePreviewELA(glacierTerrain.GetELA());
    }
    else if (ui->radioButton_brushPrec->isChecked()) {
        glacierTerrain.GetAccumRate().addGaussian(ctr, brushRadius, -precStrenth);
        glacierTerrain.updateAccumRateGPU();
		updatePreviewAccum(glacierTerrain.GetAccumRate());
    }

    Render(false);
}

void MainWindow::updateBrushRadius(double d)
{
	brushRadius = d;
	glaciersWidget->updateBrushRadius(d);
}

void MainWindow::configSimulation()
{
    double factorUdeform = std::min(2.0 - 0.1 * double(ui->sliderFactorU->value()), 1.0);
    double factorUslide = std::min(0.1 * double(ui->sliderFactorU->value()), 1.0);

    ScalarField2 configELA = ui->checkUseELAMap->isChecked() ?
        getVariableELA(ui->sb_ELA->value(), ui->sb_ELAdev->value()) :
        ScalarField2(glacierTerrain.getDomain(), glacierTerrain.numCellsX(), glacierTerrain.numCellsY(), ui->sb_ELA->value());

    ScalarField2 configBeta = ui->checkUsePrecMap->isChecked() ?
        factorBetaMaps[multiresIndex] :
        ScalarField2(glacierTerrain.getDomain(), glacierTerrain.numCellsX(), glacierTerrain.numCellsY(), 1.0);

    glacierTerrain.configSimulation(configELA, configBeta,
        ui->sb_AccumRate->value(), ui->sb_AblateRate->value(),
        factorUdeform, factorUslide);

	double dIceSteady = ui->sb_dIce->value();
	glaciersWidget->setSteadyCondition(dIceSteady);
    glaciersWidget->setMinimumSimulationYears(1);
}

void MainWindow::runGlacierSimulation()
{
    configSimulation();
    glaciersWidget->runSimulation(ui->checkPauseOnSteady->isChecked());
}

void MainWindow::pauseGlacierSimulation()
{
    glaciersWidget->pauseSimulation();
}

void MainWindow::resetGlacierSimulation()
{
    glacierTerrain.resetSimulation();
    if (ui->checkUseInitialIce->isChecked() && int(initialIceMaps.size()) > multiresIndex) {
        glacierTerrain.setIceMap(initialIceMaps[multiresIndex]);
    }
    glaciersWidget->pauseSimulation();
    glaciersWidget->resetSimulation();
    glaciersWidget->update();
}

void MainWindow::upscaleGlacier()
{
    if (multiresIndex <= 0) {
        return;
    }

    multiresIndex--;
    ScalarField2 hiresIce = glacierTerrain.remapIceSurface(bedrocksMultires[multiresIndex]);
    glacierTerrain = GlacierTerrain(bedrocksMultires[multiresIndex], hiresIce);

	double dIce = glaciersWidget->getLastdIce();
	if (bedrocksMultires.size() > 0) {
		double r = double(bedrocksMultires[multiresIndex].getSizeX()) / double(bedrocksMultires[multiresIndex+1].getSizeX());
		dIce *= r;
	}
	ui->sb_dIce->setValue(dIce);

    glaciersWidget->setGlacierTerrain(&glacierTerrain);
    configSimulation();

	ui->labelResolution->setText(QString::number(hiresIce.getSizeX()) + "x" + QString::number(hiresIce.getSizeY()) + (multiresIndex > 0 ? "" : " (max)"));

    Render(false);   
}

void MainWindow::updateSteadyCondition(double d)
{
	glaciersWidget->setSteadyCondition(d);
}

void MainWindow::saveScene()
{
    QString filename = QFileDialog::getSaveFileName(nullptr,
      "Choose a filename to save",
      "",
      "png file (*.png)");

    if (!filename.isNull()) {

        double iceMin, iceMax;
        glacierTerrain.GetIce().getRange(iceMin, iceMax);
		QImage imageIce = glacierTerrain.GetIce().CreateImage(0.0, 256.0 * 256.0, false);
		imageIce.save(QString(filename).replace(".png", "_ice_" + QString::number(iceMax) + ".png"));

		double hmin, hmax;
		ScalarField2 hf = glacierTerrain.GetHeightfield();
		hf.getRange(hmin, hmax);
		hf.CreateImage(0.0, hmax, false).mirrored(false, true).save(QString(filename).replace(".png", "_heights_" + QString::number(hmax) + ".png"));

        QImage aboveELA(imageIce.width(), imageIce.height(), QImage::Format::Format_ARGB32);
        for (int i = 0; i < aboveELA.width(); i++) {
            for (int j = 0; j < aboveELA.height(); j++) {
                aboveELA.setPixel(i, j, glacierTerrain.aboveELA(i,j) ? QColor(255,255,255).rgba() : QColor(0,0,0,255).rgba() );
            }
        }
        aboveELA.save(QString(filename).replace(".png", "_ela.png"));
    }
}

void MainWindow::loadBedrocks(const QString& filename, const std::vector<int>& downscales, double terrainKM, const std::vector<double>& elevRange)
{
	QImage bedImg = QImage(filename).mirrored(false, true);

	int w = bedImg.width();
	int h = bedImg.height();
	double terrainWidth  = w >= h ? terrainKM : terrainKM * (double(w) / double(h));
	double terrainHeight = h >= w ? terrainKM : terrainKM * (double(h) / double(w));
	Box2 terrainSize(Vector2(0), Vector2(terrainWidth*1000.0, terrainHeight*1000.0));

	ScalarField2 bedrock(terrainSize, bedImg, 0, 0.1 * (256 * 256 - 1), true);
	double hmin, hmax;
	bedrock.getRange(hmin, hmax);
	if (elevRange.size() == 1) {
		bedrock *= elevRange[0];
	}
	else if (elevRange.size() == 2) {
		for (int i = 0; i < bedrock.getSizeX() * bedrock.getSizeY(); i++) {
			bedrock[i] = elevRange[0] + (elevRange[1] - elevRange[0]) * (bedrock[i] - hmin) / (hmax - hmin);
		}
		hmin = elevRange[0];
		hmax = elevRange[1];
	}

	bedrocksMultires.clear();
	for (int s : downscales) {
		bedrocksMultires.push_back(bedrock.setResolution(bedrock.getSizeX() / s, bedrock.getSizeY() / s));
	}
	multiresIndex = int(bedrocksMultires.size()) - 1;

	ui->sb_dIce->setValue(0.001);

	glacierTerrain = GlacierTerrain(bedrocksMultires.back());
	glaciersWidget->setGlacierTerrain(&glacierTerrain);

	Render(true);

	demFilename = filename;
	ui->btnSelectDEM->setText(QFileInfo(filename).fileName());
	ui->labelResolution->setText(QString::number(bedrocksMultires.back().getSizeX())
		+ "x" + QString::number(bedrocksMultires.back().getSizeY()));
	ui->sb_DEM_km->setValue(terrainKM);
	ui->sb_DEM_hmin->setValue(hmin);
	ui->sb_DEM_hmax->setValue(hmax);
	QString ss = QString::number(downscales[0]);
	for (int i = 1; i < downscales.size(); i++) ss += "," + QString::number(downscales[i]);
	ui->le_DEM_scales->setText(ss);

	int s = std::min(ui->imglabelDEM->width(), ui->imglabelDEM->height());
	ui->imglabelDEM->setPixmap(QPixmap::fromImage(bedImg.mirrored(false, true)).scaled(s, s, Qt::KeepAspectRatio));
}

void MainWindow::loadMapsELA(const QString& filename)
{
    QImage sunImg = QImage(filename).mirrored(false, true);
    double sunavg = 0.0;
    for (int i = 0; i < sunImg.width(); i++) {
        for (int j = 0; j < sunImg.height(); j++) {
            sunavg += double(qRed(sunImg.pixel(i, j))) / 255.0;
        }
    }
    sunavg /= sunImg.width() * sunImg.height();

    for (unsigned int i = 0; i < bedrocksMultires.size(); i++) {
        const ScalarField2& bedrock = bedrocksMultires[i];
        int nx = bedrock.getSizeX();
        int ny = bedrock.getSizeY();
        QImage sunlight = sunImg.scaled(QSize(nx, ny), Qt::KeepAspectRatio, Qt::SmoothTransformation);
        ScalarField2 elaFactor(bedrock.getDomain(), nx, ny, 1.0);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                elaFactor(i, j) = double(qRed(sunlight.pixel(i, j))) / 255.0 - sunavg;
            }
        }
        factorElaDeviation.push_back(elaFactor);
    }

	updatePreviewELA(getVariableELA(ui->sb_ELA->value(), ui->sb_ELAdev->value()));
}

void MainWindow::loadMapsPrecipitation(const QString& filename)
{
    QImage precImg = QImage(filename).mirrored(false, true);

    for (unsigned int i = 0; i < bedrocksMultires.size(); i++) {
        const ScalarField2& bedrock = bedrocksMultires[i];
        int nx = bedrock.getSizeX();
        int ny = bedrock.getSizeY();
        QImage precipitation = precImg.scaled(QSize(nx, ny), Qt::KeepAspectRatio, Qt::SmoothTransformation);
        ScalarField2 precFactor(bedrock.getDomain(), ny, ny, 1.0);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                precFactor(i, j) = double(qRed(precipitation.pixel(i, j))) / 255.0;
            }
        }
        factorBetaMaps.push_back(precFactor);
    }

	updatePreviewAccum(factorBetaMaps.back());
}

void MainWindow::loadMapsIce(const QString& filename)
{
    QImage iceImg = QImage(filename); // .mirrored(false, true);

    for (unsigned int i = 0; i < bedrocksMultires.size(); i++) {
        const ScalarField2& bedrock = bedrocksMultires[i];
        int nx = bedrock.getSizeX();
        int ny = bedrock.getSizeY();
        QImage ice = iceImg.scaled(QSize(nx, ny), Qt::KeepAspectRatio, Qt::SmoothTransformation);
        ScalarField2 iceMap(bedrock.getDomain(), ice, 0, 256.0 * 256.0 - 1, false);
        initialIceMaps.push_back(iceMap);
    }

	ui->imglabelInitialIce->setPixmap(QPixmap::fromImage(iceImg.mirrored(false, true))
		.scaled(ui->imglabelELA->width(), ui->imglabelELA->height(), Qt::KeepAspectRatio));
}

ScalarField2 MainWindow::getVariableELA(double avgELA, double devELA)
{
    ScalarField2 elamap(factorElaDeviation[multiresIndex]);
    elamap *= devELA;
    elamap += avgELA;
    return elamap;
}

void MainWindow::updatePreviewELA(const ScalarField2& ela)
{
    ui->imglabelELA->setPixmap(QPixmap::fromImage(ela.CreateImage(ColorPalette::CoolWarm()).mirrored(false, true))
        .scaled(ui->imglabelELA->width(), ui->imglabelELA->height(), Qt::KeepAspectRatio));
}

void MainWindow::updatePreviewAccum(const ScalarField2& accum)
{
	ui->imglabelPrec->setPixmap(QPixmap::fromImage(accum.CreateImage(ColorPalette::CoolWarm()).mirrored(false, true))
		.scaled(ui->imglabelPrec->width(), ui->imglabelPrec->height(), Qt::KeepAspectRatio));
}
