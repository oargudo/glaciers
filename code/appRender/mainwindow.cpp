#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "vectorfield2.h"
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
    features = nullptr;

    createActions();
}

MainWindow::~MainWindow()
{
    if (features) delete features;
    delete ui;
}

void MainWindow::createActions()
{
	connect(ui->btnOpenFolder, SIGNAL(clicked()), this, SLOT(openFolder()));
	connect(ui->btnLoadTerrain, SIGNAL(clicked()), this, SLOT(loadTerrain()));
	connect(ui->btnUpdateFeatures, SIGNAL(clicked()), this, SLOT(updateFeatures()));

    connect(&featuresThread, SIGNAL(textureUpdate(const QImage&)), this, SLOT(updateTexture(const QImage&)));
    connect(&featuresThread, SIGNAL(statusUpdate(const QString&)), this, SLOT(setStatusText(const QString&)));
    connect(&featuresThread, SIGNAL(finished()), this, SLOT(finishedFeatures()));
}

void MainWindow::openFolder()
{
	QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
					"", QFileDialog::ShowDirsOnly);
	if (!dir.isEmpty()) {
		ui->mapsFolder->setText(dir);
	}
}

void MainWindow::loadTerrain()
{
	QString dir = QDir::cleanPath(ui->mapsFolder->text()) + QDir::separator();
	QString demFile = dir + "dem.png";
	QString hiFile  = dir + "dem-hires.png";
	QString iceFile  = dir + "ice.png";
	QString elaFile = dir + "ela.png";
	
    if (QFile::exists(demFile) && QFile::exists(iceFile) && QFile::exists(elaFile)) {

        statusBar()->showMessage("Loading terrain...");

        if (features) delete features;
        features = new FeaturesPainter();
        features->TEXTURE_RES = ui->sbTexResolution->value();

		double tw = ui->sbTerrainWidth->value();
		double th = ui->sbTerrainHeight->value();
        features->terrainBox = Box2(Vector2(0), Vector2(tw*1000.0, th*1000.0));

		QImage imgDemLow(demFile);
        QImage imgDemHi;
        if (QFile::exists(hiFile))
            imgDemHi.load(hiFile);
        else
            imgDemHi.load(demFile);

		// low res
        features->bedLow = ScalarField2(features->terrainBox, imgDemLow.mirrored(false, true), 0, 0.1*(256*256-1), true);
        features->iceLow = ScalarField2(features->terrainBox, QImage(iceFile), 0.0, (256 * 256) - 1.0, false);
        GlacierTerrain terrainLow(features->bedLow, features->iceLow);
        features->hfLow = features->bedLow + features->iceLow;

		// hi res
        features->bedHi = ScalarField2(features->terrainBox, imgDemHi.mirrored(false, true), 0, 0.1 * (256 * 256 - 1), true);
        features->iceHi = terrainLow.remapIceSurface(features->bedHi);
        features->hfHi = features->bedHi + features->iceHi;

		// above ELA texture
        features->accumArea = ScalarField2(features->terrainBox, QImage(elaFile), 0.0, 1.0, true);
        features->accumArea.smooth(1);
        features->accumArea = features->accumArea.setResolution(features->hfHi.getSizeX(), features->hfHi.getSizeY());

        // create base texture
        features->imgBase = features->baseTexture(features->hfHi, features->iceHi, features->accumArea);

		// update viewer
        glaciersWidget->setHeightfield(features->hfHi);
        glaciersWidget->setTexture(features->imgBase);
		glaciersWidget->update();

        statusBar()->showMessage("Terrain loaded");
	}
	else {
        statusBar()->showMessage("ERROR: Incorrect folder or missing map files");
	}
}

void MainWindow::updateFeatures()
{
    ui->btnUpdateFeatures->setEnabled(false);
    ui->btnUpdateFeatures->setText("Placing features...");
    statusBar()->showMessage("Placing features...");

    featuresThread.setUI(ui);
    featuresThread.setPainter(features);
    featuresThread.start();
}

void MainWindow::finishedFeatures()
{
    ui->btnUpdateFeatures->setText("Update features");
    ui->btnUpdateFeatures->setEnabled(true);
    statusBar()->showMessage("Done!", 5000);
}

void MainWindow::updateTexture(const QImage &img)
{
    // update viewer
    glaciersWidget->setTexture(img.mirrored(false, true));
    glaciersWidget->update();
}

void MainWindow::setStatusText(const QString& s)
{
    statusBar()->showMessage(s);
}

void FeaturesWorker::run()
{
    emit statusUpdate("Placing features...");

	// Final texture
    QImage imgFeatures(features->imgBase.mirrored(false, true));

	// declarations
	QString dir = QDir::cleanPath(ui->mapsFolder->text()) + QDir::separator();

    // replicability
    features->rne = std::minstd_rand(133742);

    // maps derived from heightfield, only need to compute once unless terrain reloaded
    if (features->flowDirField.getSizeX()*features->flowDirField.getSizeY() <= 0) {
        // load and compute additional maps that will be used during feature placement
        features->flowDirField = features->hfLow.gradientField();
        features->slopeField = ScalarField2(features->terrainBox, features->flowDirField.getSizeX(), features->flowDirField.getSizeY(), 0);
        for (int i = 0; i < features->flowDirField.getSizeX(); i++) {
            for (int j = 0; j < features->flowDirField.getSizeY(); j++) {
                double n = Norm(features->flowDirField.at(i, j));
                features->slopeField(i, j) = n;
                if (n > 1e-6) features->flowDirField(i, j) = -features->flowDirField.at(i, j) / n;
                else features->flowDirField(i, j) = Vector2(0, 0);
            }
        }
    }
    features->valleyDist = ScalarField2(features->terrainBox, QImage(dir + "valleydist.png"), 0.0, 256 * 256 - 1.0, false);
    features->glacierIds = QImage(dir + "segmentation.png");

	// alpha map for crevasses
    ScalarField2 alphaCrevasses(features->iceHi);
	alphaCrevasses.step(0, 2.0);
	alphaCrevasses = alphaCrevasses.setResolution(imgFeatures.width(), imgFeatures.height());

	// SERACS AND ICEFALLS: perturb heightfield
	if (ui->renderSeracs->isChecked()) {
        emit statusUpdate("Seracs...");

        ScalarField2 seracsField(features->terrainBox, QImage(dir + "seracs.png"), 0.0, 256 * 256 - 1.0, true);
        ScalarField2 icefalls = ScalarField2(features->terrainBox, QImage(dir + "icefallsRaw.png"), 0.0, 256 * 256 - 1.0, true);
		double seracDisp = ui->sbSeracDisp->value();

		for (int i = 0; i < seracsField.getSizeX() * seracsField.getSizeY(); i++) {
			seracsField[i] = std::min(1.0, seracsField[i] + icefalls[i]);
		}
        ScalarField2 iceSeracs = features->seracs(seracsField, seracDisp);

        ScalarField2 hfSeracs = features->bedHi + iceSeracs;
        //glaciersWidget->setHeightfield(hfSeracs, false);
        imgFeatures = features->baseTexture(hfSeracs, features->iceHi, features->accumArea).mirrored(false, true);
        emit textureUpdate(imgFeatures);
	}

	// OGIVES: ondulations, need to draw before moraines as they will cover its margins (which have interpolation artifacts on distance field)
	if (ui->renderOgives->isChecked()) {
        emit statusUpdate("Ogives...");

        ScalarField2 ogiveDistSrc = ScalarField2(features->terrainBox, QImage(dir + "ogivesDistSrc.png"), 0.0, (256. * 256.) - 1.0, false);
        ScalarField2 ogiveDistCtr = ScalarField2(features->terrainBox, QImage(dir + "ogivesDistCtr.png"), 0.0, (256. * 256.) - 1.0, false);
		QImage ogiveId = QImage(dir + "ogivesId.png");
		double scaleLen = ui->sbOgiveLength->value();
		double scaleFreq = ui->sbOgiveFreq->value();

        ScalarField2 ogiveStrength(features->terrainBox, imgFeatures.width(), imgFeatures.height(), 0);
        QImage ogiveMap = features->ogives(ogiveDistSrc, ogiveDistCtr, ogiveId, scaleLen, scaleFreq, imgFeatures, ogiveStrength);
        emit textureUpdate(imgFeatures);

		// prevent crevasses from forming on top of ogives (for better illustration)
        ogiveStrength = ogiveStrength.setResolution(features->iceLow.getSizeX(), features->iceLow.getSizeY());
        ogiveStrength.smooth(1);
		ogiveStrength.step(0.1, 0.4);
		for (int i = 0; i < alphaCrevasses.getSizeX(); i++) {
			for (int j = 0; j < alphaCrevasses.getSizeY(); j++) {
                Vector2 p = alphaCrevasses.domainCoords(i, j);
                alphaCrevasses(i, j) = std::max(0.0, std::min(1.0, alphaCrevasses.at(i, j) - ogiveStrength.value(p)));
			}
		}
	}

	// MORAINES: draw moraine paths
	if (ui->renderMoraines->isChecked()) {
        emit statusUpdate("Moraines...");

        //const Vector3 moraineColor(164 / 255.0, 151 / 255.0, 125 / 255.0);
        const Vector3 moraineColor(0.9 * 171 / 255.0, 0.9 * 167 / 255.0, 0.9 * 164 / 255.0);
		double scale = ui->sbMoraineScale->value();
        QImage moraineMap = features->moraines(moraineColor, scale, imgFeatures);
        emit textureUpdate(imgFeatures);

		// prevent crevasses from forming on top of moraines (for better illustration)
		for (int i = 0; i < alphaCrevasses.getSizeX(); i++) {
			for (int j = 0; j < alphaCrevasses.getSizeY(); j++) {
				alphaCrevasses(i, j) = std::max(0.0, std::min(1.0, alphaCrevasses.at(i, j) 
					- double(moraineMap.pixelColor(i, moraineMap.height() - 1 - j).red()) / 255.0));
			}
		}
	}

	// CREVASSES: texture bombing
	if (ui->renderCrevasses->isChecked()) {
		double modDensityCT = ui->crevDensityT->value() / 10.0;
		double modDensityCL = ui->crevDensityL->value() / 10.0;
		double modDensityCM = ui->crevDensityT->value() / 10.0;
		double crevasseScale = ui->sbCrevScale->value();
		ScalarField2 crevasseField;

        emit statusUpdate("Transverse crevasses...");
        crevasseField = ScalarField2(features->terrainBox, QImage(dir + "crevassesTrans.png"), 0.0, 256 * 256 - 1.0, true);
		crevasseField.normalize();
		crevasseField *= modDensityCT;
        for (int i = 0; i < crevasseField.getSizeX() * crevasseField.getSizeY(); i++) {
            // fewer crevasses where no much ice above ELA
            crevasseField[i] *= std::min(1.0, features->iceLow[i] / 10.0) * features->accumArea[i] + (1 - features->accumArea[i]);
		}
        QImage crevTransMap = features->crevasses(crevasseField, 0, crevasseScale, alphaCrevasses, imgFeatures);
        emit textureUpdate(imgFeatures);

        emit statusUpdate("Longitudinal crevasses...");
        crevasseField = ScalarField2(features->terrainBox, QImage(dir + "crevassesLong.png"), 0.0, 256 * 256 - 1.0, true);
		crevasseField.normalize();
		crevasseField *= modDensityCL;
        QImage crevLongMap = features->crevasses(crevasseField, 1, crevasseScale, alphaCrevasses, imgFeatures);
        emit textureUpdate(imgFeatures);

        emit statusUpdate("Marginal crevasses...");
        crevasseField = ScalarField2(features->terrainBox, QImage(dir + "crevassesMargin.png"), 0.0, 256 * 256 - 1.0, true);
		crevasseField.normalize();
		crevasseField *= modDensityCM;
        QImage crevMarginMap = features->crevasses(crevasseField, 2, crevasseScale, alphaCrevasses, imgFeatures);
        emit textureUpdate(imgFeatures);
	}

	// ICEFALLS: as crevasses
	if (ui->renderIcefalls->isChecked()) {
        emit statusUpdate("Icefalls...");
		double modDensityCI = ui->crevDensityI->value() / 10.0;
        ScalarField2 icefalls = ScalarField2(features->terrainBox, QImage(dir + "icefallsRaw.png"), 0.0, 256 * 256 - 1.0, true);
		icefalls.normalize();
		icefalls *= modDensityCI;
        QImage crevIcefallsMap = features->crevasses(icefalls, 3, 1.0, alphaCrevasses, imgFeatures);
        emit textureUpdate(imgFeatures);
	}

	// RIMAYES: fixed positions
	if (ui->renderRimayes->isChecked()) {
        emit statusUpdate("Rimayes...");
        ScalarField2 rimayeField = ScalarField2(features->terrainBox, QImage(dir + "rimayes.png"), 0.0, 256 * 256 - 1.0, true);
        QImage rimayeMap = features->rimayes(rimayeField, imgFeatures);
        emit textureUpdate(imgFeatures);
    }

    emit finished();
}
