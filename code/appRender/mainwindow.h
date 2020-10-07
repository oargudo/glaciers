#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QThread>
#include "glacierswidget.h"
#include "glacierterrain.h"
#include "featurespainter.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE


class FeaturesWorker : public QThread
{
    Q_OBJECT
    void run() override;
signals:
    void textureUpdate(const QImage& img);
    void statusUpdate(const QString& s);
    void finished();
public:
    void setUI(Ui::MainWindow* ui) { this->ui = ui; }
    void setPainter(FeaturesPainter* painter) { features = painter; }
private:
    Ui::MainWindow* ui;
    FeaturesPainter* features;
};


class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

public slots:
	void openFolder();
	void loadTerrain();
	void updateFeatures();
    void finishedFeatures();

    void updateTexture(const QImage& img);
    void setStatusText(const QString& s);

private:
    void createActions();

private:
    Ui::MainWindow *ui;	
    GlaciersWidget *glaciersWidget;
    FeaturesWorker featuresThread;

    FeaturesPainter *features;
};
#endif // MAINWINDOW_H
