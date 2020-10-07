#ifndef FEATURESPAINTER_H
#define FEATURESPAINTER_H

#include <QImage>
#include <random>
#include "core.h"
#include "scalarfield2.h"
#include "vectorfield2.h"


class FeaturesPainter {

public:
    QImage baseTexture(const ScalarField2& hf, const ScalarField2& ice, const ScalarField2& ela) const;

    QImage crevasses(const ScalarField2& crevasseField, int crevasseType, double scale, const ScalarField2& alphaCrevasses, QImage& texture);
    QImage rimayes(const ScalarField2& rimayeField, QImage& texture);
    QImage ogives(const ScalarField2& ogiveDistSrc, const ScalarField2& ogiveDistCtr, const QImage& ogiveId, double ogiveLenScale, double ogiveFreq,
        QImage& texture, ScalarField2& ogiveStrength);
    QImage moraines(const Vector3& moraineColor, double moraineScale, QImage& texture);
    ScalarField2 seracs(const ScalarField2& seracs, double displacement);

public:
    Box2 terrainBox;

    ScalarField2 bedLow, iceLow, hfLow;
    ScalarField2 bedHi, iceHi, hfHi;
    ScalarField2 accumArea;

    ScalarField2 slopeField;
    VectorField2 flowDirField;
    ScalarField2 valleyDist;
    QImage glacierIds;

    QImage imgBase;
    double TEXTURE_RES = 1.25; // m

    std::minstd_rand rne;
    std::uniform_real_distribution<double> uniformDist = std::uniform_real_distribution<double>(0.0, 1.0);

protected:
    double randUniform(double a = 0.0, double b = 1.0);

};

#endif // FEATURESPAINTER_H
