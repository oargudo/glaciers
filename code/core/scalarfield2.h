#ifndef SCALARFIELD2_H
#define SCALARFIELD2_H

#include <QVector>
#include <QImage>
#include "core.h"
#include "vectorfield2.h"


class ScalarField2
{
public:
    ScalarField2();
    ScalarField2(const Box2& domain, int nx, int ny, double v = 0.0);
    ScalarField2(const Box2& domain, const QImage& image, const double& a, const double& b, bool grayscale);

    int getSizeX() const { return nx; }
    int getSizeY() const { return ny; }
    Vector2 getCellSize() const { return cellSize; }

    Box2 getDomain() const { return domain; };
    void getRange(double& vmin, double& vmax) const;

    int vertexIndex(int i, int j) const { return i + nx*j; }
    bool validIndex(int i, int j) const { return (i >= 0) && (i < nx - 1) && (j >= 0) && (j < ny - 1); }

    void cellCoords(const Vector2& p, int& i, int& j, double& u, double& v) const;
    void cellIntegerCoords(const Vector2& p, int& i, int& j) const;
    Vector2 domainCoords(int i, int j) const { return domain.getMin() + Vector2(i * cellSize[0], j * cellSize[1]); }

    double at(int i, int j) const { return field.at(vertexIndex(i, j)); }
    double& operator()(int i, int j) { return field[vertexIndex(i, j)]; }
    double at(int i) const { return field.at(i); }
    double& operator[](int i) { return field[i]; }
    double value(const Vector2& p) const;

    Vector3 vertex(int i, int j) const;
	Vector3 normal(int i, int j) const;
	Vector2 gradient(int i, int j) const;

	VectorField2 gradientField() const;

    void fill(double d);
    void smooth(int n);
	void gaussianBlur();
    void step(const double& a, const double& b);
	void normalize();
    void addGaussian(const Vector2& center, const double& radius, const double& height);

    ScalarField2 setResolution(int nx, int ny) const;


    QImage CreateImage(bool grayscale) const;
    QImage CreateImage(double a, double b, bool grayscale) const;
    QImage CreateImage(const ColorPalette& palette) const;
    QImage CreateImage(double a, double b, const ColorPalette& palette) const;

    ScalarField2& operator+=(const ScalarField2& s);
    ScalarField2& operator+=(const double& d);
    ScalarField2& operator*=(const double& d);
    friend ScalarField2 operator+(const ScalarField2& s1, const ScalarField2& s2);

protected:

    QVector<double> field;
    int nx, ny;
    Box2 domain;
    Vector2 cellSize;

};



#endif // SCALARFIELD2_H
