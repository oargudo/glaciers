#ifndef VECTORFIELD2_H
#define VECTORFIELD2_H

#include <QVector>
#include <QImage>
#include "core.h"


class VectorField2
{
public:
	VectorField2();
	VectorField2(const Box2& domain, int nx, int ny, const Vector2& v = Vector2(0));

    int getSizeX() const { return nx; }
    int getSizeY() const { return ny; }
    Vector2 getCellSize() const { return cellSize; }

    Box2 getDomain() const { return domain; };
    void getRange(Vector2& vmin, Vector2& vmax) const;

    int vertexIndex(int i, int j) const { return i + nx*j; }
    bool validIndex(int i, int j) const { return (i >= 0) && (i < nx - 1) && (j >= 0) && (j < ny - 1); }

	void cellCoords(const Vector2& p, int& i, int& j, double& u, double& v) const;
    Vector2 domainCoords(int i, int j) const { return domain.getMin() + Vector2(i * cellSize[0], j * cellSize[1]); }

	Vector2 at(int i, int j) const { return field.at(vertexIndex(i, j)); }
	Vector2& operator()(int i, int j) { return field[vertexIndex(i, j)]; }
	Vector2 at(int i) const { return field.at(i); }
	Vector2& operator[](int i) { return field[i]; }
	Vector2 value(const Vector2& p) const;

protected:

    QVector<Vector2> field;
    int nx, ny;
    Box2 domain;
    Vector2 cellSize;

};

#endif // VECTORFIELD2_H
