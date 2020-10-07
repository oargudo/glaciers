#include "vectorfield2.h"
#include <QPoint>

VectorField2::VectorField2() : nx(0), ny(0), domain(Box2(Vector2(0), Vector2(1)))
{
}

VectorField2::VectorField2(const Box2 &domain, int nx, int ny, const Vector2& v) : nx(nx), ny(ny), domain(domain)
{
    field.fill(v, nx*ny);
    cellSize = Vector2(domain.width()/(nx-1), domain.height()/(ny-1));
}

void VectorField2::getRange(Vector2& vmin, Vector2& vmax) const
{
    vmin = vmax = field.at(0);
    for (int i = 1; i < field.size(); i++) {
        double x = field.at(i)[0];
        if (x < vmin[0]) vmin[0] = x;
        if (x > vmax[0]) vmax[0] = x;
		double y = field.at(i)[1];
		if (y < vmin[1]) vmin[1] = y;
		if (y > vmax[1]) vmax[1] = y;
    }
}

void VectorField2::cellCoords(const Vector2& p, int& i, int& j, double& u, double& v) const
{
	Vector2 q = p - domain.getMin();
	u = q[0] / cellSize[0];
	v = q[1] / cellSize[1];

	// Integer coordinates
	i = int(u);
	j = int(v);

	// Local coordinates within cell
	u -= i;
	v -= j;
}

Vector2 VectorField2::value(const Vector2& p) const
{
    double u, v;
    int i, j;
    cellCoords(p, i, j, u, v);

    if (!validIndex(i, j)) return Vector2(0.0);

	Vector2 a00 = at(i, j);
	Vector2 a01 = at(i, j + 1);
	Vector2 a10 = at(i + 1, j);
	Vector2 a11 = at(i + 1, j + 1);
	double x = Math::Bilinear(a00[0], a10[0], a11[0], a01[0], u, v);
	double y = Math::Bilinear(a00[1], a10[1], a11[1], a01[1], u, v);
	return Vector2(x, y);
}
