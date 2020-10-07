#include "scalarfield2.h"
#include <QPoint>

ScalarField2::ScalarField2() : nx(0), ny(0), domain(Box2(Vector2(0), Vector2(1)))
{
}

ScalarField2::ScalarField2(const Box2 &domain, int nx, int ny, double v) : nx(nx), ny(ny), domain(domain)
{
    field.fill(v, nx*ny);
    cellSize = Vector2(domain.width()/(nx-1), domain.height()/(ny-1));
}

ScalarField2::ScalarField2(const Box2& box, const QImage& image, const double& a, const double& b, bool grayscale)
{
    domain = box;
    nx = image.width();
    ny = image.height();
    field.resize(nx*ny);
	cellSize = Vector2(domain.width()/(nx - 1), domain.height()/(ny - 1));

    // Write Heightmap
    for (int i = 0; i < image.width(); i++) {
        for (int j = 0; j < image.height(); j++) {
            double t = 0.0;

            // Grayscale
            if (grayscale) {
                // Grayscale 16 bits
                if (image.format() == QImage::Format_Grayscale16) {
                    QColor thecolor = image.pixelColor(i, j);
                    t = thecolor.blueF();
                }
                // Grayscale 8 bits
                else {
                    QRgb color = image.pixel(i, j);
                    t = double(qGray(color)) / 255.0;
                }
            }
            // Color
            else {
                QRgb color = image.pixel(i, j);
                // Maximum value is 256^3-1
                t = double(qRed(color) << 16 | qGreen(color) << 8 | qBlue(color)) / (16777216.0 - 1.0);
            }

            field[vertexIndex(i, j)] = Math::Lerp(a, b, t);
        }
    }
}

void ScalarField2::getRange(double &vmin, double &vmax) const
{
    vmin = vmax = field.at(0);
    for (int i = 1; i < field.size(); i++) {
        double x = field.at(i);
        if (x < vmin) vmin = x;
        if (x > vmax) vmax = x;
    }
}

void ScalarField2::cellCoords(const Vector2& p, int& i, int& j, double& u, double& v) const
{
    Vector2 q = p - domain.getMin();
    u = q[0]/cellSize[0];
    v = q[1]/cellSize[1];

    // Integer coordinates
    i = int(u);
    j = int(v);

    // Local coordinates within cell
    u -= i;
    v -= j;
}

void ScalarField2::cellIntegerCoords(const Vector2 &p, int &i, int &j) const
{
    Vector2 q = p - domain.getMin();
    i = int(q[0]/cellSize[0]);
    j = int(q[1]/cellSize[1]);
}

double ScalarField2::value(const Vector2& p) const
{
    double u, v;
    int i, j;
    cellCoords(p, i, j, u, v);

    if (!validIndex(i, j)) return 0.0;
    return Math::Bilinear(at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1), u, v);
}

Vector3 ScalarField2::vertex(int i, int j) const
{
     return Vector3(domain.getMin()[0] + i*cellSize[0], domain.getMin()[1] + j*cellSize[1], at(i, j));
}

inline Vector3 triangleAreaNormal(const Vector3& p0, const Vector3& p1, const Vector3& p2) {
    return 0.5*cross(p1 - p0, p2 - p0);
}

Vector3 ScalarField2::normal(int i, int j) const
{
    Vector3 n;
    if (i == 0) {
        if (j == 0) {
            // Corner: 0/1
            n = triangleAreaNormal(vertex(i, j), vertex(i + 1, j), vertex(i + 1, j + 1))
              + triangleAreaNormal(vertex(i, j), vertex(i + 1, j + 1), vertex(i, j + 1));
        }
        else if (j == ny - 1) {
            // Corner: 5
            n = triangleAreaNormal(vertex(i, j), vertex(i, j - 1), vertex(i + 1, j));
        }
        else {
            // Edge: 0/1/5
            n = triangleAreaNormal(vertex(i, j), vertex(i + 1, j), vertex(i + 1, j + 1))
              + triangleAreaNormal(vertex(i, j), vertex(i + 1, j + 1), vertex(i, j + 1))
              + triangleAreaNormal(vertex(i, j), vertex(i, j - 1), vertex(i + 1, j));
        }
    }
    else if (i == nx - 1) {
        if (j == 0) {
            // Corner: 2
            n = triangleAreaNormal(vertex(i, j), vertex(i, j + 1), vertex(i - 1, j));
        }
        else if (j == ny - 1) {
            // Corner: 3/4
            n = triangleAreaNormal(vertex(i, j), vertex(i - 1, j - 1), vertex(i, j - 1))
              + triangleAreaNormal(vertex(i, j), vertex(i - 1, j), vertex(i - 1, j - 1));
        } else {
            // Edge: 2/3/4
            n = triangleAreaNormal(vertex(i, j), vertex(i, j + 1), vertex(i - 1, j))
              + triangleAreaNormal(vertex(i, j), vertex(i - 1, j), vertex(i - 1, j - 1))
              + triangleAreaNormal(vertex(i, j), vertex(i - 1, j - 1), vertex(i, j - 1));
      }
    }
    else {
        if (j == 0) {
            // Edge: 0/1/2
            n = triangleAreaNormal(vertex(i, j), vertex(i + 1, j), vertex(i + 1, j + 1))
              + triangleAreaNormal(vertex(i, j), vertex(i + 1, j + 1), vertex(i, j + 1))
              + triangleAreaNormal(vertex(i, j), vertex(i, j + 1), vertex(i - 1, j));
        }
        else if (j == ny - 1) {
            // Edge: 3/4/5
            n = triangleAreaNormal(vertex(i, j), vertex(i - 1, j), vertex(i - 1, j - 1))
              + triangleAreaNormal(vertex(i, j), vertex(i - 1, j - 1), vertex(i, j - 1))
              + triangleAreaNormal(vertex(i, j), vertex(i, j - 1), vertex(i + 1, j));
        }
        else {
            // Face: 0/1/2/3/4/5
            n = triangleAreaNormal(vertex(i, j), vertex(i + 1, j), vertex(i + 1, j + 1))
              + triangleAreaNormal(vertex(i, j), vertex(i + 1, j + 1), vertex(i, j + 1))
              + triangleAreaNormal(vertex(i, j), vertex(i, j + 1), vertex(i - 1, j))
              + triangleAreaNormal(vertex(i, j), vertex(i - 1, j), vertex(i - 1, j - 1))
              + triangleAreaNormal(vertex(i, j), vertex(i - 1, j - 1), vertex(i, j - 1))
              + triangleAreaNormal(vertex(i, j), vertex(i, j - 1), vertex(i + 1, j));
        }
    }

    return Normalized(n);
}

Vector2 ScalarField2::gradient(int i, int j) const
{
	Vector2 n;

	// Gradient along x axis
	if (i == 0) 
		n[0] = (at(i + 1, j) - at(i, j)) / cellSize[0];
	else if (i == nx - 1) 
		n[0] = (at(i, j) - at(i - 1, j)) / cellSize[0];
	else 
		n[0] = (at(i + 1, j) - at(i - 1, j)) / (2.0 * cellSize[0]);

	// Gradient along y axis
	if (j == 0)
		n[1] = (at(i, j + 1) - at(i, j)) / cellSize[1];
	else if (j == ny - 1)
		n[1] = (at(i, j) - at(i, j - 1)) / cellSize[1];
	else
		n[1] = (at(i, j + 1) - at(i, j - 1)) / (2.0 * cellSize[1]);

	return n;
}

VectorField2 ScalarField2::gradientField() const
{
	VectorField2 v(getDomain(), nx, ny);
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			v(i, j) = gradient(i, j);
		}
	}
	return v;
}

void ScalarField2::fill(double d)
{
    field.fill(d, nx*ny);
}

void ScalarField2::smooth(int n)
{
    // Smooth the scalar field using a discrete gaussian kernel.
    // The function uses a 3^2 approximation of the Gaussian kernel.
    QVector<double> smoothed;
    smoothed.resize(nx * ny);
    int k;

    for (int iter = 0; iter < n; iter++) {
        // Smooth center
        for (int i = 1; i < nx - 1; i++) {
            for (int j = 1; j < ny - 1; j++) {
                k = vertexIndex(i, j);
                smoothed[k] = (4.0*at(k) + 2.0*at(k - 1) + 2.0*at(k + 1) + 2.0*at(k - nx) + 2.0*at(k + nx) + at(k - 1 - nx) + at(k + 1 - nx) + at(k - 1 + nx) + at(k + 1 + nx)) / 16.0;
            }
        }

        // Smooth edges
        for (int i = 1; i < nx - 1; i++) {
            k = vertexIndex(i, 0);
            smoothed[k] = (4.0*at(k) + 2.0*at(k - 1) + 2.0*at(k + 1) + 2.0*at(k + nx) + at(k - 1 + nx) + at(k + 1 + nx)) / 12.0;
            k = vertexIndex(i, ny - 1);
            smoothed[k] = (4.0*at(k) + 2.0*at(k - 1) + 2.0*at(k + 1) + 2.0*at(k - nx) + at(k - 1 - nx) + at(k + 1 - nx)) / 12.0;
        }
        for (int j = 1; j < ny - 1; j++) {
            k = vertexIndex(0, j);
            smoothed[k] = (2.0*at(k - nx) + 4.0*at(k) + 2.0*at(k + nx) + at(k + 1 - nx) + 2.0*at(k + 1) + at(k + 1 + nx)) / 12.0;
            k = vertexIndex(nx - 1, j);
            smoothed[k] = (2.0*at(k - nx) + 4.0*at(k) + 2.0*at(k + nx) + at(k - 1 - nx) + 2.0*at(k - 1) + at(k - 1 + nx)) / 12.0;
        }

        // Corners
        k = vertexIndex(0, 0);
        smoothed[k] = (4.0*at(k) + 2.0*at(k + 1) + 2.0*at(k + nx) + 1.0*at(k + nx + 1)) / 9.0;
        k = vertexIndex(nx - 1, 0);
        smoothed[k] = (4.0*at(k) + 2.0*at(k - 1) + 2.0*at(k + nx) + 1.0*at(k + nx - 1)) / 9.0;
        k = vertexIndex(0, ny - 1);
        smoothed[k] = (4.0*at(k) + 2.0*at(k + 1) + 2.0*at(k - nx) + 1.0*at(k - nx + 1)) / 9.0;
        k = vertexIndex(nx - 1, ny - 1);
        smoothed[k] = (4.0*at(k) + 2.0*at(k - 1) + 2.0*at(k - nx) + 1.0*at(k - nx - 1)) / 9.0;
    }
    field = smoothed;
}

void ScalarField2::gaussianBlur()
{
	const int kernelSize = 9;
	double kernel[kernelSize] = { 1, 8, 28, 56, 70, 56, 28, 8, 1 };

	ScalarField2 temp(*this);
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			double v = 0;
			double w = 0;
			for (int dj = -kernelSize / 2; dj <= kernelSize / 2; dj++) {
				if (j + dj >= 0 && j + dj < ny - 1) {
					v += kernel[dj + kernelSize / 2] * at(i, j + dj);
					w += kernel[dj + kernelSize / 2];
				}
			}
			temp(i, j) = v / w;
		}
	}
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			double v = 0;
			double w = 0;
			for (int di = -kernelSize / 2; di <= kernelSize / 2; di++) {
				if (i + di >= 0 && i + di < nx - 1) {
					v += kernel[di + kernelSize / 2] * temp.at(i + di, j);
					w += kernel[di + kernelSize / 2];
				}
			}
			(*this)(i, j) = v / w;
		}
	}
}

void ScalarField2::step(const double& a, const double& b)
{
    for (int i = 0; i < field.size(); i++) {
        field[i] = Math::LinearStep(field.at(i), a, b);
    }
}

void ScalarField2::normalize()
{
	double a, b;
	getRange(a, b);
	if (a == b) {
		field.fill(1.0);
	}
	else {
		for (int i = 0; i < field.size(); i++) {
			field[i] = (field[i] - a) / (b - a);
		}
	}
}

void ScalarField2::addGaussian(const Vector2& center, const double& radius, const double& height)
{

    Box2 box (center, radius);

    int ia,ib,ja,jb;
    cellIntegerCoords(box.getMin(), ia, ja);
    cellIntegerCoords(box.getMax(), ib, jb);

    QPoint pa(ia, ja);
    QPoint pb(ib, jb);

    // Rectangle
    QRect area(pa.x(), pa.y(), pb.x() - pa.x() + 1, pb.y() - pa.y() + 1);

    // Limit to domain
    QRect mask(0, 0, nx - 1, ny - 1);
    area = area.intersected(mask);

    // Add to field
    for (int y = area.y(); y <= area.y() + area.height(); y++) {
        for (int x = area.x(); x <= area.x() + area.width(); x++) {
            // Distance between central point and current point
            double u = SquaredNorm(center - domainCoords(x, y));

            if (u < radius * radius) {
                field[vertexIndex(x, y)] += height * Math::CubicSmooth(u, radius * radius);
           }
        }
    }
}

ScalarField2 ScalarField2::setResolution(int x, int y) const
{
    // Sampled scalar field
    ScalarField2 sampled(domain, x, y);

    // Corners
    sampled(0, 0) = at(0, 0);
    sampled(0, y - 1) = at(0, ny - 1);
    sampled(x - 1, 0) = at(nx - 1, 0);
    sampled(x - 1, y - 1) = at(nx - 1, ny - 1);

    // Borders (use linear interpolation)
    for (int i = 1; i < x - 1; i++) {
        double tx = (nx - 1) * (i / double(x - 1));
        int x0 = int(floor(tx));
        int x1 = int(ceil(tx));

        sampled(i, 0) = Math::Lerp(at(x0, 0), at(x1, 0), tx - x0);
        sampled(i, y - 1) = Math::Lerp(at(x0, ny - 1), at(x1, ny - 1), tx - x0);
    }
    for (int j = 1; j < y - 1; j++) {
        double ty = (ny - 1) * (j / double(y - 1));
        int y0 = int(floor(ty));
        int y1 = int(ceil(ty));

        sampled(0, j) = Math::Lerp(at(0, y0), at(0, y1), ty - y0);
        sampled(x - 1, j) = Math::Lerp(at(nx - 1, y0), at(nx - 1, y1), ty - y0);
    }

    // Interior
    for (int i = 1; i < x - 1; i++) {
        for (int j = 1; j < y - 1; j++) {
            sampled(i, j) = value(sampled.domainCoords(i, j));
        }
    }

    return sampled;
}

QImage ScalarField2::CreateImage(bool grayscale) const
{
    double a, b;
    this->getRange(a, b);
    if (a == b) {
        b = a + 1.0;
    }
    return CreateImage(a, b, grayscale);
}

QImage ScalarField2::CreateImage(double a, double b, bool grayscale) const
{
    QImage image(nx, ny, QImage::Format_ARGB32);
    for (int i = 0; i < image.width(); i++) {
        for (int j = 0; j < image.height(); j++) {
            double x = field.at(vertexIndex(i, j));
            double y = Math::LinearStep(x, a, b);

            QColor color;
            if (grayscale) {
                int c = int(y * 255.0);
                color = QColor(c, c, c);
            }
            else {
                int c = int(y * (256.0 * 256.0 * 256.0 - 1.0));
                int cr = (c >> 16) & 255;
                int cv = (c >> 8) & 255;
                int cb = c & 255;
                color = QColor(cr, cv, cb);
            }
            image.setPixel(i, j, color.rgb());
        }
    }
    return image;
}

QImage ScalarField2::CreateImage(const ColorPalette& palette) const
{
    double a, b;
    getRange(a, b);
    if (a == b) {
        b = a + 1.0;
    }
    return CreateImage(a, b, palette);
}

QImage ScalarField2::CreateImage(double a, double b, const ColorPalette& palette) const
{
    QImage image(nx, ny, QImage::Format_ARGB32);
    for (int i = 0; i < image.width(); i++) {
        for (int j = 0; j < image.height(); j++) {
            double x = field.at(vertexIndex(i, j));
            double y = Math::LinearStep(x, a, b);
            QColor color = palette.getColor(y).toQColor();
            image.setPixel(i, j, color.rgb());
        }
    }
    return image;
}

ScalarField2& ScalarField2::operator+=(const ScalarField2& s)
{
    for (int i = 0; i < field.size(); i++) {
        field[i] += s.at(i);
    }
    return *this;
}

ScalarField2 &ScalarField2::operator+=(const double &d)
{
    for (int i = 0; i < field.size(); i++) {
        field[i] += d;
    }
    return *this;
}

ScalarField2& ScalarField2::operator*=(const double &d)
{
    for (int i = 0; i < field.size(); i++) {
        field[i] *= d;
    }
    return *this;
}

ScalarField2 operator+(const ScalarField2& s1, const ScalarField2& s2)
{
    ScalarField2 r(s1);
    for (int i = 0; i < r.field.size(); i++) {
        r.field[i] += s2.at(i);
    }
    return r;
}
