#ifndef _CORE_H_
#define _CORE_H_

#include <cmath>
#include <vector>
#include <QColor>

namespace Math {
    const double Pi = 3.14159265358979323846;

    inline double min(double a, double b) {
        return (a < b ? a : b);
    }

    inline double max(double a, double b) {
        return (a > b ? a : b);
    }

    inline double DegreeToRadian(double a) {
        return a * Math::Pi / 180.0;
    }

    inline double RadianToDegree(double a) {
        return a * 180.0 / Math::Pi;
    }

    inline double LinearStep(double x, double a, double b) {
        if (x < a) return 0;
        if (x > b) return 1;
        return (x - a) / (b - a);
    }

    inline double Clamp(double x, double a = 0, double b = 1) {
        return (x < a ? a : (x > b ? b : x));
    }

    inline double Lerp(double a, double b, double t) {
        return a + t * (b - a);
    }

    inline double Bilinear(double a00, double a10, double a11, double a01, double u, double v) {
        return (1 - u)*(1 - v)*a00 + (1 - u)*(v)*a01 + (u)*(1 - v)*a10 + (u)*(v)*a11;
    }

    inline double CubicSmooth(double x, double r) {
        return (1.0 - x / r)*(1.0 - x / r)*(1.0 - x / r);
    }

    inline int Integer(double x) {
        return x > 0.0 ? int(x) : int(x) - 1;
    }

    inline double Ridge(const double& z, const double& r) {
        if (z < r) return z;
        else return 2.0 * r - z;
    }
}

class SimplexNoise {
public:
    SimplexNoise() {}
    ~SimplexNoise() {}
    double at(double x, double y) const;
protected:
    double dot(const int* g, const double& x, const double& y) const {
        return g[0] * x + g[1] * y;
    }
    static const int perm[512];   //!< Permutation table, 256 entries duplicated once to avoid modulo computations.
    static const int grad2[8][2]; //!< Array of gradients for 2D noise.
    static const double F2, G2;   //!< Unskew factors for 2D case.
};


class Vector2 {
protected:
    double c[2];

public:
    Vector2() {
        c[0] = c[1] = 0;
    }
    explicit Vector2(double d) {
        c[0] = c[1] = d;
    }
    explicit Vector2(double d0, double d1) {
        c[0] = d0;
        c[1] = d1;
    }

    double& operator[] (int i) { return c[i]; };
    double operator[] (int i) const { return c[i]; };

    Vector2 operator- () const { return Vector2(-c[0], -c[1]); };
    friend Vector2 operator+(const Vector2& u, const Vector2& v) { return Vector2(u[0]+v[0], u[1]+v[1]); };
    friend Vector2 operator-(const Vector2& u, const Vector2& v) { return Vector2(u[0]-v[0], u[1]-v[1]); };
    friend Vector2 operator*(const Vector2& u, double a) { return Vector2(u[0]*a, u[1]*a); }
    friend Vector2 operator*(double a, const Vector2& v) { return v * a; }
    friend Vector2 operator/(const Vector2& u, double a) { return Vector2(u[0]/a, u[1]/a); }

    friend double Norm(const Vector2& u) { return sqrt(u[0]*u[0] + u[1]*u[1]); }
    friend double SquaredNorm(const Vector2& u) { return u[0]*u[0] + u[1]*u[1]; }
    friend Vector2 Normalized(const Vector2& u) { return u/Norm(u); }
};


class Vector3 {
protected:
    double c[3];

public:
    Vector3() {
        c[0] = c[1] = c[2] = 0;
    }
    explicit Vector3(double d) {
        c[0] = c[1] = c[2] = d;
    };
    explicit Vector3(double d0, double d1, double d2) {
        c[0] = d0;
        c[1] = d1;
        c[2] = d2;
    };
    explicit Vector3(const Vector2& u) {
        c[0] = u[0];
        c[1] = u[1];
        c[2] = 0;
    }

    double& operator[] (int i) { return c[i]; };
    double operator[] (int i) const { return c[i]; };
    Vector3 operator-() const { return Vector3(-c[0], -c[1], -c[2]); }

    friend Vector3 operator+(const Vector3& u, const Vector3& v) { return Vector3(u[0]+v[0], u[1]+v[1], u[2]+v[2]); };
    friend Vector3 operator-(const Vector3& u, const Vector3& v) { return Vector3(u[0]-v[0], u[1]-v[1], u[2]-v[2]); };
    friend Vector3 operator*(const Vector3& u, double a) { return Vector3(u[0]*a, u[1]*a, u[2]*a); }
    friend Vector3 operator*(double a, const Vector3& v) { return v * a; }
    friend Vector3 operator/(const Vector3& u, double a) { return Vector3(u[0]/a, u[1]/a, u[2]/a); }

    friend double Norm(const Vector3& u) { return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]); }
    friend double SquaredNorm(const Vector3& u) { return u[0]*u[0] + u[1]*u[1] + u[2]*u[2]; }
    friend Vector3 Normalized(const Vector3& u) { return u/Norm(u); }

    friend double dot(const Vector3& u, const Vector3& v) { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; }
    friend Vector3 cross(const Vector3& u, const Vector3& v) {
        return Vector3(u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]);
    }

    static Vector3 fromQColor(const QColor& c) {
        return Vector3(c.red() / 255.0, c.green() / 255.0, c.blue() / 255.0);
    }
    QColor toQColor() const {
        return QColor(int(255.0 * Math::Clamp(c[0])),
                      int(255.0 * Math::Clamp(c[1])),
                      int(255.0 * Math::Clamp(c[2])),
                      255);
    }
};


class Box2 {
protected:
    Vector2 bmin; // min
    Vector2 bmax; // max

public:
    Box2() : bmin(0), bmax(0) {};
    Box2(const Vector2& pmin, const Vector2& pmax) : bmin(pmin), bmax(pmax) {}
    Box2(const Vector2& c, double r) : bmin(c - Vector2(r)), bmax(c + Vector2(r)) {}

    Vector2 getMin() const { return bmin; }
    Vector2 getMax() const { return bmax; }
    Vector2 center() const { return 0.5*(bmin + bmax); }
    double radius() const { return 0.5 * Norm(bmax - bmin); }
    double width() const { return bmax[0] - bmin[0]; }
    double height() const { return bmax[1] - bmin[1]; }

    bool Intersect(const Vector2& s0, const Vector2& s1, double& tmin, double& tmax);
};


class Box3 {
protected:
    Vector3 bmin; // min
    Vector3 bmax; // max

public:
    Box3() : bmin(0), bmax(0) {};
    Box3(const Vector3& pmin, const Vector3& pmax) : bmin(pmin), bmax(pmax) {}

    Vector3 getMin() const { return bmin; }
    Vector3 getMax() const { return bmax; }
    Vector3 center() const { return 0.5*(bmin + bmax); }
    double radius() const { return 0.5 * Norm(bmax - bmin); }
    double width() const { return bmax[0] - bmin[0]; }
    double height() const { return bmax[1] - bmin[1]; }
    double depth() const { return bmax[2] - bmin[2]; }
};


class Ray
{
protected:
    Vector3 p; // Origin of the ray.
    Vector3 d; // Direction.

public:
    Ray() {}
    explicit Ray(const Vector3& p, const Vector3& d) : p(p), d(d) {}

    Vector3 origin() const { return p; }
    Vector3 direction() const { return d; }

    Vector3 operator()(double t) const { return p + t * d; }

    Ray reflect(const Vector3& p, const Vector3& n) { return Ray(p, n - 2 * n * dot(d, n)); }
};


class Camera {
protected:
    Vector3 eye;      // Eye
    Vector3 at;       // Look at point
    Vector3 up;       // Up vector
    double cah;       // Camera aperture horizontal
    double cav;       // Camera aperture vertical
    double nearplane; // Near plane
    double farplane;  // Far plane
    double fl;        // Focal length

public:
    Camera();
    Camera(const Vector3& eye, const Vector3& at, const Vector3& up = Vector3(0,0,1), double near = 1.0, double far = 100000.0);

    Vector3 getEye() const { return eye; }
    Vector3 getAt() const { return at; }
    Vector3 getUp() const { return up; }
    Vector3 getViewDir() const { return Normalized(at - eye); }
    double getNearPlane() const { return nearplane; }
    double getFarPlane() const { return farplane; }
    double getAngleOfViewH(double, double) const;
    double getAngleOfViewV(double, double) const;

    void setAt(const Vector3& p) { at = p; up = Vector3(0,0,1); }
    void setEye(const Vector3& p) { eye = p; }
    void setPlanes(double n, double f) { nearplane = n; farplane = f;}

    // Move camera around eye
    void upDownRound(double a);
    void leftRightRound(double a);
    void backForth(double a, bool moveAt = false);

    // Move camera in a plane
    void upDownPlane(double);
    void leftRightPlane(double);

    Ray pixelToRay(int px, int py, int w, int h) const;

    static Camera View(const Box3& box);
};


class ColorPalette {

protected:
    std::vector<Vector3> colors;
    std::vector<double> anchors;

public:
    ColorPalette() : colors({Vector3(1)}), anchors({0}) {}
    ColorPalette(const std::vector<Vector3>& c, const std::vector<double>& a) : colors(c), anchors(a) {}

    Vector3 getColor(double u) const;

    static ColorPalette CoolWarm() {
        const Vector3 Cool = Vector3(97, 130, 234) / 255.0;
        const Vector3 White = Vector3(221, 220, 219) / 255.0;
        const Vector3 Warm = Vector3(220, 94, 75) / 255.0;
        return ColorPalette({Cool, White, Warm}, {0, 0.5, 1});
    }
    static ColorPalette Relief() {
        const std::vector<Vector3> c = { Vector3(160, 220, 105) / 255.0,
                                         Vector3(1.0, 0.9, 0.45),
                                         Vector3(168/ 255.0, 155 / 255.0, 138 / 255.0),
                                         Vector3(0.95, 0.95, 0.95) };
        return ColorPalette(c, {0, 150, 250, 400});
    }
};


#endif
