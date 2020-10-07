#include "core.h"

const int SimplexNoise::perm[512] = {
  151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142,
  8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117,
  35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71,
  134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41,
  55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89,
  18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226,
  250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182,
  189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43,
  172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97,
  228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239,
  107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
  138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180,

  151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142,
  8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117,
  35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71,
  134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41,
  55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89,
  18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226,
  250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182,
  189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43,
  172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97,
  228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239,
  107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
  138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
};
const int SimplexNoise::grad2[8][2] = {
  { 1, 1 }, { -1, 1 }, { 1, -1 }, { -1, -1 },
  { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 }
};
const double SimplexNoise::F2 = 0.5 * (sqrt(3.0) - 1.0);
const double SimplexNoise::G2 = (3.0 - sqrt(3.0)) / 6.0;

double SimplexNoise::at(double x, double y) const
{
	// Noise contributions from the three corners
	double n[3] = { 0.0, 0.0, 0.0 };

	// Skew the input space to determine which simplex cell we are in

	// Hairy factor for 2D
	double s = (x + y) * F2;

	int i = Math::Integer(x + s);
	int j = Math::Integer(y + s);

	double t = (i + j) * G2;

	// Unskew the cell origin back to (x,y) space
	double X0 = i - t;
	double Y0 = j - t;

	// The x,y distances from the cell origin
	double x0 = x - X0;
	double y0 = y - Y0;

	// For the 2D case, the simplex shape is an equilateral triangle.
	// Determine which simplex we are in.

	int i1, j1; // Offsets for second (middle) corner of simplex in (i,j) coords

	if (x0 > y0) { i1 = 1; j1 = 0; } // lower triangle, XY order: (0,0)->(1,0)->(1,1)
	else { i1 = 0; j1 = 1; } // upper triangle, YX order: (0,0)->(0,1)->(1,1)

	// A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
	// a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where c = (3-sqrt(3))/6

	double x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
	double y1 = y0 - j1 + G2;
	double x2 = x0 - 1.0 + 2.0 * G2; // Offsets for last corner in (x,y) unskewed coords
	double y2 = y0 - 1.0 + 2.0 * G2;

	// Work out the hashed gradient indices of the three simplex corners
	int ii = i & 255;
	int jj = j & 255;

	int gi0 = perm[ii + perm[jj]] % 8;
	int gi1 = perm[ii + i1 + perm[jj + j1]] % 8;
	int gi2 = perm[ii + 1 + perm[jj + 1]] % 8;

	// Calculate the contribution from the three corners

	double t0 = 0.5 - x0 * x0 - y0 * y0;
	if (t0 >= 0)
	{
		t0 *= t0;
		n[0] = t0 * t0 * dot(grad2[gi0], x0, y0);
	}

	double t1 = 0.5 - x1 * x1 - y1 * y1;
	if (t1 >= 0)
	{
		t1 *= t1;
		n[1] = t1 * t1 * dot(grad2[gi1], x1, y1);
	}

	double t2 = 0.5 - x2 * x2 - y2 * y2;
	if (t2 >= 0)
	{
		t2 *= t2;
		n[2] = t2 * t2 * dot(grad2[gi2], x2, y2);
	}

	// Add contributions from each corner to get the final noise value.
	// The result is scaled to return values in the interval [-1,1].
	return 70.0 * (n[0] + n[1] + n[2]);
}


bool Box2::Intersect(const Vector2 &s0, const Vector2 &s1, double &tmin, double &tmax)
{
    const double epsilon = 1.0e-5;

    tmin = -1e16;
    tmax = 1e16;

    const Vector2& a = bmin;
    const Vector2& b = bmax;
    Vector2 p = s0;
    Vector2 d = s1 - s0;

    double t;
    // Ox
    if (d[0] < -epsilon) {
        t = (a[0] - p[0]) / d[0];
        if (t < tmin)
            return false;
        if (t <= tmax)
            tmax = t;
        t = (b[0] - p[0]) / d[0];
        if (t >= tmin) {
            if (t > tmax)
                return false;
            tmin = t;
        }
    }
    else if (d[0] > epsilon) {
        t = (b[0] - p[0]) / d[0];
        if (t < tmin)
            return false;
        if (t <= tmax)
            tmax = t;
        t = (a[0] - p[0]) / d[0];
        if (t >= tmin) {
            if (t > tmax)
                return false;
            tmin = t;
        }
    }
    else if (p[0]<a[0] || p[0]>b[0])
        return false;

    // Oy
    if (d[1] < -epsilon) {
        t = (a[1] - p[1]) / d[1];
        if (t < tmin)
            return false;
        if (t <= tmax)
            tmax = t;
        t = (b[1] - p[1]) / d[1];
        if (t >= tmin) {
            if (t > tmax)
                return false;
            tmin = t;
        }
    }
    else if (d[1] > epsilon) {
        t = (b[1] - p[1]) / d[1];
        if (t < tmin)
            return false;
        if (t <= tmax)
            tmax = t;
        t = (a[1] - p[1]) / d[1];
        if (t >= tmin) {
            if (t > tmax)
                return false;
            tmin = t;
        }
    }
    else if (p[1]<a[1] || p[1]>b[1])
        return false;

    return true;
}

Camera::Camera()
{
    Camera::eye = Vector3(0.0);
    Camera::at = Vector3(0.0, 1.0, 0.0);
    Camera::up = Vector3(0.0, 0.0, 1.0);

    // Near and far planes
    Camera::nearplane = 1.0;
    Camera::farplane = 1000.0;

    // Aperture
    Camera::cah = 0.980;
    Camera::cav = 0.735;
    Camera::fl = 35.0;
}

Camera::Camera(const Vector3& eye, const Vector3& at, const Vector3& up, double near, double far)
{
    Camera::eye = eye;
    Camera::at = at;
    Camera::up = up;

    // Near and far planes
    Camera::nearplane = near;
    Camera::farplane = far;

    // Aperture
    Camera::cah = 0.980;
    Camera::cav = 0.735;
    Camera::fl = 35.0;
}

double Camera::getAngleOfViewH(double, double) const
{
    return 2.0 * atan(cah * 25.4 * 0.5 / fl);
}

double Camera::getAngleOfViewV(double w, double h) const
{
    double avh = getAngleOfViewH(w, h);
    return 2.0 * atan(tan(avh / 2.0) * double(h) / double(w));
}

void Camera::upDownRound(double a)
{
    Vector3 z = at - eye;
    double length = Norm(z);
    z = z/length;
    Vector3 left = Normalized(cross(up, z));

    // Rotate
    z = z * cos(a) + up * sin(a);

    // Update Vector
    up = cross(z, left);
    eye = at - z * length;
}

void Camera::leftRightRound(double a)
{
    Vector3 e = eye - at;
    Vector3 left = cross(up, e);
    e = Vector3(e[0] * cos(a) - e[1] * sin(a), e[0] * sin(a) + e[1] * cos(a), e[2]);
    left = Vector3(left[0] * cos(a) - left[1] * sin(a), left[0] * sin(a) + left[1] * cos(a), 0.0);
    up = Normalized(cross(left, -e));
    eye = at + e;
}

void Camera::backForth(double a, bool moveAt)
{
    Vector3 z = at - eye;
    double length = Norm(z);
    z = z/length;
    eye = eye + a * z;
    if (moveAt) {
        at = at + a * z;
    }
}

void Camera::upDownPlane(double a)
{
    Vector3 z = at - eye;
    double length = Norm(z);
    z = z/length;
    Vector3 left = Normalized(cross(Vector3(0, 0, 1), z));

    eye = eye + a * cross(z, left);
    at = at + a * cross(z, left);
}

void Camera::leftRightPlane(double a)
{
    Vector3 z = at - eye;
    z[2] = 0.0;
    double length = Norm(z);
    z = z/length;
    Vector3 left = Normalized(cross(Vector3(0, 0, 1), z));

    eye = eye + a * left;
    at = at + a * left;
}

Ray Camera::pixelToRay(int px, int py, int w, int h) const
{
    // Get coordinates
    Vector3 view = getViewDir();
    Vector3 horizontal = Normalized(cross(view, up));
    Vector3 vertical = Normalized(cross(horizontal, view));

    double length = 1.0;

    // Convert to radians
    double rad = getAngleOfViewV(w, h);  // fov

    double vLength = tan(rad / 2.0) * length;
    double hLength = vLength * (double(w) / double(h));
    vertical = vertical*vLength;
    horizontal = horizontal*hLength;

    // Translate mouse coordinates so that the origin lies in the center of the view port
    double x = px - w / 2.0;
    double y = h / 2.0 - py;

    // Scale mouse coordinates so that half the view port width and height becomes 1.0
    x /= w / 2.0;
    y /= h / 2.0;

    // Direction is a linear combination to compute intersection of picking ray with view port plane
    return Ray(eye, Normalized(view * length + horizontal * x + vertical * y));
}

Camera Camera::View(const Box3& box)
{
    Vector3 v = 0.5*(box.getMax() - box.getMin());
	v[2] = 0;
	double r = Norm(v);
	v = 2*v;
	v[2] = -r;
    return Camera(box.center() - v, box.center(), Vector3(0.0, 0.0, 1.0), r, 3*r);
}


Vector3 ColorPalette::getColor(double t) const {
    if (colors.size() == 0) return Vector3(1);
    if (colors.size() == 1) return colors[0];

    if (t < anchors.front()) return colors.front();
    if (t > anchors.back()) return colors.back();
    for (int i = 0; i < int(colors.size() - 1); i++) {
        if (t < anchors[i+1]) {
            double s = Math::LinearStep(t, anchors[i], anchors[i+1]);
            return (1-s)*colors[i] + s*colors[i+1];
        }
    }

    return Vector3(1);
}
