#include "Quaternion.h"

// Quaternion::Quaternion()
//{

//}

Quaternion operator*(float s, const Quaternion &q) { return q * s; }

Quaternion::Quaternion() : Quaternion(1.f, 0.f, 0.f, 0.f) {}

Quaternion::Quaternion(float w, float x, float y, float z) : w(w), x(x), y(y), z(z) {}

Quaternion Quaternion::AxisAngle(const Vector3 &axis, float angle) {
    if (!std::isfinite(angle) || !axis.isValid() || axis.norm2() == 0)
        throw std::invalid_argument("Invalid rotation axis/angle");
    auto n = axis.normalized();
    float h = angle * .5f, s = std::sin(h);
    return Quaternion(std::cos(h), n.x() * s, n.y() * s, n.z() * s);
}

Vector3 Quaternion::rotate(const Vector3 &v) const {
    double n2 = double(w) * w + double(x) * x + double(y) * y + double(z) * z;
    if (n2 == 0 || !std::isfinite(n2))
        throw std::domain_error("Invalid rotation quaternion");
    if (!v.isValid())
        return Vector3::invalid;
    if (n2 < 1. / 16 || n2 > 16 ||
        std::max({std::abs(v.x()), std::abs(v.y()), std::abs(v.z())}) >
            std::numeric_limits<float>::max() / 128) {
        const double tx = 2 * (double(y) * v.z() - double(z) * v.y());
        const double ty = 2 * (double(z) * v.x() - double(x) * v.z());
        const double tz = 2 * (double(x) * v.y() - double(y) * v.x());
        auto represent = [](double value) {
            const double limit = std::numeric_limits<float>::max();
            if (!std::isfinite(value) ||
                std::abs(value) > limit * (1 + 32 * std::numeric_limits<double>::epsilon()))
                throw std::overflow_error("Rotated vector is outside float range");
            return float(std::clamp(value, -limit, limit));
        };
        return Vector3(represent(v.x() + (w * tx + y * tz - z * ty) / n2),
                       represent(v.y() + (w * ty + z * tx - x * tz) / n2),
                       represent(v.z() + (w * tz + x * ty - y * tx) / n2));
    }
    const float tx = 2 * (y * v.z() - z * v.y()), ty = 2 * (z * v.x() - x * v.z()),
                tz = 2 * (x * v.y() - y * v.x());
    const float scale = std::abs(n2 - 1) < 1e-6 ? 1.f : float(1 / n2);
    return Vector3(v.x() + scale * (w * tx + y * tz - z * ty), v.y() + scale * (w * ty + z * tx - x * tz),
                   v.z() + scale * (w * tz + x * ty - y * tx), v.valid);
}

Quaternion Quaternion::slerp(const Quaternion &q1, const Quaternion &q2, float t) {
    if (!std::isfinite(t))
        throw std::invalid_argument("Nonfinite interpolation parameter");
    auto a = q1.normalized(), b = q2.normalized();
    float d = a.dot(b);
    if (d < 0) {
        b = -b;
        d = -d;
    }
    d = std::clamp(d, -1.f, 1.f);
    if (d > 1 - 1e-6f)
        return (a + (b - a) * t).normalized();
    float angle = std::acos(d), denom = std::sin(angle);
    return (a * (std::sin((1 - t) * angle) / denom) + b * (std::sin(t * angle) / denom)).normalized();
}

float Quaternion::dot(const Quaternion &b) const { return w * b.w + x * b.x + y * b.y + z * b.z; }

Quaternion Quaternion::operator*(float s) const { return Quaternion(w * s, x * s, y * s, z * s); }

Quaternion Quaternion::operator/(float s) const {
    if (s == 0 || !std::isfinite(s))
        throw std::domain_error("Invalid quaternion divisor");
    return *this * (1.f / s);
}

Quaternion Quaternion::operator+(const Quaternion &b) const {
    return Quaternion(w + b.w, x + b.x, y + b.y, z + b.z);
}

Quaternion Quaternion::operator-(const Quaternion &b) const {
    return Quaternion(w - b.w, x - b.x, y - b.y, z - b.z);
}

Quaternion Quaternion::operator-() const { return *this * -1.f; }

Vector3 Quaternion::toVector3() const { return Vector3(x, y, z); }

float Quaternion::norm() const {
    return float(std::hypot(std::hypot(double(w), double(x)), std::hypot(double(y), double(z))));
}

Quaternion Quaternion::normalized() const {
    double n = std::hypot(std::hypot(double(w), double(x)), std::hypot(double(y), double(z)));
    if (n == 0 || !std::isfinite(n))
        throw std::domain_error("Cannot normalize zero/nonfinite quaternion");
    return Quaternion(w / n, x / n, y / n, z / n);
}

Quaternion Quaternion::conjugate() const { return Quaternion(w, -x, -y, -z); }

Quaternion Quaternion::inverse() const {
    double n = double(w) * w + double(x) * x + double(y) * y + double(z) * z;
    if (n == 0 || !std::isfinite(n))
        throw std::domain_error("Cannot invert zero/nonfinite quaternion");
    for (double value : {w / n, -x / n, -y / n, -z / n})
        if (std::abs(value) > std::numeric_limits<float>::max())
            throw std::overflow_error("Quaternion inverse is outside float range");
    return Quaternion(w / n, -x / n, -y / n, -z / n);
}

Quaternion Quaternion::lerp(const Quaternion &q1, const Quaternion &q2, float t) {
    return (1 - t) * q1 + t * q2;
}

Quaternion Quaternion::squad(const Quaternion &q1, const Quaternion &q2, const Quaternion &a,
                             const Quaternion &b, float t) {
    Quaternion c = slerp(q1, q2, t);
    Quaternion d = slerp(a, b, t);
    return slerp(c, d, 2.0f * t * (1.0f - t));
}

Quaternion Quaternion::log() const {
    const double length = std::hypot(std::hypot(double(w), double(x)), std::hypot(double(y), double(z)));
    if (length == 0 || !std::isfinite(length))
        throw std::domain_error("Invalid quaternion logarithm");
    const double v = std::hypot(double(x), double(y), double(z));
    if (v == 0)
        return Quaternion(std::log(length), w < 0 ? float(M_PI) : 0.f, 0, 0);
    const double scale = std::atan2(v, double(w)) / v;
    return Quaternion(std::log(length), x * scale, y * scale, z * scale);
}

Quaternion Quaternion::exp() const {
    if (!std::isfinite(w) || !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
        throw std::domain_error("Invalid quaternion exponential");
    const double length = std::hypot(double(x), double(y), double(z)), magnitude = std::exp(double(w));
    if (!std::isfinite(magnitude) || magnitude > std::numeric_limits<float>::max())
        throw std::overflow_error("Quaternion exponential overflow");
    const double scale = length == 0 ? magnitude : magnitude * std::sin(length) / length;
    return Quaternion(magnitude * std::cos(length), x * scale, y * scale, z * scale);
}

Vector3 Quaternion::toAxisAngle(float &angle) const {
    auto q = normalized();
    float c = std::clamp(q.w, -1.f, 1.f);
    angle = 2 * std::acos(c);
    float s = std::sqrt(std::max(0.f, 1 - c * c));
    return s < 1e-6f ? Vector3(1, 0, 0) : Vector3(q.x / s, q.y / s, q.z / s);
}

Quaternion Quaternion::fromAxisAngle(const Vector3 &axis, float angle) { return AxisAngle(axis, angle); }

Matrix Quaternion::toRotationMatrix() const {
    auto q = normalized();
    const float w = q.w, x = q.x, y = q.y, z = q.z;
    return Matrix(std::vector<std::vector<float>>{
        {1 - 2 * y * y - 2 * z * z, 2 * x * y - 2 * z * w, 2 * x * z + 2 * y * w},
        {2 * x * y + 2 * z * w, 1 - 2 * x * x - 2 * z * z, 2 * y * z - 2 * x * w},
        {2 * x * z - 2 * y * w, 2 * y * z + 2 * x * w, 1 - 2 * x * x - 2 * y * y}});
}

Quaternion Quaternion::fromRotationMatrix(const Matrix &m) {
    if (m.rows() != 3 || m.cols() != 3 || m.size() != 9)
        throw std::invalid_argument("Rotation matrix must be 3x3");
    for (float x : m.data())
        if (!std::isfinite(x))
            throw std::invalid_argument("Nonfinite rotation matrix");
    double t = m.trace();
    Quaternion q;
    if (t > 0) {
        double s = 2 * std::sqrt(t + 1);
        q = Quaternion(s / 4, (m(1, 2) - m(2, 1)) / s, (m(2, 0) - m(0, 2)) / s, (m(0, 1) - m(1, 0)) / s);
    } else {
        int i = 0;
        for (int k = 1; k < 3; ++k)
            if (m(k, k) > m(i, i))
                i = k;
        int j = (i + 1) % 3, k = (i + 2) % 3;
        double s = 2 * std::sqrt(std::max(0., 1. + m(i, i) - m(j, j) - m(k, k)));
        if (s == 0)
            throw std::invalid_argument("Invalid rotation matrix");
        float v[3];
        v[i] = s / 4;
        v[j] = (m(j, i) + m(i, j)) / s;
        v[k] = (m(k, i) + m(i, k)) / s;
        q = Quaternion((m(j, k) - m(k, j)) / s, v[0], v[1], v[2]);
    }
    return q.normalized();
}

Quaternion Quaternion::operator*(const Quaternion &b) const {
    return Quaternion(w * b.w - x * b.x - y * b.y - z * b.z, // new w
                      w * b.x + x * b.w + y * b.z - z * b.y, // new x
                      w * b.y - x * b.z + y * b.w + z * b.x, // new y
                      w * b.z + x * b.y - y * b.x + z * b.w  // new z
    );
}
