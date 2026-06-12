#include "BSpline.h"

#include "Curves.h"

BSpline::BSpline()
{}
BSpline::BSpline(const std::vector<Vector3> &points, bool clamped)
    : points(points)
{
    size_t knotSize = points.size() + degree + 1;
    this->knots = std::vector<float>(knotSize);

    if (clamped) {
        for (size_t i = 0; i < knotSize; i++) {
            if (int(i) < degree + 1) {
                knots[i] = 0;
            } else if (i > knotSize - (degree + 2)) {
                knots[i] = 1;
            } else {
                float t = float(i - (degree)) / float((knotSize - 1) - (degree) * 2);
                knots[i] = t;
            }
        }
    } else {
        for (size_t i = 0; i < knotSize; i++) {
            float t = float(i) / float(knotSize - 1);
            knots[i] = t;
        }
    }
}

BSpline::BSpline(const std::vector<Vector3>& points, const std::vector<float>& knots)
    : points(points), knots(knots), degree(knots.size() - (points.size() + 1))
{
}

BSpline::BSpline(const Curve &curve)
    : BSpline(curve.getPath(), true)
{}

std::vector<Vector3> BSpline::getPath(int numberOfPoints) const
{
    if (numberOfPoints < 0)
        return this->points;

    std::vector<Vector3> points;
    for (int i = 0; i < numberOfPoints; i++) {
        float t = float(i) / float(numberOfPoints - 1);
        points.push_back(getPoint(t));
        if (!points.back().isValid()) {
            std::cout << "Error on t=" << t << std::endl;
        }
    }
    return points;
}

/*Vector3 BSpline::getPoint(float _x) const
{
    float x = std::clamp(_x, 0.f, 1.f);
    size_t nbPoints = numPoints();

    if (nbPoints == 0) return Vector3::invalid;
    else if (nbPoints == 1) return points[0];
    else if (nbPoints == 2) return points[0] * (1.f - x) + points[1] * x;

    float res = 1.f / (float)(numSegments());
    int iFloor = int(x / res);
    int iCeil = iFloor + 1;
    float resFloor = iFloor * res;
    float resCeil = iCeil * res;
    float x_prime = map(x, resFloor, resCeil, 0.f, 1.f);
    iFloor += 1;
    iCeil += 1;
    if (_x >= 1.f) {
        x_prime = 1.f;
        iFloor = numSegments();
        iCeil = iFloor+1;
    }

    auto P0 = (iFloor > 0 || closed ? points[pointIndex(iFloor - 1)] : points[pointIndex(iFloor)]);
    auto P1 = points[pointIndex(iFloor)];
    auto P2 = points[pointIndex(iCeil)];
    auto P3 = (iCeil + 1 < nbPoints || closed ? points[pointIndex(iCeil + 1)] : points[pointIndex(iCeil)]);

    return ((-P0 + P1 * 3.f - P2 * 3.f + P3) * x_prime * x_prime * x_prime
            + (P0 * 3.f - P1 * 6.f + P2 * 3.f) * x_prime * x_prime
            + (P0 * -3.f + P2 * 3.f) * x_prime
            + (P0 + P1 * 4.f + P2)) / 6.f;
}*/

Vector3 BSpline::getPoint(float t) const
{
    const int p = degree;
    const int n = int(points.size()) - 1;

    if (points.empty()) return Vector3::invalid;
    if (points.size() == 1) return points[0];

    const auto& U = knots;

    // Valid knot vector size:
    // knots.size() == points.size() + degree + 1

    const float uMin = U[p];
    const float uMax = U[n + 1];

    float u = std::clamp(t, uMin, uMax);

    // Special case only because the span convention is U[k] <= u < U[k + 1]
    if (u == uMax) {
        u = std::nextafter(uMax, uMin);
    }

    int k = p;

    for (int i = p; i <= n; ++i) {
        if (U[i] <= u && u < U[i + 1]) {
            k = i;
            break;
        }
    }

    std::vector<Vector3> d(p + 1);

    for (int j = 0; j <= p; ++j) {
        d[j] = points[k - p + j];
    }

    for (int r = 1; r <= p; ++r) {
        for (int j = p; j >= r; --j) {
            const int i = k - p + j;

            const float denom = U[i + p + 1 - r] - U[i];

            float alpha = 0.f;
            if (denom != 0.f) {
                alpha = (u - U[i]) / denom;
            }

            d[j] = d[j - 1] * (1.f - alpha) + d[j] * alpha;
        }
    }

    return d[p];
}


Vector3 BSpline::getDerivative(float x, bool normalize) const
{
    return Vector3::invalid;
}

Vector3 BSpline::getSecondDerivative(float x, bool normalize) const
{
    return Vector3::invalid;
}

float BSpline::estimateClosestTime(const Vector3 &pos) const
{
    return toBezier(*this).estimateClosestTime(pos);
}

Vector3 BSpline::estimateClosestPos(const Vector3 &pos) const
{
    return getPoint(estimateClosestTime(pos));
}

float BSpline::estimateSqrDistanceFrom(const Vector3 &pos) const
{
    return (estimateClosestPos(pos) - pos).norm2();
}

float BSpline::length() const
{
    float length = 0;
    const auto& points = getPath(numPoints());
    for (size_t i = 0; i < points.size() - 1; i++) {
        length += (points[i] - points[i+1]).norm();
    }
    return length;
}

size_t BSpline::getIndex(int i) {
    return (i + numPoints()) % numPoints();
}

BSpline& BSpline::setPoint(int i, const Vector3 &newPos)
{
    this->points[pointIndex(i)] = newPos;
    return *this;
}

std::vector<Vector3> BSpline::getPoints() const {
    return this->points;
}

Vector3 &BSpline::get(int i) {
    return this->points[pointIndex(i)];
}

Vector3 BSpline::get(int i) const {
    return this->points[pointIndex(i)];
}

size_t BSpline::numSegments() const
{
    return numPoints() - (closed ? 0 : 3);
}

std::vector<Vector3> BSpline::solveLinearSystem(std::vector<std::vector<float>> A, std::vector<Vector3> b)
{
    const int n = int(b.size());

    for (int i = 0; i < n; ++i)
    {
        int pivot = i;
        for (int r = i + 1; r < n; ++r)
        {
            if (std::abs(A[r][i]) > std::abs(A[pivot][i]))
                pivot = r;
        }

        std::swap(A[i], A[pivot]);
        std::swap(b[i], b[pivot]);

        float div = A[i][i];
        if (std::abs(div) < 1e-8f)
            continue;

        for (int c = i; c < n; ++c)
            A[i][c] /= div;

        b[i] /= div;

        for (int r = 0; r < n; ++r)
        {
            if (r == i)
                continue;

            float factor = A[r][i];

            for (int c = i; c < n; ++c)
                A[r][c] -= factor * A[i][c];

            b[r] -= factor * b[i];
        }
    }

    return b;
}

float BSpline::basis(int i, int p, float u, const std::vector<float>& U)
{
    if (p == 0)
    {
        if ((U[i] <= u && u < U[i + 1]) ||
            (u == 1.f && U[i + 1] == 1.f))
            return 1.f;

        return 0.f;
    }

    float left = 0.f;
    float right = 0.f;

    float denom1 = U[i + p] - U[i];
    if (denom1 != 0.f)
        left = ((u - U[i]) / denom1) * basis(i, p - 1, u, U);

    float denom2 = U[i + p + 1] - U[i + 1];
    if (denom2 != 0.f)
        right = ((U[i + p + 1] - u) / denom2) * basis(i + 1, p - 1, u, U);

    return left + right;
}

std::vector<float> BSpline::makeUniformClampedKnots(int numCtrlPoints, int degree)
{
    const int p = degree;
    const int n = numCtrlPoints - 1;
    const int m = n + p + 1;

    std::vector<float> U(m + 1, 0.f);

    for (int i = 0; i <= p; ++i)
        U[i] = 0.f;

    for (int i = m - p; i <= m; ++i)
        U[i] = 1.f;

    const int numInterior = m - 2 * p;

    for (int j = 1; j < numInterior; ++j)
        U[p + j] = float(j) / float(numInterior);

    return U;
}

BSpline& BSpline::resamplePoints(int newNbPoints)
{
    const int p = degree;

    if (newNbPoints < 0 || newNbPoints == numPoints())
        newNbPoints = numPoints();

    if (newNbPoints < p + 1)
        newNbPoints = p + 1;

    const int n = newNbPoints - 1;

    BSpline old = *this;

    std::vector<float> newKnots = makeUniformClampedKnots(newNbPoints, p);

    std::vector<float> params(newNbPoints);
    for (int i = 0; i < newNbPoints; ++i)
        params[i] = float(i) / float(newNbPoints - 1);

    std::vector<Vector3> samples(newNbPoints);
    for (int i = 0; i < newNbPoints; ++i)
        samples[i] = old.getPoint(params[i]);

    std::vector<std::vector<float>> A(newNbPoints, std::vector<float>(newNbPoints, 0.f));

    for (int row = 0; row < newNbPoints; ++row)
    {
        float u = params[row];

        for (int col = 0; col < newNbPoints; ++col)
            A[row][col] = basis(col, p, u, newKnots);
    }

    std::vector<Vector3> newPoints = solveLinearSystem(A, samples);

    knots = std::move(newKnots);
    points = std::move(newPoints);

    return *this;
}

BSpline& BSpline::reverseVertices()
{
    std::reverse(this->points.begin(), this->points.end());
    return *this;
}

std::pair<Vector3, Vector3> BSpline::AABBox() const
{
    Vector3 mini = Vector3::max();
    Vector3 maxi = Vector3::min();

    for (auto& p : points) {
        mini.x() = std::min(mini.x(), p.x());
        mini.y() = std::min(mini.y(), p.y());
        mini.z() = std::min(mini.z(), p.z());
        maxi.x() = std::max(maxi.x(), p.x());
        maxi.y() = std::max(maxi.y(), p.y());
        maxi.z() = std::max(maxi.z(), p.z());
    }
    return {mini, maxi};
}

BSpline& BSpline::scale(const Vector3 &factor)
{
    for (auto& p : points) {
        p *= factor;
    }
    return *this;
}

BSpline& BSpline::translate(const Vector3 &translation)
{
    for (auto& p : points) {
        p += translation;
    }
    return *this;
}

BSpline& BSpline::removeDuplicates()
{
    return *this;
}

size_t BSpline::size() const {
    return numPoints();
}

size_t BSpline::numPoints() const
{
    return points.size();
}

std::vector<Vector3>::const_iterator BSpline::begin() const {
    return points.begin();
}

std::vector<Vector3>::const_iterator BSpline::end() const {
    return points.end();
}

std::vector<Vector3>::iterator BSpline::begin() {
    return points.begin();
}

std::vector<Vector3>::iterator BSpline::end() {
    return points.end();
}

std::string BSpline::toString() const
{
    std::ostringstream out;
    out << "B-Spline with " << this->points.size() << " points (" << (closed ? "closed" : "not closed") << ") :\n";
    for (auto& p : this->points)
        out << "- " << p << "\n";
    return out.str();
}

BSpline& BSpline::close()
{
    if (closed) return *this;
    Curve::close();
    // if (this->points.size() > 0 && this->points.front() != this->points.back()) {
    const auto P0 = get(-2);
    const auto P1 = get(-1);
    const auto P2 = get(0);
    const auto P3 = get(1);
    // points.insert(points.begin(), P1);
    // points.insert(points.begin(), P0);
    // points.push_back(P2);
    // points.push_back(P3);
    // this->points.insert(points.begin(), points.back());
    // this->points.push_back(points[0]);
    // this->points.push_back(points[1]);
        // this->points.push_back(points[2]);
    // }
    return *this;
}

BSpline &BSpline::reset() {
    points.clear();
    return *this;
}

std::vector<std::shared_ptr<Curve> > BSpline::slice(float t) const
{
    if (t <= 0.f || t >= 1.f) return { std::make_shared<BSpline>(*this) };
    auto p1 = *this;
    // Adding a clamped point
    while (p1.knotMultiplicity(t) < p1.degree + 1) {
        p1.insertKnot(t);
    }
    // Position of the knot in the array
    size_t k = std::distance(p1.knots.begin(), std::find(p1.knots.begin(), p1.knots.end(), t));
    auto p2 = p1;

    // Removing right-side and left-side curves from p1 and p2
    p1.knots.erase(p1.knots.begin() + k + p1.degree + 1, p1.knots.end());
    p1.points.erase(p1.points.begin() + k, p1.points.end());
    p2.knots.erase(p2.knots.begin(), p2.knots.begin() + k);
    p2.points.erase(p2.points.begin(), p2.points.begin() + k);

    // Remapping knots between  0 and 1
    float min1 = p1.knots.front(), min2 = p2.knots.front();
    float max1 = p1.knots.back(), max2 = p2.knots.back();

    for (auto& k : p1.knots)
        k = map(k, min1, max1, 0.f, 1.f);
    for (auto& k : p2.knots)
        k = map(k, min2, max2, 0.f, 1.f);

    return {std::make_shared<BSpline>(p1), std::make_shared<BSpline>(p2)};
}
    /*
    size_t lowerIndex = this->timeToLowerIndex(t);
    std::vector<Vector3> firstPoints(lowerIndex + 2);
    std::vector<Vector3> lastPoints(numPoints() - lowerIndex);

    for (size_t i = 0; i < numSegments() + 1; i++) {
        if (i <= lowerIndex) {
            firstPoints[i] = points[i];
        } else {
            lastPoints[i - lowerIndex] = points[pointIndex(i)];
        }
    }
    // float u = ::map(t, indexToTime(lowerIndex), indexToTime(lowerIndex + 1), 0.f, 1.f);
    // const auto p = points[lowerIndex] * (1.f - u) + points[pointIndex(lowerIndex + 1)] * u; //getPoint(t);
    firstPoints.back() = points[pointIndex(lowerIndex + 1)]; // p;
    lastPoints.front() = points[lowerIndex]; //p;
    BSpline p1(firstPoints);
    BSpline p2(lastPoints);

    // p1.addPoint(p2.get(1));
    // p2.insertPoint(0, p1.get(-2));
    return {std::make_shared<BSpline>(p1), std::make_shared<BSpline>(p2)};

}*/

BSpline& BSpline::setDegree(int newDegree)
{
    this->degree = newDegree;
    return *this;
}

std::vector<float> BSpline::uniqueInternalKnots() const
{
    std::vector<float> res(this->knots.size() - 2 * (degree + 1));
    for (size_t i = 0; i < res.size(); i++) {
        res[i] = knots[degree + 1 + i];
    }
    return res;
}

std::vector<float> BSpline::uniqueInternalKnots(float minKnot, float maxKnot) const
{
    std::vector<float> res;
    for (size_t i = 0; i < knots.size(); i++) {
        if (minKnot < knots[i] && knots[i] < maxKnot)
            res.push_back(knots[i]);
    }
    return res;
}

int BSpline::knotMultiplicity(float u) const
{
    return std::count(knots.begin(), knots.end(), u);
}

size_t BSpline::insertKnot(float u)
{
    const int p = degree;
    const int n = int(points.size()) - 1;

    auto oldPoints = points;
    auto oldKnots = knots;

    // Find span k such that oldKnots[k] <= u < oldKnots[k + 1]
    int k = p;

    if (u >= oldKnots[n + 1]) {
        k = n;
    } else {
        for (int i = p; i <= n; ++i) {
            if (oldKnots[i] <= u && u < oldKnots[i + 1]) {
                k = i;
                break;
            }
        }
    }

    // Insert knot
    knots.insert(knots.begin() + k + 1, u);

    // New control point count = old + 1
    std::vector<Vector3> newPoints(oldPoints.size() + 1);

    // Unaffected points before insertion region
    for (int i = 0; i <= k - p; ++i) {
        newPoints[i] = oldPoints[i];
    }

    // Unaffected points after insertion region
    for (int i = k; i <= n; ++i) {
        newPoints[i + 1] = oldPoints[i];
    }

    // Recomputed local points
    for (int i = k - p + 1; i <= k; ++i) {
        float denom = oldKnots[i + p] - oldKnots[i];

        float alpha = 0.f;
        if (denom != 0.f) {
            alpha = (u - oldKnots[i]) / denom;
        }

        newPoints[i] =
            oldPoints[i - 1] * (1.f - alpha)
            + oldPoints[i] * alpha;
    }

    points = std::move(newPoints);

    return k;
}

BSpline& BSpline::removeKnotAtKnotIndex(int knotIndex, float tolerance)
{
    const int p = degree;
    const int n = static_cast<int>(points.size()) - 1;

    if (p < 1 || n < p)
        return *this;

    if (knotIndex < 0 || knotIndex >= static_cast<int>(knots.size()))
        return *this;

    const float u = knots[knotIndex];

    // Find multiplicity and last occurrence of u.
    int r = knotIndex;
    while (r + 1 < static_cast<int>(knots.size()) &&
           std::abs(knots[r + 1] - u) < 1e-6f)
    {
        ++r;
    }

    int firstOccurrence = r;
    while (firstOccurrence - 1 >= 0 &&
           std::abs(knots[firstOccurrence - 1] - u) < 1e-6f)
    {
        --firstOccurrence;
    }

    const int s = r - firstOccurrence + 1; // multiplicity

    const int first = r - p;
    const int last  = r - s;

    if (first < 1 || last + 1 > n)
        return *this;

    std::vector<Vector3> temp(last - first + 2);

    temp[0] = points[first - 1];
    temp[last + 1 - first] = points[last + 1];

    int i = first;
    int j = last;
    int ii = 1;
    int jj = last - first;

    bool removable = true;

    while (j - i > 0)
    {
        const float alphaI =
            (u - knots[i]) /
            (knots[i + p + 1] - knots[i]);

        const float alphaJ =
            (u - knots[j]) /
            (knots[j + p + 1] - knots[j]);

        if (std::abs(alphaI) < 1e-6f ||
            std::abs(1.f - alphaJ) < 1e-6f)
        {
            removable = false;
            break;
        }

        temp[ii] =
            (points[i] - (1.f - alphaI) * temp[ii - 1]) / alphaI;

        temp[jj] =
            (points[j] - alphaJ * temp[jj + 1]) / (1.f - alphaJ);

        ++i;
        --j;
        ++ii;
        --jj;
    }

    if (!removable)
        return *this;

    // Final consistency check.
    if (j - i == 0)
    {
        const float alpha =
            (u - knots[i]) /
            (knots[i + p + 1] - knots[i]);

        const Vector3 reconstructed =
            (1.f - alpha) * temp[ii - 1] + alpha * temp[jj + 1];

        if ((reconstructed - points[i]).length() > tolerance)
            return *this;
    }
    else
    {
        if ((temp[ii - 1] - temp[jj + 1]).length() > tolerance)
            return *this;
    }

    // Commit: replace affected control-point block with the recovered one.
    std::vector<Vector3> newPoints;

    for (int k = 0; k < first - 1; ++k)
        newPoints.push_back(points[k]);

    for (const Vector3& q : temp)
        newPoints.push_back(q);

    for (int k = last + 2; k <= n; ++k)
        newPoints.push_back(points[k]);

    points = std::move(newPoints);

    knots.erase(knots.begin() + r);

    return *this;
}

BSpline& BSpline::removeKnot(int pointIndex, float tolerance)
{
    const int p = degree;

    if (p < 1)
        return *this;

    int knotIndex = pointIndex + p - 1;

    if (knotIndex < 0 || knotIndex >= static_cast<int>(knots.size()))
        return *this;

    return removeKnotAtKnotIndex(knotIndex, tolerance);
}

Vector3& BSpline::operator[](size_t i)
{
    return this->points[pointIndex(i)];
}

void BSpline::addPoint(const Vector3 &newPoint) {
    this->points.push_back(newPoint);
}

BSpline &BSpline::insertPoint(int i, const Vector3 &newPos) {
    this->points.insert(points.begin() + pointIndex(i), newPos); return *this;
}

BSpline &BSpline::removePoint(int i) {
    this->points.erase(points.begin() + pointIndex(i)); return *this;
}
const Vector3& BSpline::operator[](size_t i) const
{
    return this->points[(i + size()) % size()];
}

