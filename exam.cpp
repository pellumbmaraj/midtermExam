#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>

using namespace std;

const int p = 2029; // Field modulus
const int a = 6;
const int b = 5;

struct Point {
    long long x, y;
    bool at_infinity;

    Point() : x(0), y(0), at_infinity(true) {}
    Point(long long x_, long long y_) : x(x_), y(y_), at_infinity(false) {}

    bool operator==(const Point& other) const {
        if (at_infinity && other.at_infinity) return true;
        if (at_infinity || other.at_infinity) return false;
        return x == other.x && y == other.y;
    }
};

// Modular inverse using Extended Euclidean Algorithm
long long modinv(long long a, long long mod) {
    long long t = 0, newt = 1;
    long long r = mod, newr = a;
    while (newr != 0) {
        long long quotient = r / newr;
        tie(t, newt) = make_tuple(newt, t - quotient * newt);
        tie(r, newr) = make_tuple(newr, r - quotient * newr);
    }
    if (r > 1) return -1; // not invertible
    if (t < 0) t += mod;
    return t;
}

// Point addition on elliptic curve
Point add(Point P, Point Q) {
    if (P.at_infinity) return Q;
    if (Q.at_infinity) return P;

    if (P.x == Q.x && (P.y + Q.y) % p == 0) {
        return Point(); // Point at infinity
    }

    long long lambda;
    if (P == Q) {
        long long num = (3 * P.x * P.x + a) % p;
        long long den = modinv(2 * P.y % p, p);
        if (den == -1) return Point();
        lambda = (num * den) % p;
    } else {
        long long num = (Q.y - P.y + p) % p;
        long long den = modinv((Q.x - P.x + p) % p, p);
        if (den == -1) return Point();
        lambda = (num * den) % p;
    }

    long long x3 = (lambda * lambda - P.x - Q.x + p + p) % p;
    long long y3 = (lambda * (P.x - x3 + p) - P.y + p) % p;

    return Point(x3, y3);
}

// Scalar multiplication
Point multiply(Point P, int k) {
    Point result; // point at infinity
    while (k > 0) {
        if (k % 2 == 1) result = add(result, P);
        P = add(P, P);
        k /= 2;
    }
    return result;
}

// Check if a point lies on the curve
bool on_curve(Point P) {
    if (P.at_infinity) return true;
    long long lhs = (P.y * P.y) % p;
    long long rhs = (P.x * P.x * P.x + a * P.x + b) % p;
    return lhs == rhs;
}

int main() {
    Point R(1616, 829);
    vector<Point> matching_points;

    for (long long x = 0; x < p; ++x) {
        long long rhs = (x * x * x + a * x + b) % p;
        for (long long y = 0; y < p; ++y) {
            if ((y * y) % p == rhs) {
                Point Q(x, y);
                if (multiply(Q, 3) == R) {
                    matching_points.push_back(Q);
                }
            }
        }
    }

    cout << "Found " << matching_points.size() << " point(s) Q such that 3Q = R:\n";
    for (auto& pt : matching_points) {
        cout << "(" << pt.x << ", " << pt.y << ")\n";
    }

    return 0;
}
