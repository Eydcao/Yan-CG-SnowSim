#include <cmath>
#include <iostream>
#include <random>

#undef M_PI
#define M_PI 3.141592653589793f

static const float EPSILON = 0.001;
static const float kInfinity = std::numeric_limits<float>::max();
static const float FRAMERATE = 60.;
// Adaptive timestep adjustment
static const float CFL = .04;
static const float MAX_TIMESTEP = 5.e-4;
// Percentage that FLIP TAKES, PIC is 1 - that
static const float FLIP_PERCENT = .95;
// Percentage that should be implicit vs explicit
static const float IMPLICIT_RATIO = 0;
static const int MAX_IMPLICIT_ITERS = 30;
static const float MAX_IMPLICIT_ERR = 1.e4;
static const float MIN_IMPLICIT_ERR = 1.e-4;
static const float GRAVITY = -9.8;

inline float clamp(const float &lo, const float &hi, const float &v)
{
    return std::max(lo, std::min(hi, v));
}

inline bool solveQuadratic(const float &a, const float &b, const float &c,
                           float &x0, float &x1)
{
    float discr = b * b - 4 * a * c;
    if (discr < 0)
        return false;
    else if (discr == 0)
        x0 = x1 = -0.5 * b / a;
    else
    {
        float q = (b > 0) ? -0.5 * (b + sqrt(discr)) : -0.5 * (b - sqrt(discr));
        x0 = q / a;
        x1 = c / q;
    }
    if (x0 > x1) std::swap(x0, x1);
    return true;
}

// inline float get_random_float()
// {
//     std::random_device dev;
//     std::mt19937 rng(dev());
//     std::uniform_real_distribution<float> dist(0.f, 1.f); // distribution in
//     range [1, 6]

//     return dist(rng);
// }

// inline void UpdateProgress(float progress)
// {
//     int barWidth = 70;

//     std::cout << "[";
//     int pos = barWidth * progress;
//     for (int i = 0; i < barWidth; ++i) {
//         if (i < pos) std::cout << "=";
//         else if (i == pos) std::cout << ">";
//         else std::cout << " ";
//     }
//     std::cout << "] " << int(progress * 100.0) << " %\r";
//     std::cout.flush();
// };
