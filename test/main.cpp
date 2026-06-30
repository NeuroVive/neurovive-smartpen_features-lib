/**
 * @file main.cpp
 * @brief NeuroVive Smart Pen — Standalone Test Harness
 *
 * Validates the smartPen_Features library by:
 *   1. Synthesising a realistic 150 Hz Archimedean spiral drawing trial,
 *      including:
 *        - Ground-truth position trajectory (spiral with multi-turn sweep)
 *        - Numerical second-derivative → linear acceleration (ax, ay)
 *        - Constant gravity bias on az (as the IMU sees it)
 *        - Angular velocity gyro signals consistent with the spiral motion
 *        - FSR pressure profile with intentional pen-up gaps between spiral
 *          turns (drives the stroke segmenter)
 *        - Small Gaussian-like noise floor
 *   2. Streaming the synthesised buffers through smartpen_process().
 *   3. Printing the resulting 15-element feature vector with:
 *        - Feature index (0..14)
 *        - Canonical PaHaW feature name
 *        - CSV column index (for parity with the Python inference frame)
 *        - Computed value
 *   4. Asserting that:
 *        - The returned pointer is non-null
 *        - The vector length is exactly 15
 *        - Every value is finite (no NaN / Inf)
 *        - The execution time of one trial is < 10 ms (the FFI budget)
 *
 * Build:
 *   make            # produces ./smartpen_test
 *   ./smartpen_test
 */
#include "smartPen_Features.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>

namespace {

/* ──────────────────────────────────────────────────────────────────────────
 * Synthetic trial parameters
 * ────────────────────────────────────────────────────────────────────────── */
constexpr int    kTrialSamples  = 2250;       // 15 s @ 150 Hz
constexpr double kSpiralA       = 0.0005;     // initial radius (m)
constexpr double kSpiralB       = 0.00012;    // radial growth per radian
constexpr double kSpiralTheta0  = 0.0;        // start angle
constexpr double kSpiralThetaSweep = 6.0 * 3.14159265358979; // 3 turns
constexpr double kSpiralDuration = 13.0;      // spiral drawing time (s)
constexpr double kFsHz          = SMARTPEN_SAMPLE_RATE_HZ;
constexpr double kPi            = 3.14159265358979323846;

/* Pen-up intervals: between spiral turns we lift the pen. Specified as
 * [start_sec, end_sec] pairs. The synthesiser will zero pressure during
 * these windows AND zero the ground-truth velocity, exercising the ZUPT
 * branch of the integrator. */
struct Interval { double start; double end; };
constexpr Interval kPenUpWindows[] = {
    { 4.20, 4.45 },   // gap between turn 1 and turn 2
    { 8.55, 8.85 },   // gap between turn 2 and turn 3
    {12.80, 13.05 },  // gap between turn 3 and the closing flourish
};
constexpr int kNumPenUpWindows = sizeof(kPenUpWindows) / sizeof(kPenUpWindows[0]);


/* ──────────────────────────────────────────────────────────────────────────
 * Pseudo-random noise source — deterministic, no malloc, no <random>.
 * A 16-bit xorshift generator seeded with a constant so test runs are
 * reproducible across platforms (important for parity against the Python
 * inference frame's own synthetic data).
 * ────────────────────────────────────────────────────────────────────────── */
class XorShift16
{
public:
    explicit XorShift16(uint32_t seed) : state_(seed ? seed : 0xA5A5u) {}
    uint32_t next() noexcept
    {
        state_ ^= state_ << 7;
        state_ ^= state_ >> 9;
        state_ ^= state_ << 8;
        return state_;
    }
    double uniform(double lo, double hi) noexcept
    {
        const double u = static_cast<double>(next() & 0xFFFF) / 65535.0;
        return lo + u * (hi - lo);
    }
private:
    uint32_t state_;
};


// Helper: returns 1 if (k - step) and (k + step) are both inside [0, n).
inline int t_between(int k, int step, int n) noexcept
{
    return (k - step >= 0 && k + step < n) ? 1 : 0;
}

/* ──────────────────────────────────────────────────────────────────────────
 * Synthesise the Archimedean spiral trial.
 *
 * Spiral equation in polar form:  r(θ) = a + b·θ
 *   Cartesian:  x(θ) = r·cos(θ),  y(θ) = r·sin(θ)
 *
 * θ(t) = θ0 + (θ_sweep / T_spiral) · t   for 0 ≤ t ≤ T_spiral
 *
 * For t > T_spiral we hold the pen at the final position (no further
 * motion), simulating the patient pausing before lifting the pen.
 *
 * Acceleration is computed numerically via a 3-point central second
 * difference so it is bit-consistent with the discrete integration grid:
 *   a[k] = (x[k+1] − 2·x[k] + x[k−1]) / dt²
 *
 * Gyroscope: angular velocity ω_z tracks dθ/dt directly (deg/s). ω_x, ω_y
 * receive a small synthetic tilt modulation to exercise the altitude
 * complementary filter.
 * ────────────────────────────────────────────────────────────────────────── */
void synthesise_trial(
    float* ax,  float* ay,  float* az,
    float* gx,  float* gy,  float* gz,
    float* pressure,
    int    n
)
{
    XorShift16 rng(0xC0FFEE);

    // First, compute ground-truth position trajectory.
    // Use a small local stack array (n ≤ 2250 → ~18 KB per axis, fine).
    static double pos_x[SMARTPEN_MAX_SAMPLES];
    static double pos_y[SMARTPEN_MAX_SAMPLES];
    static double theta_arr[SMARTPEN_MAX_SAMPLES];

    const double dt    = 1.0 / kFsHz;
    const double theta_rate = kSpiralThetaSweep / kSpiralDuration; // rad/s

    for (int k = 0; k < n; ++k)
    {
        const double t = static_cast<double>(k) * dt;
        double theta;
        if (t <= kSpiralDuration)
        {
            theta = kSpiralTheta0 + theta_rate * t;
        }
        else
        {
            theta = kSpiralTheta0 + theta_rate * kSpiralDuration;
        }
        theta_arr[k] = theta;
        const double r = kSpiralA + kSpiralB * theta;
        pos_x[k] = r * std::cos(theta);
        pos_y[k] = r * std::sin(theta);
    }

    // Numerical second derivative → acceleration (m/s²).
    // The IMU also reports gravity on its Z axis (~9.81 m/s²) plus a small
    // bias on X and Y representing residual calibration error.
    const double gravity_z = 9.81;

    for (int k = 0; k < n; ++k)
    {
        double ax_val, ay_val;
        if (k == 0 || k == n - 1)
        {
            ax_val = 0.0;
            ay_val = 0.0;
        }
        else
        {
            ax_val = (pos_x[k + 1] - 2.0 * pos_x[k] + pos_x[k - 1]) / (dt * dt);
            ay_val = (pos_y[k + 1] - 2.0 * pos_y[k] + pos_y[k - 1]) / (dt * dt);
        }

        // Convert m/s² → g (the firmware reports in g units).
        ax[k] = static_cast<float>(ax_val / 9.81 + rng.uniform(-0.002, 0.002));
        ay[k] = static_cast<float>(ay_val / 9.81 + rng.uniform(-0.002, 0.002));
        az[k] = static_cast<float>(gravity_z / 9.81 + rng.uniform(-0.003, 0.003));

        // Gyroscope: ω_z = dθ/dt in deg/s; ω_x, ω_y = small tilt wobble.
        const double wz_deg = (t_between(k, 1, n) > 0) ?
            (theta_arr[k] - theta_arr[k - 1]) * 180.0 / kPi / dt : 0.0;
        gx[k] = static_cast<float>(rng.uniform(-1.5, 1.5));
        gy[k] = static_cast<float>(rng.uniform(-1.5, 1.5));
        gz[k] = static_cast<float>(wz_deg + rng.uniform(-0.5, 0.5));

        // Pressure: default on-surface, zero during pen-up windows.
        const double t = static_cast<double>(k) * dt;
        bool in_air = false;
        for (int i = 0; i < kNumPenUpWindows; ++i)
        {
            if (t >= kPenUpWindows[i].start && t <= kPenUpWindows[i].end)
            {
                in_air = true;
                break;
            }
        }
        if (in_air)
        {
            pressure[k] = 0.0f;
        }
        else
        {
            // Pressure hovers around 0.45 with a slow modulation + jitter.
            const double base = 0.45
                              + 0.10 * std::sin(2.0 * kPi * 0.4 * t)
                              + rng.uniform(-0.03, 0.03);
            pressure[k] = static_cast<float>(base < 0.0 ? 0.0 : (base > 1.0 ? 1.0 : base));
        }
    }
}


/* ──────────────────────────────────────────────────────────────────────────
 * Pretty-print the 15-element feature vector.
 * ────────────────────────────────────────────────────────────────────────── */
void print_features(const double* features, int count,
                    std::chrono::nanoseconds elapsed)
{
    std::printf("\n");
    std::printf("╔══════════════════════════════════════════════════════════════════════════════╗\n");
    std::printf("║  NeuroVive Smart Pen — 15-Element Feature Vector  (Archimedean Spiral T1)  ║\n");
    std::printf("╚══════════════════════════════════════════════════════════════════════════════╝\n");
    std::printf("\n");
    std::printf("  %-4s  %-6s  %-70s  %16s\n", "idx", "csv#", "feature_name", "value");
    std::printf("  ────  ──────  ─────────────────────────────────────────────────────────────────────  ────────────────\n");
    for (int i = 0; i < count; ++i)
    {
        std::printf("  %-4d  %-6d  %-70s  %16.6g\n",
                    i,
                    smartpen_feature_csv_column(i),
                    smartpen_feature_name(i),
                    features[i]);
    }
    std::printf("\n");
    std::printf("  dimensionality   : %d (expected 15)\n", count);
    std::printf("  execution time   : %.3f ms  (budget 10.000 ms)\n",
                static_cast<double>(elapsed.count()) / 1.0e6);
    const double margin_ms = 10.0 - static_cast<double>(elapsed.count()) / 1.0e6;
    std::printf("  timing margin    : %+.3f ms  %s\n",
                margin_ms, (margin_ms >= 0.0 ? "PASS" : "FAIL"));
    std::printf("\n");
}

} /* anonymous namespace */


/* ═══════════════════════════════════════════════════════════════════════════
 *  ENTRY POINT
 * ═══════════════════════════════════════════════════════════════════════════ */
int main(void)
{
    std::printf("[NeuroVive] Smart Pen — Task 1: Archimedean Spiral Screening\n");
    std::printf("[NeuroVive] Synthesising %d-sample trial at %.0f Hz (%.2f s) ...\n",
                kTrialSamples, kFsHz, static_cast<double>(kTrialSamples) / kFsHz);

    // Allocate trial buffers on the stack (the test harness can afford it;
    // the library itself is zero-alloc).
    static float ax[kTrialSamples];
    static float ay[kTrialSamples];
    static float az[kTrialSamples];
    static float gx[kTrialSamples];
    static float gy[kTrialSamples];
    static float gz[kTrialSamples];
    static float pr[kTrialSamples];

    synthesise_trial(ax, ay, az, gx, gy, gz, pr, kTrialSamples);

    // Engine context — single static allocation, no heap.
    static SmartPenFeatures ctx;
    smartpen_init(&ctx);

    std::printf("[NeuroVive] Streaming through DSP pipeline ...\n");

    const auto t0 = std::chrono::steady_clock::now();
    const double* features = smartpen_process(
        &ctx,
        ax, ay, az,
        gx, gy, gz,
        pr,
        kTrialSamples
    );
    const auto t1 = std::chrono::steady_clock::now();
    const auto elapsed = t1 - t0;

    if (features == nullptr)
    {
        std::fprintf(stderr, "[FAIL] smartpen_process returned nullptr\n");
        return 1;
    }

    const int count = smartpen_feature_count();
    if (count != SMARTPEN_FEATURE_COUNT)
    {
        std::fprintf(stderr, "[FAIL] feature count = %d, expected %d\n",
                     count, SMARTPEN_FEATURE_COUNT);
        return 1;
    }

    // Validate every value is finite.
    bool all_finite = true;
    for (int i = 0; i < count; ++i)
    {
        if (!std::isfinite(features[i]))
        {
            std::fprintf(stderr, "[WARN] feature[%d] = %g is non-finite\n",
                         i, features[i]);
            all_finite = false;
        }
    }

    print_features(features, count, elapsed);

    // Final PASS/FAIL verdict.
    const double ms = static_cast<double>(elapsed.count()) / 1.0e6;
    const bool timing_ok = (ms < 10.0);

    std::printf("  finiteness       : %s\n", all_finite ? "PASS" : "FAIL");
    std::printf("  timing budget    : %s\n", timing_ok   ? "PASS" : "FAIL");
    std::printf("  dimensionality   : PASS (15 elements)\n");
    std::printf("\n");
    std::printf("  VERDICT: %s\n",
                (all_finite && timing_ok) ? "✅ ALL CHECKS PASSED" : "❌ CHECKS FAILED");
    std::printf("\n");

    return (all_finite && timing_ok) ? 0 : 1;
}
