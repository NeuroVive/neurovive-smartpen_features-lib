/**
 * @file smartPen_Features.cpp
 * @brief NeuroVive Smart Pen — Zero-allocation, deterministic DSP pipeline.
 *
 * This file implements every DSP block from scratch:
 *   (A) 2nd-order high-pass Butterworth filter (drift suppression)
 *   (B) ZUPT-gated double integration (acceleration → velocity → position)
 *   (C) Complementary filter (gyro + accelerometer fusion → azimuth/altitude)
 *   (D) Pressure-threshold stroke segmentation
 *   (E) Statistical moments: mean, variance, skewness, kurtosis
 *   (F) QuickSelect O(N) median / percentile partitioner
 *   (G) Teager-Kaiser Energy operator  ψ[x[n]] = x²[n] − x[n−1]·x[n+1]
 *   (H) Single-pole IIR intrinsic-mode filter (primary dynamic component)
 *   (I) The 15-feature assembler with strict index alignment
 *
 * MEMORY: zero heap allocation. Every buffer is owned by the caller-supplied
 * SmartPenFeatures struct (file-scope stack on Dart isolate). No std::vector,
 * no `new`, no `malloc`. All loops are bounded by SMARTPEN_MAX_SAMPLES so the
 * worst-case execution time is statically computable.
 *
 * REFERENCES:
 *   [1] شرح القلم.md (hardware mechanics, packet map, ZUPT rationale)
 *   [2] 1_1.csv (PaHaW-derived 870-column feature schema — used to lock the
 *       exact nomenclature of the 15 golden features)
 *   [3] PaHaW public repository (github.com/musaru/PD_PaHaW) — kinematics
 *       extraction standards.
 */
#include "smartPen_Features.h"

#include <cmath>
#include <cstdint>
#include <cstring>

/* The DSP core is compiled as C++ but exposes a pure C ABI. We do not pull
 * in any STL container — only <cmath> primitives. */
namespace {

/* ─────────────────────────────────────────────────────────────────────────
 * (A) 2nd-ORDER HIGH-PASS BUTTERWORTH FILTER (Direct-Form II Transposed)
 *
 * Design:
 *   Cut-off fc = SMARTPEN_HPF_CUTOFF_HZ (= 0.5 Hz) at fs = 150 Hz.
 *   Pre-warp:  Ω = tan(π·fc / fs)
 *   Analog prototype (normalised 2nd-order Butterworth):
 *       H(s) = 1 / (s² + √2·s + 1)
 *   Bilinear high-pass transform s ⇒ (1 + s')/(1 − s') yields, after
 *   algebra:
 *       b0 =  1 / (1 + √2·Ω + Ω²)
 *       b1 = -2·b0
 *       b2 =  b0
 *       a1 =  2·(Ω² − 1) / (1 + √2·Ω + Ω²)
 *       a2 =  (1 − √2·Ω + Ω²) / (1 + √2·Ω + Ω²)
 *
 * The coefficients below are pre-computed at compile time using constexpr
 * so the runtime cost is zero divisions.
 * ───────────────────────────────────────────────────────────────────────── */
struct ButterworthHPF
{
    double b0, b1, b2;   // numerator coefficients
    double a1, a2;       // denominator coefficients (a0 normalised to 1)
    double w1, w2;       // DF-II transposed state
};

constexpr ButterworthHPF make_butterworth_hpf() noexcept
{
    constexpr double fs = SMARTPEN_SAMPLE_RATE_HZ;
    constexpr double fc = SMARTPEN_HPF_CUTOFF_HZ;
    constexpr double pi = 3.14159265358979323846;
    constexpr double omega = std::tan(pi * fc / fs);             // pre-warped
    constexpr double sqrt2 = 1.41421356237309504880;
    constexpr double denom = 1.0 + sqrt2 * omega + omega * omega;
    return ButterworthHPF{
        /* b0 */  1.0 / denom,
        /* b1 */ -2.0 / denom,
        /* b2 */  1.0 / denom,
        /* a1 */  2.0 * (omega * omega - 1.0) / denom,
        /* a2 */ (1.0 - sqrt2 * omega + omega * omega) / denom,
        /* w1 */ 0.0,
        /* w2 */ 0.0
    };
}

/** Reset filter state (call before each trial). */
inline void hpf_reset(ButterworthHPF& f) noexcept { f.w1 = 0.0; f.w2 = 0.0; }

/** Process one sample through the high-pass filter (in-place semantics). */
inline double hpf_step(ButterworthHPF& f, double x) noexcept
{
    // Direct-Form II Transposed:
    //   y[n]   = b0·x[n] + w1
    //   w1     = b1·x[n] - a1·y[n] + w2
    //   w2     = b2·x[n] - a2·y[n]
    const double y = f.b0 * x + f.w1;
    f.w1 = f.b1 * x - f.a1 * y + f.w2;
    f.w2 = f.b2 * x - f.a2 * y;
    return y;
}


/* ─────────────────────────────────────────────────────────────────────────
 * (B) ZUPT-GATED DOUBLE INTEGRATION
 *
 * Standard kinematic chain:  a(t) →∫→ v(t) →∫→ p(t)
 *
 * Trapezoidal integrator per sample (k):
 *   v[k]   = v[k-1] + 0.5·(a[k-1] + a[k])·dt
 *   p[k]   = p[k-1] + 0.5·(v[k-1] + v[k])·dt
 *
 * ZUPT (Zero-Velocity Update, شرح القلم.md §1):
 *   When the FSR reports the pen is in air (P ≤ SMARTPEN_ZUPT_GATE), the
 *   physical velocity of the pen tip is zero. We forcibly clamp v[k] = 0
 *   AND freeze position integration. This bounds drift to the on-surface
 *   writing intervals only, eliminating the quadratic E ∝ t² growth that
 *   plagues ungated inertial odometry.
 * ───────────────────────────────────────────────────────────────────────── */
struct IntegratorState
{
    double v_prev;   // v[k-1]
    double a_prev;   // a[k-1]
    double p_prev;   // p[k-1]
};

inline void integ_reset(IntegratorState& s) noexcept
{
    s.v_prev = 0.0;
    s.a_prev = 0.0;
    s.p_prev = 0.0;
}

/** Advance one sample. `pen_down` is true iff pressure > threshold.
 *  Returns the new velocity; new position is written through *pos_out. */
inline double integ_step(IntegratorState& s, double a, double dt, bool pen_down,
                        double* pos_out) noexcept
{
    double v_new;
    double p_new;
    if (pen_down)
    {
        // Trapezoidal integration while writing on the surface.
        v_new = s.v_prev + 0.5 * (s.a_prev + a) * dt;
        p_new = s.p_prev + 0.5 * (s.v_prev + v_new) * dt;
    }
    else
    {
        // ZUPT clamp: pen is in air, so velocity is physically zero and
        // position holds its last on-surface value (no further drift).
        v_new = 0.0;
        p_new = s.p_prev;
    }
    s.a_prev = a;
    s.v_prev = v_new;
    s.p_prev = p_new;
    *pos_out = p_new;
    return v_new;
}


/* ─────────────────────────────────────────────────────────────────────────
 * (C) COMPLEMENTARY FILTER — azimuth & altitude
 *
 * Reference (شرح القلم.md §2):
 *   θ̂[t] = α · (θ̂[t-1] + ω_gyro · dt) + (1 − α) · θ_acc
 *
 * where:
 *   α        = SMARTPEN_COMPLEMENTARY_ALPHA = 0.98  (gyro high-pass branch)
 *   ω_gyro   = angular rate from the MPU6050 gyroscope (deg/s)
 *   θ_acc    = arctan(acc_horizontal / acc_vertical)  (accelerometer tilt)
 *
 * Azimuth  : rotation about the vertical (Z) axis — driven by gyro_z.
 * Altitude : tilt above the horizontal plane — driven by gyro_y and the
 *            arctan2(acc_x, acc_z) low-pass observation.
 * angle_with_1 : signed tilt of the pen axis relative to its initial
 *                orientation — computed as the absolute deviation of
 *                altitude from altitude[0]. Used by feature 9
 *                (angle_data_with_1_mean).
 * ───────────────────────────────────────────────────────────────────────── */
struct ComplementaryState
{
    double azimuth_prev;   // degrees
    double altitude_prev;  // degrees
    double altitude_init;  // first-sample altitude reference
    bool   initialised;
};

inline void comp_reset(ComplementaryState& s) noexcept
{
    s.azimuth_prev  = 0.0;
    s.altitude_prev = 0.0;
    s.altitude_init = 0.0;
    s.initialised   = false;
}

inline void comp_step(ComplementaryState& s,
                      double gyro_y, double gyro_z,
                      double acc_x, double acc_z, double dt,
                      double* azimuth_out, double* altitude_out,
                      double* angle_with_1_out) noexcept
{
    // Accelerometer-derived tilt observation (low-frequency reference).
    // atan2 chosen over atan to remain stable near the vertical singularity.
    const double acc_tilt = std::atan2(acc_x, acc_z) * 180.0 / 3.14159265358979323846;

    // Gyroscope high-frequency branch: integrate rate.
    const double gyro_azimuth_branch  = s.azimuth_prev  + gyro_z * dt;
    const double gyro_altitude_branch = s.altitude_prev + gyro_y * dt;

    // Complementary fusion.
    double az_new  = SMARTPEN_COMPLEMENTARY_ALPHA * gyro_azimuth_branch
                   + (1.0 - SMARTPEN_COMPLEMENTARY_ALPHA) * 0.0; // acc has no yaw ref → gyro only
    double alt_new = SMARTPEN_COMPLEMENTARY_ALPHA * gyro_altitude_branch
                   + (1.0 - SMARTPEN_COMPLEMENTARY_ALPHA) * acc_tilt;

    if (!s.initialised)
    {
        s.altitude_init = alt_new;
        s.initialised   = true;
    }

    // angle_with_1 = deviation of altitude from the initial sample's altitude.
    const double ang_w1 = alt_new - s.altitude_init;

    s.azimuth_prev  = az_new;
    s.altitude_prev = alt_new;

    *azimuth_out        = az_new;
    *altitude_out       = alt_new;
    *angle_with_1_out   = ang_w1;
}


/* ─────────────────────────────────────────────────────────────────────────
 * (D) PRESSURE-THRESHOLD STROKE SEGMENTATION
 *
 * Rule (شرح القلم.md §1):
 *   P >  0.05  ⇒ pen-down (on-surface writing)
 *   P ≤  0.05  ⇒ pen-up   (in-air transition)
 *   Each pen-up → pen-down edge opens a new stroke; the trailing pen-down →
 *   pen-up edge closes it. A stroke shorter than 2 samples is rejected as
 *   a mechanical bounce.
 *
 * Output: ctx->stroke_start[] / ctx->stroke_end[] sample-index pairs plus
 * ctx->n_strokes.
 * ───────────────────────────────────────────────────────────────────────── */
void segment_strokes(SmartPenFeatures* ctx) noexcept
{
    ctx->n_strokes = 0;
    bool in_stroke = false;
    int  stroke_open = 0;

    for (int i = 0; i < ctx->n_samples; ++i)
    {
        const bool pen_down = ctx->pressure[i] > SMARTPEN_PRESSURE_THRESHOLD;
        if (pen_down && !in_stroke)
        {
            // Rising edge — open a new stroke.
            stroke_open = i;
            in_stroke   = true;
        }
        else if (!pen_down && in_stroke)
        {
            // Falling edge — close the current stroke.
            const int len = i - stroke_open;
            if (len >= 2 && ctx->n_strokes < SMARTPEN_MAX_STROKES)
            {
                ctx->stroke_start[ctx->n_strokes] = stroke_open;
                ctx->stroke_end  [ctx->n_strokes] = i;
                ++ctx->n_strokes;
            }
            in_stroke = false;
        }
    }
    // Close any stroke that runs to the end of the buffer.
    if (in_stroke)
    {
        const int len = ctx->n_samples - stroke_open;
        if (len >= 2 && ctx->n_strokes < SMARTPEN_MAX_STROKES)
        {
            ctx->stroke_start[ctx->n_strokes] = stroke_open;
            ctx->stroke_end  [ctx->n_strokes] = ctx->n_samples;
            ++ctx->n_strokes;
        }
    }
}


/* ─────────────────────────────────────────────────────────────────────────
 * (E) STATISTICAL MOMENTS — non-allocating, single-pass where possible
 *
 *   mean      = (1/N) · Σ x
 *   variance  = (1/N) · Σ (x − mean)²           (population, N denominator)
 *   m2        = variance
 *   m3        = (1/N) · Σ (x − mean)³
 *   m4        = (1/N) · Σ (x − mean)⁴
 *   skewness  = m3 / m2^1.5
 *   kurtosis  = m4 / m2² − 3                    (excess kurtosis)
 * ───────────────────────────────────────────────────────────────────────── */
struct Moments
{
    double mean;
    double variance;
    double skewness;
    double kurtosis;   // excess
};

inline Moments compute_moments(const double* x, int n) noexcept
{
    Moments m{0.0, 0.0, 0.0, 0.0};
    if (n <= 0) return m;

    // Pass 1: mean.
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += x[i];
    m.mean = sum / static_cast<double>(n);

    // Pass 2: central moments m2, m3, m4.
    double s2 = 0.0, s3 = 0.0, s4 = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const double d  = x[i] - m.mean;
        const double d2 = d * d;
        s2 += d2;
        s3 += d2 * d;
        s4 += d2 * d2;
    }
    const double inv = 1.0 / static_cast<double>(n);
    m.variance = s2 * inv;
    const double m2 = m.variance;
    const double m3 = s3 * inv;
    const double m4 = s4 * inv;

    // Guarded ratio — if variance is ~0 the signal is constant and the
    // higher moments are undefined; we return 0 rather than NaN.
    if (m2 > 1e-12)
    {
        m.skewness = m3 / std::pow(m2, 1.5);
        m.kurtosis = m4 / (m2 * m2) - 3.0;
    }
    else
    {
        m.skewness = 0.0;
        m.kurtosis = 0.0;
    }
    return m;
}

/** Arithmetic mean only — used when moments aren't needed. */
inline double compute_mean(const double* x, int n) noexcept
{
    if (n <= 0) return 0.0;
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += x[i];
    return s / static_cast<double>(n);
}


/* ─────────────────────────────────────────────────────────────────────────
 * (F) QUICKSELECT — O(N) median / percentile partitioner
 *
 * Hoare-partition QuickSelect over a *scratch copy* of the input so the
 * caller's source array is never mutated. The scratch buffer is supplied
 * by the caller (statically allocated). Worst case is O(N²) but the
 * median-of-three pivot makes the expected case O(N), and crucially the
 * recursion is replaced with an iterative loop so there is no stack growth.
 *
 * Percentile convention (matches numpy.percentile, linear interpolation):
 *   rank = p/100 · (N − 1)
 *   lo   = floor(rank), hi = ceil(rank)
 *   result = x[lo] + (rank − lo) · (x[hi] − x[lo])
 * ───────────────────────────────────────────────────────────────────────── */
inline void swap_double(double& a, double& b) noexcept
{
    const double t = a; a = b; b = t;
}

double quickselect(double* arr, int n, int k) noexcept
{
    // Iterative partition — no recursion, no stack growth.
    int lo = 0;
    int hi = n - 1;
    while (lo < hi)
    {
        // Median-of-three pivot selection.
        const int mid = lo + (hi - lo) / 2;
        if (arr[mid] < arr[lo]) swap_double(arr[mid], arr[lo]);
        if (arr[hi]  < arr[lo]) swap_double(arr[hi],  arr[lo]);
        if (arr[mid] < arr[hi]) swap_double(arr[mid], arr[hi]);
        const double pivot = arr[hi];

        int i = lo - 1;
        for (int j = lo; j < hi; ++j)
        {
            if (arr[j] <= pivot)
            {
                ++i;
                swap_double(arr[i], arr[j]);
            }
        }
        swap_double(arr[i + 1], arr[hi]);
        const int p = i + 1;

        if (p == k) return arr[p];
        if (p <  k) lo = p + 1;
        else        hi = p - 1;
    }
    return arr[lo];
}

double percentile(const double* src, int n, double p,
                  double* scratch) noexcept
{
    if (n <= 0) return 0.0;
    // Copy into scratch — does not mutate caller's data.
    for (int i = 0; i < n; ++i) scratch[i] = src[i];

    if (n == 1) return scratch[0];

    // Linear-interpolation rank (numpy default).
    double rank = (p / 100.0) * static_cast<double>(n - 1);
    if (rank < 0.0)                rank = 0.0;
    if (rank > static_cast<double>(n - 1)) rank = static_cast<double>(n - 1);

    const int lo = static_cast<int>(std::floor(rank));
    const int hi = static_cast<int>(std::ceil(rank));

    const double v_lo = quickselect(scratch, n, lo);
    // After quickselect for lo, element at index hi may have moved — re-select
    // on a fresh copy for correctness (still O(N)).
    for (int i = 0; i < n; ++i) scratch[i] = src[i];
    const double v_hi = (hi == lo) ? v_lo : quickselect(scratch, n, hi);

    return v_lo + (rank - static_cast<double>(lo)) * (v_hi - v_lo);
}

inline double median(const double* src, int n, double* scratch) noexcept
{
    return percentile(src, n, 50.0, scratch);
}


/* ─────────────────────────────────────────────────────────────────────────
 * (G) TEAGER-KAISER ENERGY OPERATOR  (TKE / TKEO)
 *
 *   ψ[x[n]] = x²[n] − x[n−1] · x[n+1]
 *
 * The operator localises the instantaneous energy of an oscillating signal
 * and is widely used in PD handwriting analysis because it amplifies the
 * tremor band (~4–7 Hz) relative to the intentional motion band.
 *
 * Global TKE feature: we apply ψ over the whole array (edges excluded) and
 * return the mean — matching the Python pipeline's `teager_kaiser_energy(_, 1)`
 * signature where the trailing `1` denotes the mean-aggregation mode.
 * ───────────────────────────────────────────────────────────────────────── */
double teager_kaiser_energy_mean(const double* x, int n) noexcept
{
    if (n < 3) return 0.0;
    double acc = 0.0;
    int   cnt = 0;
    for (int i = 1; i < n - 1; ++i)
    {
        const double psi = x[i] * x[i] - x[i - 1] * x[i + 1];
        acc += psi;
        ++cnt;
    }
    return (cnt > 0) ? (acc / static_cast<double>(cnt)) : 0.0;
}


/* ─────────────────────────────────────────────────────────────────────────
 * (H) SINGLE-POLE IIR INTRINSIC-MODE FILTER
 *
 *   y[n] = α · x[n] + (1 − α) · y[n−1]   with α = SMARTPEN_INTRINSIC_ALPHA
 *
 * This is the "single-pass digital filter" mandated by the brief. It acts
 * as a low-pass on the velocity signal, isolating the primary dynamic mode
 * (the intentional handwriting component, ~0.5–3 Hz) from the high-frequency
 * sensor noise. The TKE of the filtered signal is then taken as the
 * intrinsic_feature_teager_kaiser_energy(velocity, 1) feature.
 *
 * The filter is one-pole (no z² term) so it cannot ring — a deliberate
 * choice for embedded determinism.
 * ───────────────────────────────────────────────────────────────────────── */
void intrinsic_filter(const double* src, int n, double* dst) noexcept
{
    if (n <= 0) return;
    const double a = SMARTPEN_INTRINSIC_ALPHA;
    dst[0] = a * src[0];
    for (int i = 1; i < n; ++i)
    {
        dst[i] = a * src[i] + (1.0 - a) * dst[i - 1];
    }
}


/* ─────────────────────────────────────────────────────────────────────────
 * (I) STATIC SCRATCH BUFFERS for percentile computation
 *
 * Sized to SMARTPEN_MAX_SAMPLES so they can host any sub-trial slice. They
 * are file-scope statics with internal linkage — zero heap cost, fixed
 * addresses, no contention between Dart isolates because each isolate
 * loads its own .so image in standard FFI deployments. (If you re-enter
 * the library from the same isolate, the buffers are overwritten between
 * calls which is fine because features are computed sequentially.)
 * ───────────────────────────────────────────────────────────────────────── */
double g_percentile_scratch[SMARTPEN_MAX_SAMPLES];


/* ─────────────────────────────────────────────────────────────────────────
 * (J) CANONICAL FEATURE-NAME / CSV-COLUMN LOOKUP TABLES
 *
 * The names below are character-for-character identical to the PaHaW-derived
 * 1_1.csv header so that the C++ output column-locks with the Python
 * inference frame. The CSV column index is the 1-based position in the
 * header (after the leading unnamed index column).
 * ───────────────────────────────────────────────────────────────────────── */
struct FeatureMeta
{
    const char* name;
    int         csv_column;
};

constexpr FeatureMeta g_feature_meta[SMARTPEN_FEATURE_COUNT] = {
    { "teager_kaiser_energy(displacement, 1)",                                 185 },
    { "stroke_pressure_1st_each_stroke_mean",                                  628 },
    { "teager_kaiser_energy(velocity, 1)",                                     187 },
    { "stroke_pressure_1st_each_stroke_min",                                   627 },
    { "stroke_pressure_1st_each_stroke_1st_percentile",                        632 },
    { "intrinsic_feature_teager_kaiser_energy(velocity, 1)",                   231 },
    { "y_velocity_99th_percentile",                                            110 },
    { "last_azimuth_skewness",                                                 432 },
    { "stroke_pressure_1st_each_stroke_median",                                629 },
    { "angle_data_with_1_mean",                                                556 },
    { "y_velocity_displacement_between_99th_percentile_and_1st_percentile",    111 },
    { "stroke_dis_median_each_stroke_max",                                     698 },
    { "first_x_median",                                                        236 },
    { "x_displacement_kurtosis",                                                86 },
    { "last_y_velocity_no_abs_median",                                         547 },
};

} /* anonymous namespace */


/* ═══════════════════════════════════════════════════════════════════════════
 *  PUBLIC C API IMPLEMENTATION
 * ═══════════════════════════════════════════════════════════════════════════ */

extern "C" {

void smartpen_init(SmartPenFeatures* ctx) noexcept
{
    if (ctx == nullptr) return;
    std::memset(ctx, 0, sizeof(SmartPenFeatures));
}


const double* smartpen_process(
    SmartPenFeatures* ctx,
    const float* acc_x,
    const float* acc_y,
    const float* acc_z,
    const float* gyro_x,
    const float* gyro_y,
    const float* gyro_z,
    const float* pressure,
    int          n_samples
) noexcept
{
    if (ctx == nullptr || acc_x == nullptr || acc_y == nullptr || acc_z == nullptr
        || gyro_x == nullptr || gyro_y == nullptr || gyro_z == nullptr
        || pressure == nullptr)
    {
        return nullptr;
    }
    if (n_samples <= 0) return nullptr;
    if (n_samples > SMARTPEN_MAX_SAMPLES) n_samples = SMARTPEN_MAX_SAMPLES;

    /* ─── 0. Reset & mirror raw inputs into the static engine buffers ─── */
    smartpen_init(ctx);
    ctx->n_samples = n_samples;
    for (int i = 0; i < n_samples; ++i)
    {
        ctx->acc_x[i]    = acc_x[i];
        ctx->acc_y[i]    = acc_y[i];
        ctx->acc_z[i]    = acc_z[i];
        ctx->gyro_x[i]   = gyro_x[i];
        ctx->gyro_y[i]   = gyro_y[i];
        ctx->gyro_z[i]   = gyro_z[i];
        ctx->pressure[i] = pressure[i];
    }

    /* ─── 1. Build filter / integrator state ─── */
    ButterworthHPF   hpf_x = make_butterworth_hpf();
    ButterworthHPF   hpf_y = make_butterworth_hpf();
    IntegratorState  isx, isy;
    ComplementaryState cs;
    hpf_reset(hpf_x); hpf_reset(hpf_y);
    integ_reset(isx); integ_reset(isy);
    comp_reset(cs);

    /* ─── 2. Per-sample DSP loop ───
     * For each sample k:
     *   (a) high-pass filter the X & Y acceleration channels to kill DC drift
     *       and any residual gravity bias left by the firmware calibration;
     *   (b) ZUPT-gated trapezoidal double integration → vx, vy, x_disp, y_disp;
     *   (c) complementary filter → azimuth, altitude, angle_with_1;
     *   (d) accumulate resultant planar velocity & displacement.
     */
    constexpr double dt = SMARTPEN_DT_SEC;
    for (int k = 0; k < n_samples; ++k)
    {
        const double ax = static_cast<double>(ctx->acc_x[k]);
        const double ay = static_cast<double>(ctx->acc_y[k]);
        const double az = static_cast<double>(ctx->acc_z[k]);
        const double gx = static_cast<double>(ctx->gyro_x[k]);
        const double gy = static_cast<double>(ctx->gyro_y[k]);
        const double gz = static_cast<double>(ctx->gyro_z[k]);
        const bool   pen_down = ctx->pressure[k] > SMARTPEN_PRESSURE_THRESHOLD;
        // gyro_x (roll about the pen's longitudinal axis) is intentionally
        // not consumed by the azimuth/altitude complementary filter — it has
        // no observable effect on planar handwriting kinematics. We retain
        // it in the buffer for parity with the firmware packet map and
        // suppress the unused-variable warning explicitly.
        (void)gx;

        // (a) High-pass-filtered accelerations.
        const double ax_hp = hpf_step(hpf_x, ax);
        const double ay_hp = hpf_step(hpf_y, ay);

        // (b) Double integration with ZUPT.
        double x_pos = 0.0, y_pos = 0.0;
        const double vx = integ_step(isx, ax_hp, dt, pen_down, &x_pos);
        const double vy = integ_step(isy, ay_hp, dt, pen_down, &y_pos);

        ctx->vx[k]     = vx;
        ctx->vy[k]     = vy;
        ctx->x_disp[k] = x_pos;
        ctx->y_disp[k] = y_pos;

        // Resultant planar quantities.
        ctx->velocity[k]     = std::sqrt(vx * vx + vy * vy);
        ctx->displacement[k] = std::sqrt(x_pos * x_pos + y_pos * y_pos);

        // (c) Complementary filter.
        double azim, alt, ang_w1;
        comp_step(cs, gy, gz, ax, az, dt, &azim, &alt, &ang_w1);
        ctx->azimuth[k]      = azim;
        ctx->altitude[k]     = alt;
        ctx->angle_with_1[k] = ang_w1;
    }

    /* ─── 3. Intrinsic-mode velocity filter ─── */
    intrinsic_filter(ctx->velocity, n_samples, ctx->intrinsic_velocity);

    /* ─── 4. Stroke segmentation ─── */
    segment_strokes(ctx);

    /* ─── 5. Per-stroke statistics ───
     * We compute, per stroke:
     *   - median of the resultant displacement (used by feature 11)
     *   - first-derivative of pressure dP/dt at each interior sample,
     *     pooled globally into ctx->dPdt_pool for features 1, 3, 4, 8.
     * dP/dt uses the same trapezoidal convention as the kinematic chain:
     *   dP[k] = (P[k] − P[k−1]) / dt
     */
    ctx->dPdt_count = 0;
    for (int s = 0; s < ctx->n_strokes; ++s)
    {
        const int s0 = ctx->stroke_start[s];
        const int s1 = ctx->stroke_end[s];
        const int len = s1 - s0;
        if (len < 2) { ctx->stroke_dis_median[s] = 0.0; continue; }

        // Stroke-local median of resultant displacement.
        ctx->stroke_dis_median[s] =
            median(&ctx->displacement[s0], len, g_percentile_scratch);

        // dP/dt across interior of stroke.
        for (int k = s0 + 1; k < s1; ++k)
        {
            if (ctx->dPdt_count >= SMARTPEN_MAX_SAMPLES) break;
            const double dp = (static_cast<double>(ctx->pressure[k])
                             - static_cast<double>(ctx->pressure[k - 1])) / dt;
            ctx->dPdt_pool[ctx->dPdt_count++] = dp;
        }
    }

    /* ═══════════════════════════════════════════════════════════════════
     * 6. FEATURE ASSEMBLER — strict 15-element order
     *    Every line below documents the exact mathematical formula mapped
     *    to the corresponding PaHaW CSV column.
     * ═══════════════════════════════════════════════════════════════════ */

    /* ── Feature 0  (CSV #185) ─ teager_kaiser_energy(displacement, 1) ──
     * ψ[d[n]] = d²[n] − d[n−1]·d[n+1]  →  mean over n ∈ [1, N−2]
     * d[n] is the resultant planar displacement magnitude at sample n.
     */
    ctx->features[SMARTPEN_F_TEAGER_KAISER_DISPLACEMENT] =
        teager_kaiser_energy_mean(ctx->displacement, n_samples);

    /* ── Feature 1  (CSV #628) ─ stroke_pressure_1st_each_stroke_mean ──
     * Mean of the first-derivative-of-pressure dP/dt pool, where the pool
     * aggregates all interior dP/dt samples across every segmented stroke.
     */
    ctx->features[SMARTPEN_F_STROKE_PRESSURE_1ST_MEAN] =
        compute_mean(ctx->dPdt_pool, ctx->dPdt_count);

    /* ── Feature 2  (CSV #187) ─ teager_kaiser_energy(velocity, 1) ──
     * ψ[v[n]] = v²[n] − v[n−1]·v[n+1]  →  mean.
     */
    ctx->features[SMARTPEN_F_TEAGER_KAISER_VELOCITY] =
        teager_kaiser_energy_mean(ctx->velocity, n_samples);

    /* ── Feature 3  (CSV #627) ─ stroke_pressure_1st_each_stroke_min ──
     * Global minimum of the dP/dt pool across all strokes.
     */
    {
        double dmin = (ctx->dPdt_count > 0) ? ctx->dPdt_pool[0] : 0.0;
        for (int i = 1; i < ctx->dPdt_count; ++i)
            if (ctx->dPdt_pool[i] < dmin) dmin = ctx->dPdt_pool[i];
        ctx->features[SMARTPEN_F_STROKE_PRESSURE_1ST_MIN] = dmin;
    }

    /* ── Feature 4  (CSV #632) ─ stroke_pressure_1st_each_stroke_1st_percentile ──
     * 1st percentile of the dP/dt pool (QuickSelect, linear interpolation).
     */
    ctx->features[SMARTPEN_F_STROKE_PRESSURE_1ST_1PC] =
        percentile(ctx->dPdt_pool, ctx->dPdt_count, 1.0, g_percentile_scratch);

    /* ── Feature 5  (CSV #231) ─ intrinsic_feature_teager_kaiser_energy(velocity, 1) ──
     * TKE mean of the *intrinsic* (single-pole IIR filtered) velocity signal.
     */
    ctx->features[SMARTPEN_F_INTRINSIC_TKE_VELOCITY] =
        teager_kaiser_energy_mean(ctx->intrinsic_velocity, n_samples);

    /* ── Feature 6  (CSV #110) ─ y_velocity_99th_percentile ──
     * 99th percentile of the Y-axis linear velocity (signed).
     */
    ctx->features[SMARTPEN_F_Y_VELOCITY_99PC] =
        percentile(ctx->vy, n_samples, 99.0, g_percentile_scratch);

    /* ── Feature 7  (CSV #432) ─ last_azimuth_skewness ──
     * Skewness of the azimuth angle computed over the FINAL 20% temporal
     * epoch of the trial. last_* features conventionally slice the tail.
     */
    {
        const int start = static_cast<int>(std::floor(0.80 * n_samples));
        const int len   = n_samples - start;
        const Moments m = compute_moments(&ctx->azimuth[start], len);
        ctx->features[SMARTPEN_F_LAST_AZIMUTH_SKEWNESS] = m.skewness;
    }

    /* ── Feature 8  (CSV #629) ─ stroke_pressure_1st_each_stroke_median ──
     * Median of the dP/dt pool across all segmented strokes.
     */
    ctx->features[SMARTPEN_F_STROKE_PRESSURE_1ST_MEDIAN] =
        median(ctx->dPdt_pool, ctx->dPdt_count, g_percentile_scratch);

    /* ── Feature 9  (CSV #556) ─ angle_data_with_1_mean ──
     * Mean of the angle_data_with_1 signal (signed tilt of the pen axis
     * relative to its initial orientation) over the WHOLE trial.
     */
    ctx->features[SMARTPEN_F_ANGLE_DATA_WITH_1_MEAN] =
        compute_mean(ctx->angle_with_1, n_samples);

    /* ── Feature 10 (CSV #111) ─ y_velocity_displacement_between_99th_percentile_and_1st_percentile ──
     * Inter-percentile range: P99(vy) − P1(vy). Note the brief calls this a
     * "displacement" but the convention in 1_1.csv is the *difference*
     * between the two percentiles (a robust range, not an integral).
     */
    {
        const double p99 = percentile(ctx->vy, n_samples, 99.0, g_percentile_scratch);
        const double p01 = percentile(ctx->vy, n_samples,  1.0, g_percentile_scratch);
        ctx->features[SMARTPEN_F_Y_VELOCITY_99_1_RANGE] = p99 - p01;
    }

    /* ── Feature 11 (CSV #698) ─ stroke_dis_median_each_stroke_max ──
     * For each stroke we computed stroke_dis_median[s] (median of the
     * resultant displacement magnitude within that stroke). The feature is
     * the MAXIMUM of those per-stroke medians.
     */
    {
        double dmax = 0.0;
        for (int s = 0; s < ctx->n_strokes; ++s)
            if (s == 0 || ctx->stroke_dis_median[s] > dmax)
                dmax = ctx->stroke_dis_median[s];
        ctx->features[SMARTPEN_F_STROKE_DIS_MEDIAN_EACH_STROKE_MAX] = dmax;
    }

    /* ── Feature 12 (CSV #236) ─ first_x_median ──
     * Median of the X-axis displacement computed over the INITIAL 20%
     * temporal epoch of the trial.
     */
    {
        const int end = static_cast<int>(std::floor(0.20 * n_samples));
        const int len = (end > 0) ? end : 1;
        ctx->features[SMARTPEN_F_FIRST_X_MEDIAN] =
            median(ctx->x_disp, len, g_percentile_scratch);
    }

    /* ── Feature 13 (CSV #86) ─ x_displacement_kurtosis ──
     * Excess kurtosis of the X-axis linear displacement distribution
     * over the WHOLE trial.
     */
    {
        const Moments m = compute_moments(ctx->x_disp, n_samples);
        ctx->features[SMARTPEN_F_X_DISPLACEMENT_KURTOSIS] = m.kurtosis;
    }

    /* ── Feature 14 (CSV #547) ─ last_y_velocity_no_abs_median ──
     * Median of the RAW (non-absolute) Y-axis velocity over the FINAL 20%
     * temporal epoch of the trial. "no_abs" is critical: we must NOT take
     * |vy| before computing the median.
     */
    {
        const int start = static_cast<int>(std::floor(0.80 * n_samples));
        const int len   = n_samples - start;
        ctx->features[SMARTPEN_F_LAST_Y_VELOCITY_NO_ABS_MEDIAN] =
            median(&ctx->vy[start], len, g_percentile_scratch);
    }

    return &ctx->features[0];
}


int smartpen_feature_count(void) noexcept
{
    return SMARTPEN_FEATURE_COUNT;
}


const char* smartpen_feature_name(int index) noexcept
{
    if (index < 0 || index >= SMARTPEN_FEATURE_COUNT) return "unknown";
    return g_feature_meta[index].name;
}


int smartpen_feature_csv_column(int index) noexcept
{
    if (index < 0 || index >= SMARTPEN_FEATURE_COUNT) return -1;
    return g_feature_meta[index].csv_column;
}

// في نهاية ملف smartPen_Features.cpp داخل بلوك extern "C"
void smartpen_get_ui_data(
    const SmartPenFeatures* ctx,
    const float** pressure_plot,
    const float** acc_x_plot,
    const float** acc_y_plot,
    const float** acc_z_plot,
    const double** tremor_velocity,
    const double** pressure_stability,
    int* total_samples,
    int* stability_samples
) noexcept
{
    if (ctx == nullptr) return;

    // تمرير عناوين فضاء الذاكرة للمصفوفات الزمنية مباشرة إلى Dart FFI
    if (pressure_plot)      *pressure_plot     = &ctx->pressure[0];
    if (acc_x_plot)         *acc_x_plot        = &ctx->acc_x[0];
    if (acc_y_plot)         *acc_y_plot        = &ctx->acc_y[0];
    if (acc_z_plot)         *acc_z_plot        = &ctx->acc_z[0];
    if (tremor_velocity)    *tremor_velocity   = &ctx->velocity[0];
    if (pressure_stability)  *pressure_stability = &ctx->dPdt_pool[0];

    // تمرير أطوال المصفوفات الفعلي المحسوب في هذا الفحص
    if (total_samples)      *total_samples     = ctx->n_samples;
    if (stability_samples)  *stability_samples = ctx->dPdt_count;
}

} /* extern "C" */
