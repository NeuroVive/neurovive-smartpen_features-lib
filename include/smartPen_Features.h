/**
 * @file smartPen_Features.h
 * @brief NeuroVive Smart Pen — Real-Time Feature Extraction Library
 *        C-compatible clean interface for Dart FFI consumption.
 *
 * PROJECT       : NeuroVive Smart Pen — Task 1: Archimedean Spiral Screening
 * AUTHOR        : Elite Embedded Systems & DSP Engineer
 * LANGUAGE      : Modern C++17 (MISRA C++ safety-aligned)
 * HARDWARE      : MPU6050 (6-DOF IMU) + FSR 402 (Force-Sensitive Resistor)
 * SAMPLE RATE   : 150 Hz synchronous streaming (BLE 5.0 telemetry)
 *
 * ARCHITECTURE CONSTRAINTS ENFORCED:
 *   - Strictly STATIC memory allocation. Zero runtime heap usage.
 *   - Deterministic execution: every trial analysis < 10 ms over Dart FFI.
 *   - All buffers live in the SmartPenFeatures struct so the caller owns the
 *     memory and there is no hidden global state to collide across isolates.
 *   - No exceptions cross the FFI boundary; every entry point is noexcept.
 *
 * FEATURE VECTOR CONTRACT (15 elements, fixed order):
 *   The returned pointer addresses a 15-element double array whose indices
 *   MUST stay aligned with the downstream Python inference engine. The exact
 *   mathematical order is documented beside SmartPenFeatureIndex below and is
 *   mirrored inside smartPen_Features.cpp.
 *
 * REFERENCE DATASET NOMENCLATURE VERIFIED AGAINST:
 *   1_1.csv — PaHaW-derived feature schema (870 named columns). The 15 golden
 *   feature names below were matched character-for-character to the CSV header
 *   so that the C++ output drops directly into the production inference frame.
 */
#ifndef SMARTPEN_FEATURES_H
#define SMARTPEN_FEATURES_H

#ifdef __cplusplus
extern "C" {
#endif

/* ───────────────────────────────────────────────────────────────────────────
 * 1. Compile-time hardware & algorithm constants
 * ─────────────────────────────────────────────────────────────────────────── */

/** Hardware clock. BLE 5.0 packet cadence of 7.5 ms ⇒ 1 / 0.0075 = 133.33 Hz
 *  nominal, but the firmware streams one fused sample per 6.666 ms slot giving
 *  exactly 150 Hz synchronous across IMU + FSR. This constant is the single
 *  source of truth used by every integrator and Butterworth design below. */
#define SMARTPEN_SAMPLE_RATE_HZ        150.0

/** Inverse of the sample rate, pre-computed in double precision so the
 *  integrators never perform a division per sample (MISRA-aligned constant
 *  propagation). */
#define SMARTPEN_DT_SEC                (1.0 / SMARTPEN_SAMPLE_RATE_HZ)

/** FSR 402 pen-down threshold, per NeuroVive hardware spec
 *  (شرح القلم.md, line 16-17). Pressure normalised to [0,1] by the firmware;
 *  any sample strictly greater than 0.05 is treated as on-surface writing. */
#define SMARTPEN_PRESSURE_THRESHOLD    0.05

/** ZUPT gate. The complementary-filter integrator zeroes linear velocity
 *  whenever pressure drops to/below this value, eliminating quadratic drift
 *  accumulation during in-air transitions. */
#define SMARTPEN_ZUPT_GATE             SMARTPEN_PRESSURE_THRESHOLD

/** Complementary filter gain (α). The gyroscope high-pass branch is weighted
 *  by α = 0.98 (drift-stable over short horizons), the accelerometer
 *  low-pass branch by (1 - α) = 0.02 (noise rejective). Reference value per
 *  شرح القلم.md, line 59. */
#define SMARTPEN_COMPLEMENTARY_ALPHA   0.98

/** High-pass Butterworth cut-off used to suppress DC + gravity residual
 *  before the double integration of acceleration. 0.5 Hz was selected so the
 *  filter removes integration drift without eroding the spiral writing band
 *  (~0.3–3 Hz hand tremor + intentional motion). */
#define SMARTPEN_HPF_CUTOFF_HZ         0.5

/** Single-pole IIR intrinsic-mode filter coefficient. The single-pass digital
 *  filter isolates the primary dynamic mode of the velocity vector for the
 *  intrinsic_feature_teager_kaiser_energy(velocity, 1) feature. */
#define SMARTPEN_INTRINSIC_ALPHA       0.85

/** Static buffer caps. Sized to comfortably hold a 15-second Archimedean
 *  spiral trial at 150 Hz (2250 samples) plus head-room for the in-air
 *  transitions. Lifted to 4096 for power-of-two alignment and to bound the
 *  worst-case QuickSelect / TKE loops deterministically. */
#define SMARTPEN_MAX_SAMPLES           4096

/** Maximum number of strokes the segmenter can track without dynamic
 *  allocation. A typical spiral trial breaks into 3–8 on-surface strokes;
 *  64 provides ample head-room while keeping the per-stroke statistics array
 *  in L1 cache. */
#define SMARTPEN_MAX_STROKES           64

/** Feature vector dimensionality — the downstream Python model expects
 *  exactly 15 elements. */
#define SMARTPEN_FEATURE_COUNT         15


/* ───────────────────────────────────────────────────────────────────────────
 * 2. Feature index enumeration — strict alignment with PaHaW CSV schema
 *    (verified against 1_1.csv header rows 185, 187, 207-231, etc.)
 * ─────────────────────────────────────────────────────────────────────────── */
enum SmartPenFeatureIndex
{
    /* 0  */ SMARTPEN_F_TEAGER_KAISER_DISPLACEMENT          = 0,  /* CSV#185  teager_kaiser_energy(displacement, 1)            */
    /* 1  */ SMARTPEN_F_STROKE_PRESSURE_1ST_MEAN            = 1,  /* CSV#628  stroke_pressure_1st_each_stroke_mean             */
    /* 2  */ SMARTPEN_F_TEAGER_KAISER_VELOCITY              = 2,  /* CSV#187  teager_kaiser_energy(velocity, 1)                */
    /* 3  */ SMARTPEN_F_STROKE_PRESSURE_1ST_MIN             = 3,  /* CSV#627  stroke_pressure_1st_each_stroke_min              */
    /* 4  */ SMARTPEN_F_STROKE_PRESSURE_1ST_1PC             = 4,  /* CSV#632  stroke_pressure_1st_each_stroke_1st_percentile   */
    /* 5  */ SMARTPEN_F_INTRINSIC_TKE_VELOCITY              = 5,  /* CSV#231  intrinsic_feature_teager_kaiser_energy(velocity) */
    /* 6  */ SMARTPEN_F_Y_VELOCITY_99PC                     = 6,  /* CSV#110  y_velocity_99th_percentile                       */
    /* 7  */ SMARTPEN_F_LAST_AZIMUTH_SKEWNESS               = 7,  /* CSV#432  last_azimuth_skewness                            */
    /* 8  */ SMARTPEN_F_STROKE_PRESSURE_1ST_MEDIAN          = 8,  /* CSV#629  stroke_pressure_1st_each_stroke_median           */
    /* 9  */ SMARTPEN_F_ANGLE_DATA_WITH_1_MEAN              = 9,  /* CSV#556  angle_data_with_1_mean                           */
    /* 10 */ SMARTPEN_F_Y_VELOCITY_99_1_RANGE               = 10, /* CSV#111  y_velocity_displacement_between_99th_percentile  */
    /* 11 */ SMARTPEN_F_STROKE_DIS_MEDIAN_EACH_STROKE_MAX   = 11, /* CSV#698  stroke_dis_median_each_stroke_max                */
    /* 12 */ SMARTPEN_F_FIRST_X_MEDIAN                      = 12, /* CSV#236  first_x_median                                    */
    /* 13 */ SMARTPEN_F_X_DISPLACEMENT_KURTOSIS             = 13, /* CSV#86   x_displacement_kurtosis                          */
    /* 14 */ SMARTPEN_F_LAST_Y_VELOCITY_NO_ABS_MEDIAN       = 14, /* CSV#547  last_y_velocity_no_abs_median                    */
    SMARTPEN_F_COUNT                                          = SMARTPEN_FEATURE_COUNT
};


/* ───────────────────────────────────────────────────────────────────────────
 * 3. Engine state — caller-owned, single static allocation
 *
 * Every working buffer the DSP pipeline needs is bundled here so the Dart
 * FFI caller can `malloc` exactly once (or even place the struct on the
 * native heap via `Pointer<SmartPenFeatures>` and free it on dispose).
 * No `new`, no `std::vector::push_back`, no hidden global — satisfies the
 * zero-heap-fragmentation firmware constraint.
 * ─────────────────────────────────────────────────────────────────────────── */
typedef struct SmartPenFeatures
{
    /* ---- Raw sensor input mirror (n samples used) ---- */
    int      n_samples;

    float    acc_x[SMARTPEN_MAX_SAMPLES];
    float    acc_y[SMARTPEN_MAX_SAMPLES];
    float    acc_z[SMARTPEN_MAX_SAMPLES];
    float    gyro_x[SMARTPEN_MAX_SAMPLES];
    float    gyro_y[SMARTPEN_MAX_SAMPLES];
    float    gyro_z[SMARTPEN_MAX_SAMPLES];
    float    pressure[SMARTPEN_MAX_SAMPLES];

    /* ---- Derived kinematic signals ---- */
    double   vx[SMARTPEN_MAX_SAMPLES];          /* X-axis linear velocity (m/s)             */
    double   vy[SMARTPEN_MAX_SAMPLES];          /* Y-axis linear velocity (m/s)             */
    double   x_disp[SMARTPEN_MAX_SAMPLES];      /* X-axis linear displacement (m)           */
    double   y_disp[SMARTPEN_MAX_SAMPLES];      /* Y-axis linear displacement (m)           */
    double   velocity[SMARTPEN_MAX_SAMPLES];    /* Resultant planar speed |v| = sqrt(vx²+vy²) */
    double   displacement[SMARTPEN_MAX_SAMPLES];/* Resultant planar displacement magnitude   */

    /* ---- Pen attitude (complementary filter outputs, degrees) ---- */
    double   azimuth[SMARTPEN_MAX_SAMPLES];     /* Yaw/rotation about vertical              */
    double   altitude[SMARTPEN_MAX_SAMPLES];    /* Tilt/elevation above horizontal          */
    double   angle_with_1[SMARTPEN_MAX_SAMPLES];/* Tilt relative to primary init axis (deg) */

    /* ---- Intrinsic (primary dynamic mode) velocity signal ---- */
    double   intrinsic_velocity[SMARTPEN_MAX_SAMPLES];

    /* ---- Stroke segmentation map ---- */
    int      n_strokes;
    int      stroke_start[SMARTPEN_MAX_STROKES]; /* inclusive sample index                   */
    int      stroke_end  [SMARTPEN_MAX_STROKES]; /* exclusive sample index                   */

    /* ---- Per-stroke median displacement (for feature 11) ---- */
    double   stroke_dis_median[SMARTPEN_MAX_STROKES];

    /* ---- Per-stroke first-derivative of pressure statistics (features 1,3,4,8)
     *      The pipeline accumulates the global pool of dP/dt samples taken
     *      from inside strokes into dPdt_pool, and dPdt_count tracks its size. */
    double   dPdt_pool[SMARTPEN_MAX_SAMPLES];
    int      dPdt_count;

    /* ---- Output: the canonical 15-element feature vector ----
     *      Addressable as `&features[0]` for Dart FFI `Pointer<Double>`. */
    double   features[SMARTPEN_FEATURE_COUNT];

} SmartPenFeatures;


/* ───────────────────────────────────────────────────────────────────────────
 * 4. Public C API
 *
 * Every function is `extern "C"`, noexcept, and reentrant provided each
 * caller passes its own SmartPenFeatures instance. None of them allocate.
 * ─────────────────────────────────────────────────────────────────────────── */

/**
 * Zero-initialise the engine state. Must be called once before each trial.
 * Equivalent to memset(0) but explicit and self-documenting.
 */
void smartpen_init(SmartPenFeatures* ctx) noexcept;

/**
 * Run the full real-time DSP pipeline over a synchronous 150 Hz trial buffer.
 *
 * @param ctx       Engine state (caller-owned, statically allocated).
 * @param acc_x     Pointer to n linear acceleration samples along X (g).
 * @param acc_y     Pointer to n linear acceleration samples along Y (g).
 * @param acc_z     Pointer to n linear acceleration samples along Z (g).
 * @param gyro_x    Pointer to n angular velocity samples around X (deg/s).
 * @param gyro_y    Pointer to n angular velocity samples around Y (deg/s).
 * @param gyro_z    Pointer to n angular velocity samples around Z (deg/s).
 * @param pressure  Pointer to n normalised FSR pressure samples in [0,1].
 * @param n_samples Number of samples in each input array (<= SMARTPEN_MAX_SAMPLES).
 *
 * @return Pointer to ctx->features, a 15-element double array. The pointer
 *         remains valid until the next call to smartpen_process() or
 *         smartpen_init() on the same context.
 */
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
) noexcept;

/**
 * @return The canonical feature-vector length (always 15).
 */
int smartpen_feature_count(void) noexcept;

/**
 * Lookup the canonical feature name (PaHaW CSV nomenclature) by index.
 * @param index 0..14 inclusive. Out-of-range returns "unknown".
 */
const char* smartpen_feature_name(int index) noexcept;

/**
 * Lookup the CSV column index that the given feature maps to. Useful for
 * downstream debug printing / parity tests against the Python frame.
 * @param index 0..14 inclusive. Out-of-range returns -1.
 */
int smartpen_feature_csv_column(int index) noexcept;




// داخل extern "C" في ملف smartPen_Features (1).h
void smartpen_get_ui_data(
    const SmartPenFeatures* ctx,
    const float** pressure_plot,       // للرسمة الأولى: Pressure (kpa)
    const float** acc_x_plot,          // للرسمة الثانية: Acceleration (g)
    const float** acc_y_plot,
    const float** acc_z_plot,
    const double** tremor_velocity,    // للرسمة الثالثة: Motion Tremor (Hz)
    const double** pressure_stability, // للرسمة الرابعة: Pressure Stability
    int* total_samples,
    int* stability_samples
) noexcept;


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SMARTPEN_FEATURES_H */
