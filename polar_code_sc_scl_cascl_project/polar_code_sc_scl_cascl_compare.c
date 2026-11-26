#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200809L
#endif
/*
 * polar_code_sc_scl_cascl_compare.c
 *
 * SC vs SCL(L=2,4,8) vs CRC-aided SCL(L=2,4,8) BER/FER comparison.
 *
 * Runtime code parameters supplied by the GUI:
 *   N is a power of two from 64 to 1024; R selects K=round(N*R) non-frozen positions.
 *   CA-SCL reserves 8 of the K non-frozen slots for CRC-8; payload K_INFO=K-8.
 *   Arrays are allocated for the maximum supported N, while the algorithms use runtime N/K.
 *
 * All three decoder families (SC, SCL, CA-SCL) are implemented with ONE
 * unified recursive list-decoder:
 *   - SC    == list decoding with list size L = 1 (the fork-then-keep-best-1
 *             rule is exactly equivalent to a normal hard SC decision).
 *   - SCL   == list decoding with L in {2,4,8}, pick lowest path metric.
 *   - CA-SCL== same list decoding with L in {2,4,8}, but among the L final
 *             candidates we pick the lowest-path-metric one that PASSES the
 *             CRC check; if none pass, we fall back to the global PM-minimum
 *             (the usual CA-SCL fallback rule).
 *
 * Frozen-mask construction uses the standard Gaussian-Approximation (GA)
 * recursion (phi/phi_inv), matched to each Eb/No sweep point (same
 * convention as the project's earlier scd/scl "Baseline" simulations).
 *
 * Build:
gcc -O2 -Wall -static -o compare.exe polar_code_sc_scl_cascl_compare.c -lm
 * Run:
./compare.exe

 *
 * NOTE ON RUNTIME: this does 7 Eb/No points... no, NUM_SNR_POINTS x 7
 * decoder configs Monte Carlo runs, and CA-SCL(L=8) is the most expensive.
 *
 * [수정] MAX_FRAMES 상한을 완전히 없앴습니다. 각 (Eb/No, 디코더) 조합은
 * 이제 TARGET_ERRORS(기본 5000)개의 비트 에러를 실제로 관측할 때까지
 * 프레임 수 제한 없이 계속 돕니다. 예전엔 MAX_FRAMES에 먼저 걸려서
 * 에러를 하나도 못 보고 BER/FER이 그대로 0.00e+00으로 찍히는 포인트가
 * 생겼는데(특히 Eb/No가 높고 CA-SCL(L=8)처럼 에러율이 낮은 조합), 이제는
 * 시간이 얼마가 걸리든 무조건 5000개를 채울 때까지 기다립니다. 대신 한
 * 포인트가 몇 시간씩 걸릴 수 있으므로, 50000프레임마다(그리고 최소 15초
 * 간격으로) 진행 상황(프레임 수/현재까지 에러 수/경과 시간)을 콘솔에
 * 찍어서 멈춘 게 아니라 계속 도는 중이라는 걸 확인할 수 있게 했습니다.
 * 결과는 (기존과 동일하게) 포인트 하나가 끝날 때마다 바로 파일에 저장/
 * flush되므로, 중간에 중단해도 그때까지 끝난 포인트는 남아있습니다.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#if defined(_WIN32)
  #include <windows.h>
  #include <direct.h>
  #define MKDIR(d) _mkdir(d)
#else
  #include <sys/stat.h>
  #include <sys/types.h>
  #define MKDIR(d) mkdir((d), 0755)
#endif

/* ------------------------------------------------------------------ */
/* Code parameters                                                     */
/* ------------------------------------------------------------------ */
#define N          1024        /* array capacity; runtime block length is g_n */
#define K          1024        /* array capacity; runtime info length is g_k */
#define N_STAGES   10          /* maximum log2(N) */
#define CRC_LEN    8
#define CRC_POLY   0x07        /* CRC-8, x^8 + x^2 + x + 1              */
#define K_INFO     (g_k - CRC_LEN)
#define MAX_LIST   8

/* Eb/No sweep */
#define EBNO_START      -1.0
#define EBNO_STEP       0.25
#define NUM_SNR_POINTS  29      /* 0.0 .. 3.5 dB */

/*
 * Monte Carlo stopping rule.
 *
 * MAX_FRAMES 상한은 없앴습니다 - 아래 run_point()의 while 루프는
 * bit_errors가 TARGET_ERRORS에 도달할 때까지 프레임 수 제한 없이 계속
 * 돕니다. CA-SCL(L=8)처럼 Eb/No가 높아질수록 에러가 거의 안 나는 조합은
 * 한 포인트에 몇 시간씩 걸릴 수 있지만(사용자 확인: 몇 시간 걸려도
 * 무방), 그 대신 모든 포인트가 반드시 TARGET_ERRORS개의 실제 에러
 * 관측치를 갖게 되어 BER/FER이 0으로 끊기거나 그래프가 요동치는 문제가
 * 없어집니다.
 */
#define TARGET_ERRORS   5000

/* Runtime settings supplied by the native GUI. */
#define MAX_RUNTIME_SNR_POINTS 201
static double g_ebno_start = EBNO_START;
static double g_ebno_step = EBNO_STEP;
static int g_num_snr_points = NUM_SNR_POINTS;
static long g_target_errors = TARGET_ERRORS;
static int g_n = 1024;
static int g_k = 512;
static int g_n_stages = 10;
static char g_output_dir[4096] = "result";

static bool valid_code_parameters(int n, double r, int *k_out, int *depth_out) {
    if (n < 64 || n > 1024 || (n & (n - 1)) != 0 || !isfinite(r) || r < 0.15 || r > 0.90)
        return false;
    int k = (int)lround((double)n * r);
    if (k <= CRC_LEN || k > n) return false;
    int depth = 0;
    for (int v = n; v > 1; v >>= 1) ++depth;
    if (k_out) *k_out = k;
    if (depth_out) *depth_out = depth;
    return true;
}

/* 진행 상황 출력 주기: 최소 이만큼 프레임이 지날 때마다 "최소 PRINT_MIN_SECONDS초
   간격"으로 진행 상황을 찍는다 (한 포인트가 오래 걸려도 멈춘 게 아니라는 걸
   확인할 수 있도록). */
#define PRINT_EVERY_FRAMES   50000
#define PRINT_MIN_SECONDS    15.0

static const int LIST_SIZES[3] = {2, 4, 8};
#define NUM_LIST_SIZES 3
#define NUM_CONFIGS (1 + NUM_LIST_SIZES * 2)   /* SC + 3xSCL + 3xCASCL = 7 */

typedef enum { MODE_SC, MODE_SCL, MODE_CASCL } DecoderMode;

typedef struct {
    char name[32];
    DecoderMode mode;
    int l_limit;
} DecoderConfig;

typedef struct {
    double ber;
    double fer;
    long frames;
    long bit_errors;
} RunResult;

/* ------------------------------------------------------------------ */
/* RNG: Xorshift64 (consistent with the rest of the project)           */
/* ------------------------------------------------------------------ */
static uint64_t rng_state;

static uint64_t xorshift64(void) {
    uint64_t x = rng_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    rng_state = x;
    return x;
}

static double rand_double(void) {
    /* 53-bit resolution double in (0,1) */
    return (double)(xorshift64() >> 11) * (1.0 / 9007199254740992.0);
}

static void rand_gaussian(double *z0, double *z1) {
    double u1, u2;
    do { u1 = rand_double(); } while (u1 <= 1e-12);
    u2 = rand_double();
    double r = sqrt(-2.0 * log(u1));
    *z0 = r * cos(2.0 * M_PI * u2);
    *z1 = r * sin(2.0 * M_PI * u2);
}

static inline uint8_t rand_bit(void) {
    return (uint8_t)(xorshift64() & 1ULL);
}

/* ------------------------------------------------------------------ */
/* Gaussian Approximation (GA) phi / phi_inv                           */
/* ------------------------------------------------------------------ */
static double phi(double x) {
    if (x <= 0.0) return 1.0;
    double v;
    if (x > 10.0) {
        v = sqrt(M_PI / x) * exp(-x / 4.0) * (1.0 - 10.0 / (7.0 * x));
    } else {
        v = exp(-0.4527 * pow(x, 0.86) + 0.0218);
    }
    if (v > 1.0) v = 1.0;
    if (v < 0.0) v = 0.0;
    return v;
}

static double phi_inv(double y) {
    if (y >= 1.0) return 0.0;
    if (y <= 0.0) return 200.0;
    double lo = 0.0, hi = 200.0;
    for (int it = 0; it < 60; it++) {
        double mid = 0.5 * (lo + hi);
        double v = phi(mid);
        if (v > y) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

typedef struct { int idx; double val; } IdxVal;

static int cmp_idxval_desc(const void *a, const void *b) {
    const IdxVal *A = (const IdxVal *)a, *B = (const IdxVal *)b;
    if (A->val > B->val) return -1;
    if (A->val < B->val) return 1;
    return 0;
}

/*
 * Builds frozen[N] (1 = frozen, 0 = info/non-frozen) using GA, designed for
 * the given Eb/No (dB) assuming runtime code rate R_design=K/N; CA-SCL
 * still uses this same K-position mask -- its extra CRC
 * overhead is handled later by splitting the K non-frozen slots, not by
 * changing which slots are non-frozen).
 *
 * ga_recursive() mirrors scl_decode_recursive()'s split EXACTLY (top-down:
 * left half = f/minus, right half = g/plus, then recurse into each half
 * separately) rather than a bottom-up/iterative butterfly sweep. Those two
 * orderings are NOT interchangeable for a recursive (non-bit-reversed)
 * encoder/decoder like this one -- an iterative bottom-up GA sweep silently
 * produces a mask barely better than random (confirmed during development
 * via a genie-aided per-position error-rate measurement: ~71% overlap with
 * the true best-K positions and ~400x worse average BER). Mirroring the
 * decoder's own recursion gives >99% overlap and matches the true optimum.
 */
static void ga_recursive(double *x, int n2) {
    if (n2 == 1) return;
    int half = n2 / 2;
    double xm[512], xp[512];
    for (int i = 0; i < half; i++) {
        double a = x[i], b = x[i + half];
        double pa = phi(a), pb = phi(b);
        double arg = 1.0 - (1.0 - pa) * (1.0 - pb);
        if (arg < 1e-300) arg = 1e-300;
        if (arg > 1.0) arg = 1.0;
        xm[i] = phi_inv(arg);
        xp[i] = a + b;
    }
    for (int i = 0; i < half; i++) { x[i] = xm[i]; x[i + half] = xp[i]; }
    ga_recursive(x, half);
    ga_recursive(x + half, half);
}

static void build_frozen_mask(double ebno_db, double r_design, int *frozen) {
    double ebno_lin = pow(10.0, ebno_db / 10.0);
    double sigma2 = 1.0 / (2.0 * r_design * ebno_lin);

    double x[N];
    for (int i = 0; i < g_n; i++) x[i] = 2.0 / sigma2;

    ga_recursive(x, g_n);

    IdxVal iv[N];
    for (int i = 0; i < g_n; i++) { iv[i].idx = i; iv[i].val = x[i]; }
    qsort(iv, (size_t)g_n, sizeof(IdxVal), cmp_idxval_desc);

    for (int i = 0; i < g_n; i++) frozen[i] = 1;
    for (int i = 0; i < g_k; i++) frozen[iv[i].idx] = 0;
}

/* ------------------------------------------------------------------ */
/* Polar encoder (recursive, Arikan butterfly)                         */
/* ------------------------------------------------------------------ */
static void polar_encode(const uint8_t *u, uint8_t *xout, int n2) {
    if (n2 == 1) { xout[0] = u[0]; return; }
    int half = n2 / 2;
    uint8_t u1[512] = {0}, u2[512] = {0};
    uint8_t x1[512], x2[512];
    for (int i = 0; i < half; i++) {
        u1[i] = u[i] ^ u[i + half];
        u2[i] = u[i + half];
    }
    polar_encode(u1, x1, half);
    polar_encode(u2, x2, half);
    for (int i = 0; i < half; i++) {
        xout[i] = x1[i];
        xout[i + half] = x2[i];
    }
}

/* ------------------------------------------------------------------ */
/* CRC-8 (bit-serial LFSR, poly 0x07) -- used identically for          */
/* generation (encoder side) and verification (CA-SCL decoder side)    */
/* ------------------------------------------------------------------ */
static uint8_t crc8_compute(const uint8_t *bits, int len) {
    uint8_t crc = 0x00;
    for (int i = 0; i < len; i++) {
        uint8_t bit_in = bits[i] & 1;
        uint8_t msb = (crc >> 7) & 1;
        crc = (uint8_t)(crc << 1);
        if (bit_in ^ msb) crc ^= CRC_POLY;
    }
    return crc;
}

/* ------------------------------------------------------------------ */
/* Unified SC / SCL / CA-SCL recursive list decoder                     */
/* ------------------------------------------------------------------ */
static inline double softplus(double v) {
    /* log(1+exp(v)), numerically stable */
    if (v > 0.0) return v + log1p(exp(-v));
    return log1p(exp(v));
}

static inline double f_combine_minsum(double l1, double l2) {
    double sgn = ((l1 < 0.0) != (l2 < 0.0)) ? -1.0 : 1.0;
    double minabs = fmin(fabs(l1), fabs(l2));
    return sgn * minabs;
}

/*
 * Recursively decodes a block of size n2 starting at absolute index
 * base_idx, operating on up to MAX_LIST parallel paths at once.
 *
 * llr[p][0..n2-1]   : this path's combined LLRs for this subtree (indices
 *                      are LOCAL to the subtree, i.e. 0..n2-1)
 * ubits[p][0..N-1]  : this path's full decided-bit history (GLOBAL array,
 *                      absolute indices); only base_idx..base_idx+n2-1 get
 *                      written by this call (directly at the leaves)
 * pm[p]             : this path's accumulated path metric (lower = better)
 * path_count        : in/out -- number of active paths
 * origin_out[newp]  : out -- which ENTRY path index (0..path_count_in-1)
 *                      the surviving/duplicated path newp descends from.
 *                      Needed by the caller to correctly re-index its own
 *                      (already-computed) llr array when building the
 *                      g-function inputs for the right subtree.
 */
static void scl_decode_recursive(double llr[][N], uint8_t ubits[][N], int n2,
                                  int base_idx, const int *frozen,
                                  double pm[MAX_LIST], int *path_count,
                                  int l_limit, int origin_out[MAX_LIST]) {
    int pc = *path_count;

    if (n2 == 1) {
        int idx = base_idx;
        if (frozen[idx]) {
            for (int p = 0; p < pc; p++) {
                ubits[p][idx] = 0;
                pm[p] += softplus(-llr[p][0]);
                origin_out[p] = p;
            }
            return; /* path_count unchanged */
        }

        int cand_n = pc * 2;
        double cand_pm[2 * MAX_LIST];
        int cand_parent[2 * MAX_LIST];
        int cand_bit[2 * MAX_LIST];
        for (int p = 0; p < pc; p++) {
            cand_parent[2 * p] = p; cand_bit[2 * p] = 0;
            cand_pm[2 * p] = pm[p] + softplus(-llr[p][0]);
            cand_parent[2 * p + 1] = p; cand_bit[2 * p + 1] = 1;
            cand_pm[2 * p + 1] = pm[p] + softplus(llr[p][0]);
        }

        int order[2 * MAX_LIST];
        for (int i = 0; i < cand_n; i++) order[i] = i;
        int keep = (cand_n < l_limit) ? cand_n : l_limit;
        for (int a = 0; a < keep; a++) {
            int best = a;
            for (int b = a + 1; b < cand_n; b++)
                if (cand_pm[order[b]] < cand_pm[order[best]]) best = b;
            int t = order[a]; order[a] = order[best]; order[best] = t;
        }

        double pm_new[MAX_LIST];
        static uint8_t ubits_new[MAX_LIST][N];
        for (int a = 0; a < keep; a++) {
            int c = order[a];
            int par = cand_parent[c];
            pm_new[a] = cand_pm[c];
            memcpy(ubits_new[a], ubits[par], (size_t)g_n);
            ubits_new[a][idx] = (uint8_t)cand_bit[c];
            origin_out[a] = par;
        }
        for (int a = 0; a < keep; a++) {
            pm[a] = pm_new[a];
            memcpy(ubits[a], ubits_new[a], (size_t)g_n);
        }
        *path_count = keep;
        return;
    }

    int half = n2 / 2;
    double llr_minus[MAX_LIST][N];
    for (int p = 0; p < pc; p++)
        for (int i = 0; i < half; i++)
            llr_minus[p][i] = f_combine_minsum(llr[p][i], llr[p][i + half]);

    int origin_left[MAX_LIST];
    scl_decode_recursive(llr_minus, ubits, half, base_idx, frozen, pm,
                          path_count, l_limit, origin_left);
    int pc_after_left = *path_count;

    /*
     * g-function needs the RE-ENCODED combination of the left branch's
     * decoded bits (the "partial sum"), not the raw decoded bits
     * themselves -- they only coincide at the leaf-adjacent level
     * (half == 1). Using the raw bits here silently corrupts every
     * decode where half >= 2 (verified by an exhaustive N=4 round-trip
     * test), so always re-encode before building llr_plus.
     */
    double llr_plus[MAX_LIST][N];
    uint8_t xleft_hat[N];
    for (int p = 0; p < pc_after_left; p++) {
        int par = origin_left[p];
        polar_encode(&ubits[p][base_idx], xleft_hat, half);
        for (int i = 0; i < half; i++) {
            double l1 = llr[par][i];
            double l2 = llr[par][i + half];
            uint8_t ub = xleft_hat[i];
            llr_plus[p][i] = l2 + (ub ? -l1 : l1);
        }
    }

    int origin_right[MAX_LIST];
    scl_decode_recursive(llr_plus, ubits, half, base_idx + half, frozen, pm,
                          path_count, l_limit, origin_right);
    int pc_after_right = *path_count;

    for (int p = 0; p < pc_after_right; p++)
        origin_out[p] = origin_left[origin_right[p]];
}

/* ------------------------------------------------------------------ */
/* Monte Carlo core (shared by SC / SCL / CA-SCL)                      */
/* ------------------------------------------------------------------ */
static double g_llr[MAX_LIST][N];
static uint8_t g_ubits[MAX_LIST][N];

static RunResult run_point(DecoderConfig cfg, double ebno_db, const int *frozen,
                            uint64_t seed) {
    double r_eff = (cfg.mode == MODE_CASCL) ? ((double)K_INFO / (double)g_n)
                                             : ((double)g_k / (double)g_n);
    double sigma2 = 1.0 / (2.0 * r_eff * pow(10.0, ebno_db / 10.0));
    double sigma = sqrt(sigma2);

    long frames = 0;
    long bit_errors = 0;
    long frame_errors = 0;
    long score_len = (cfg.mode == MODE_CASCL) ? K_INFO : g_k;

    rng_state = seed;
    if (rng_state == 0) rng_state = 1;

    uint8_t sentInfo[K];
    uint8_t u[N];
    uint8_t x[N];
    double rx[N];

    time_t t_run_start = time(NULL);
    time_t t_last_print = t_run_start;

    /* [수정] MAX_FRAMES 상한 없이, TARGET_ERRORS개의 비트 에러를 실제로
       볼 때까지 무조건 돈다. */
    while (bit_errors < g_target_errors) {
        frames++;

        if (cfg.mode == MODE_CASCL) {
            for (int i = 0; i < K_INFO; i++) sentInfo[i] = rand_bit();
            uint8_t crc = crc8_compute(sentInfo, K_INFO);
            for (int b = 0; b < CRC_LEN; b++)
                sentInfo[K_INFO + b] = (uint8_t)((crc >> (CRC_LEN - 1 - b)) & 1);
        } else {
            for (int i = 0; i < g_k; i++) sentInfo[i] = rand_bit();
        }

        int info_idx = 0;
        for (int i = 0; i < g_n; i++)
            u[i] = frozen[i] ? 0 : sentInfo[info_idx++];

        polar_encode(u, x, g_n);

        for (int i = 0; i < g_n; i += 2) {
            double b0 = x[i] ? -1.0 : 1.0;
            double b1 = (i + 1 < g_n) ? (x[i + 1] ? -1.0 : 1.0) : 0.0;
            double z0, z1;
            rand_gaussian(&z0, &z1);
            rx[i] = b0 + sigma * z0;
            if (i + 1 < g_n) rx[i + 1] = b1 + sigma * z1;
        }

        for (int i = 0; i < g_n; i++) g_llr[0][i] = (2.0 / sigma2) * rx[i];

        double pm[MAX_LIST];
        pm[0] = 0.0;
        int path_count = 1;
        int origin_dummy[MAX_LIST];

        scl_decode_recursive(g_llr, g_ubits, g_n, 0, frozen, pm, &path_count,
                              cfg.l_limit, origin_dummy);

        int best_path = 0;
        if (cfg.mode == MODE_CASCL) {
            int found = -1;
            double best_pm = 1e300;
            for (int p = 0; p < path_count; p++) {
                uint8_t infoCand[K];
                int ii = 0;
                for (int i = 0; i < g_n; i++)
                    if (!frozen[i]) infoCand[ii++] = g_ubits[p][i];
                uint8_t crc_calc = crc8_compute(infoCand, K_INFO);
                uint8_t crc_rx = 0;
                for (int b = 0; b < CRC_LEN; b++)
                    crc_rx = (uint8_t)((crc_rx << 1) | infoCand[K_INFO + b]);
                if (crc_calc == crc_rx && pm[p] < best_pm) {
                    best_pm = pm[p];
                    found = p;
                }
            }
            if (found >= 0) {
                best_path = found;
            } else {
                double bpm = 1e300;
                for (int p = 0; p < path_count; p++)
                    if (pm[p] < bpm) { bpm = pm[p]; best_path = p; }
            }
        } else {
            double bpm = 1e300;
            for (int p = 0; p < path_count; p++)
                if (pm[p] < bpm) { bpm = pm[p]; best_path = p; }
        }

        uint8_t decodedInfo[K];
        int ii = 0;
        for (int i = 0; i < g_n; i++)
            if (!frozen[i]) decodedInfo[ii++] = g_ubits[best_path][i];

        long err_this_frame = 0;
        for (int i = 0; i < score_len; i++)
            if (decodedInfo[i] != sentInfo[i]) err_this_frame++;

        bit_errors += err_this_frame;
        if (err_this_frame > 0) frame_errors++;

        /* 오래 걸리는 포인트에서 "멈춘 게 아니다"를 보여주는 진행률 출력.
           프레임 카운트 모듈로 체크는 가볍고, 실제 출력은 최소
           PRINT_MIN_SECONDS초 간격으로만 하도록 추가로 제한한다. */
        if (frames % PRINT_EVERY_FRAMES == 0) {
            time_t now = time(NULL);
            if (difftime(now, t_last_print) >= PRINT_MIN_SECONDS) {
                printf("    ...%-10s @ %.2fdB running: frames=%-10ld errors=%ld/%ld  elapsed=%.0fs\n",
                       cfg.name, ebno_db, frames, bit_errors,
                       (long)g_target_errors, difftime(now, t_run_start));
                fflush(stdout);
                t_last_print = now;
            }
        }
    }

    RunResult r;

r.frames = frames;
r.bit_errors = bit_errors;

if (frames <= 0 || score_len <= 0 ||
    bit_errors < 0 || frame_errors < 0 ||
    frame_errors > frames ||
    bit_errors < frame_errors) {

    fprintf(stderr,
        "[ERROR] Invalid counters: decoder=%s, EbNo=%.2f, "
        "frames=%ld, bit_errors=%ld, frame_errors=%ld, score_len=%ld\n",
        cfg.name, ebno_db, frames, bit_errors,
        frame_errors, score_len);

    exit(EXIT_FAILURE);
}

r.ber = (double)bit_errors / ((double)frames * (double)score_len);
r.fer = (double)frame_errors / (double)frames;

fprintf(stderr,
    "[CHECK] %s EbNo=%.2f bit_errors=%ld frames=%ld BER=%.9e\n",
    cfg.name, ebno_db, bit_errors, frames, r.ber);

return r;
}

/* ------------------------------------------------------------------ */
/* Result data files; the native GUI renders the curves with Win32 GDI. */
/* ------------------------------------------------------------------ */
int main(int argc, char **argv) {
    /* Arguments: [snr_start] [snr_stop] [snr_step] [target_bit_errors] [N] [R] [output_dir]. */
    double snr_start = -1.0, snr_stop = 3.5, snr_step = 0.5;
    long target_errors = 1000;
    char *end = NULL;
    if (argc >= 2) { snr_start = strtod(argv[1], &end); if (!end || *end) { fprintf(stderr, "Invalid SNR start.\n"); return 2; } }
    if (argc >= 3) { snr_stop = strtod(argv[2], &end); if (!end || *end) { fprintf(stderr, "Invalid SNR stop.\n"); return 2; } }
    if (argc >= 4) { snr_step = strtod(argv[3], &end); if (!end || *end) { fprintf(stderr, "Invalid SNR step.\n"); return 2; } }
    if (argc >= 5) { target_errors = strtol(argv[4], &end, 10); if (!end || *end || target_errors < 1 || target_errors > 100000000L) { fprintf(stderr, "Target errors must be 1..100000000.\n"); return 2; } }
    if (argc >= 6) {
        long parsed_n = strtol(argv[5], &end, 10);
        if (!end || *end) { fprintf(stderr, "N must be a power of two from 64 to 1024.\n"); return 2; }
        double r_target = 0.5;
        if (argc >= 7) { r_target = strtod(argv[6], &end); if (!end || *end) { fprintf(stderr, "Invalid R.\n"); return 2; } }
        int new_k = 0, new_depth = 0;
        if (!valid_code_parameters((int)parsed_n, r_target, &new_k, &new_depth)) {
            fprintf(stderr, "Use N=64,128,256,512,1024 and R in 0.15..0.90; rounded K must exceed CRC length 8.\n"); return 2;
        }
        g_n = (int)parsed_n; g_k = new_k; g_n_stages = new_depth;
    }
    if (argc >= 8) {
        if (strlen(argv[7]) >= sizeof(g_output_dir)) { fprintf(stderr, "Output path is too long.\n"); return 2; }
        strcpy(g_output_dir, argv[7]);
    }
    if (argc > 8 || !isfinite(snr_start) || !isfinite(snr_stop) || !isfinite(snr_step) ||
        snr_start < -20.0 || snr_stop > 30.0 || snr_stop < snr_start || snr_step <= 0.0) {
        fprintf(stderr, "Usage: compare.exe [snr_start snr_stop snr_step target_errors [N [R [output_dir]]]], valid range -20..30 dB.\n");
        return 2;
    }
    g_ebno_start = snr_start;
    g_ebno_step = snr_step;
    g_num_snr_points = (int)floor((snr_stop - snr_start) / snr_step + 1e-9) + 1;
    g_target_errors = target_errors;
    if (g_num_snr_points < 1 || g_num_snr_points > MAX_RUNTIME_SNR_POINTS) {
        fprintf(stderr, "SNR sweep may contain at most %d points.\n", MAX_RUNTIME_SNR_POINTS);
        return 2;
    }
    if (MKDIR(g_output_dir) != 0) {
#if defined(_WIN32)
        DWORD attr = GetFileAttributesA(g_output_dir);
        if (attr == INVALID_FILE_ATTRIBUTES || !(attr & FILE_ATTRIBUTE_DIRECTORY)) { fprintf(stderr, "Could not create output directory: %s\n", g_output_dir); return 1; }
#else
        struct stat st; if (stat(g_output_dir, &st) != 0 || !S_ISDIR(st.st_mode)) { fprintf(stderr, "Could not create output directory: %s\n", g_output_dir); return 1; }
#endif
    }

    DecoderConfig configs[NUM_CONFIGS];
    int nc = 0;
    strcpy(configs[nc].name, "SC"); configs[nc].mode = MODE_SC; configs[nc].l_limit = 1; nc++;
    for (int i = 0; i < NUM_LIST_SIZES; i++) {
        sprintf(configs[nc].name, "SCL_L%d", LIST_SIZES[i]);
        configs[nc].mode = MODE_SCL;
        configs[nc].l_limit = LIST_SIZES[i];
        nc++;
    }
    for (int i = 0; i < NUM_LIST_SIZES; i++) {
        sprintf(configs[nc].name, "CASCL_L%d", LIST_SIZES[i]);
        configs[nc].mode = MODE_CASCL;
        configs[nc].l_limit = LIST_SIZES[i];
        nc++;
    }

    FILE *fps[NUM_CONFIGS];
    for (int c = 0; c < nc; c++) {
        char fn[4096];
        snprintf(fn, sizeof(fn), "%s/result_%.31s.txt", g_output_dir, configs[c].name);
        fps[c] = fopen(fn, "w");
        if (!fps[c]) { fprintf(stderr, "Could not open result file for writing: %s\n", fn); for (int j = 0; j < c; ++j) fclose(fps[j]); return 1; }
        fprintf(fps[c], "# EbNo_dB BER FER frames\n");
    }

    printf("N=%d, K=%d (actual R=%.5f), CA-SCL K_INFO=%d (R_eff=%.5f), CRC-%d (poly 0x%02X)\n",
           g_n, g_k, (double)g_k / g_n, K_INFO, (double)K_INFO / g_n, CRC_LEN, CRC_POLY);
    printf("Output directory: %s\n", g_output_dir);
    printf("TARGET_ERRORS=%ld per point, no MAX_FRAMES cap (runs until target is met)\n",
           (long)g_target_errors);
    printf("Decoders: ");
    for (int c = 0; c < nc; c++) printf("%s ", configs[c].name);
    printf("\n\n");

    time_t t_start = time(NULL);
    int total_points = g_num_snr_points * nc;
    int done = 0;

    for (int s = 0; s < g_num_snr_points; s++) {
        double ebno = g_ebno_start + s * g_ebno_step;
        int frozen[N];
        build_frozen_mask(ebno, (double)g_k / (double)g_n, frozen);

        uint64_t seed = 0x243F6A8885A308D3ULL + (uint64_t)s * 9999991ULL;

        for (int c = 0; c < nc; c++) {
            RunResult r = run_point(configs[c], ebno, frozen, seed);
            fprintf(fps[c], "%.2f %.6e %.6e %ld\n", ebno, r.ber, r.fer, r.frames);
            fflush(fps[c]);

            done++;
            double elapsed = difftime(time(NULL), t_start);
            printf("[%d/%d] %-10s @ %.2fdB  BER=%.3e  FER=%.3e  frames=%-10ld elapsed=%.0fs\n",
                   done, total_points, configs[c].name, ebno, r.ber, r.fer, r.frames, elapsed);
            fflush(stdout);
        }
    }

    for (int c = 0; c < nc; c++) fclose(fps[c]);

    printf("\n[DONE] BER/FER 데이터가 %s 폴더에 저장되었습니다.\n", g_output_dir);
    /* The GUI reads result_<decoder>.txt files from g_output_dir and draws the plots with Win32 GDI. */
    return 0;
}