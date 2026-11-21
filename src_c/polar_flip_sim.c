/*
 * polar_flip_sim.c
 * Polar code 시뮬레이터: SC / SC + (CRC 실패 시 d_min 패턴 뒤집기) / CRC-aided SCL 비교
 *
 * 컴파일:  gcc -O2 -o polar_flip_sim polar_flip_sim.c -lm
 * 사용법:
 *   ./polar_flip_sim N num den [L] [EbN0_start] [EbN0_stop] [step] [max_frames] [target_err] [seed]
 *       예) ./polar_flip_sim 256 1 3                      (L=8, 1.0~3.0dB, 0.5dB 간격)
 *           ./polar_flip_sim 512 1 2 8 1.0 3.0 0.5 200000 100 1
 *   ./polar_flip_sim check N num den                      (자체 검증 모드)
 *
 * 구성
 *   - 코드: G = F^{(x)n} (natural order), F=[1 0;1 1]. K = N*R + 24 (CRC-24C 포함)
 *   - frozen set: polar_dmin.c 와 같은 beta-expansion (5G NR 방식)
 *   - 정보 벡터(오름차순 인덱스 순): [메시지 N*R 비트 | CRC 24비트]
 *   - 제안 방식(flip): SC 결과의 CRC가 실패하면, 최소 weight(d_min) 코드워드에 해당하는
 *       u-도메인 오류 패턴 e를 뒤집어서 CRC가 통과하는지 확인한다.
 *       CRC는 선형이므로 "e를 뒤집고 CRC 재확인" = "CRC 신드롬(û) == CRC 신드롬(e)인 e를 찾기".
 *       그래서 패턴마다 24비트 신드롬을 미리 저장해 둔다 (이것이 메모리가 필요한 이유).
 *   - FER / BER은 메시지 비트(CRC 제외) 기준
 *   - 결과는 ../results/perf_N{N}_R{num}_{den}.csv (폴더가 없으면 현재 폴더)에 저장
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#ifdef _WIN32
#include <windows.h>        /* 윈도우 콘솔에서 한글이 깨지지 않도록 UTF-8 출력 설정에 사용 */
#endif

#define MAXN   4096
#define MAXLOG 12
#define MAXL   32
#define CRC_LEN 24
#define CRC_POLY 0xB2B117u          /* CRC-24C (5G NR) */
#define PI_CONST 3.14159265358979323846

static int N, n, K, Kmsg;
static uint8_t frozen[MAXN];
static int infoIdx[MAXN];           /* 정보 비트의 u 인덱스 (오름차순) */

/* ===================== 난수 ===================== */
static uint64_t rng_s = 88172645463325252ULL;
static inline uint64_t rng_next(void)
{
    rng_s ^= rng_s >> 12; rng_s ^= rng_s << 25; rng_s ^= rng_s >> 27;
    return rng_s * 2685821657736338717ULL;
}
static inline double rng_u01(void) { return ((rng_next() >> 11) + 0.5) * (1.0 / 9007199254740992.0); }
static double rng_gauss(void)
{
    static int has = 0; static double spare;
    if (has) { has = 0; return spare; }
    double u1 = rng_u01(), u2 = rng_u01();
    double r = sqrt(-2.0 * log(u1));
    spare = r * sin(2.0 * PI_CONST * u2); has = 1;
    return r * cos(2.0 * PI_CONST * u2);
}

/* ===================== frozen set (beta-expansion) ===================== */
typedef struct { int idx; double z; } Rel;
static int cmp_rel(const void *a, const void *b)
{
    const Rel *x = (const Rel *)a, *y = (const Rel *)b;
    if (x->z > y->z) return -1;
    if (x->z < y->z) return 1;
    return y->idx - x->idx;
}
static void construct_code(void)
{
    Rel *r = (Rel *)malloc(sizeof(Rel) * N);
    double beta = pow(2.0, 0.25);
    for (int i = 0; i < N; i++) {
        double pw = 0.0;
        for (int j = 0; j < n; j++) if ((i >> j) & 1) pw += pow(beta, j);
        r[i].idx = i; r[i].z = pw;
    }
    qsort(r, N, sizeof(Rel), cmp_rel);
    for (int i = 0; i < N; i++) frozen[i] = 1;
    for (int k = 0; k < K; k++) frozen[r[k].idx] = 0;
    int c = 0;
    for (int i = 0; i < N; i++) if (!frozen[i]) infoIdx[c++] = i;
    free(r);
}

/* ===================== CRC / 인코더 ===================== */
/* 비트열 b[0..len-1]을 (MSB 먼저) x^24 곱한 다항식을 g(x)로 나눈 나머지. 선형(init=0, xorout=0) */
static uint32_t crc_rem(const uint8_t *b, int len)
{
    uint32_t r = 0;
    for (int i = 0; i < len; i++) {
        uint32_t fb = ((r >> 23) & 1u) ^ b[i];
        r = (r << 1) & 0xFFFFFFu;
        if (fb) r ^= CRC_POLY;
    }
    return r;
}

static void polar_transform(uint8_t *x)       /* x = x * F^{(x)n} (자기 자신이 역변환) */
{
    for (int len = 1; len < N; len <<= 1)
        for (int i = 0; i < N; i += 2 * len)
            for (int j = i; j < i + len; j++) x[j] ^= x[j + len];
}

/* ===================== SC 디코더 (exact LLR) ===================== */
static inline double f_exact(double a, double b)
{
    double aa = fabs(a), bb = fabs(b);
    double mn = aa < bb ? aa : bb;
    double s = ((a < 0) != (b < 0)) ? -1.0 : 1.0;
    return s * mn + log1p(exp(-fabs(a + b))) - log1p(exp(-fabs(a - b)));
}

static double  *scL[MAXLOG + 1];
static uint8_t *scXl[MAXLOG + 1], *scXr[MAXLOG + 1];
static uint8_t sc_uhat[MAXN];
static double  sc_Ldec[MAXN];

static void sc_alloc(void)
{
    for (int d = 0; d < n; d++) {
        int h = N >> (d + 1);
        scL[d]  = (double *)malloc(sizeof(double) * h);
        scXl[d] = (uint8_t *)malloc(h);
        scXr[d] = (uint8_t *)malloc(h);
    }
}

static void sc_rec(int m, int depth, const double *L, int base, uint8_t *xout)
{
    if (m == 1) {
        double l = L[0];
        uint8_t u = frozen[base] ? 0 : (l < 0);
        sc_uhat[base] = u; sc_Ldec[base] = l; xout[0] = u;
        return;
    }
    int h = m >> 1;
    double *buf = scL[depth];
    uint8_t *xl = scXl[depth], *xr = scXr[depth];
    for (int i = 0; i < h; i++) buf[i] = f_exact(L[i], L[i + h]);
    sc_rec(h, depth + 1, buf, base, xl);
    for (int i = 0; i < h; i++) buf[i] = L[i + h] + (xl[i] ? -L[i] : L[i]);
    sc_rec(h, depth + 1, buf, base + h, xr);
    for (int i = 0; i < h; i++) { xout[i] = xl[i] ^ xr[i]; xout[i + h] = xr[i]; }
}

static void sc_decode(const double *llr)
{
    static uint8_t xtmp[MAXN];
    sc_rec(N, 0, llr, 0, xtmp);
}

/* ===================== SCL 디코더 (Tal-Vardy lazy copy, LLR) ===================== */
static int Lmax = 8;
static double  *arrP[MAXLOG + 1][MAXL];
static uint8_t *arrC[MAXLOG + 1][MAXL];
static int pathToArr[MAXLOG + 1][MAXL];
static int inactiveArr[MAXLOG + 1][MAXL], inactiveArrCnt[MAXLOG + 1];
static int arrRef[MAXLOG + 1][MAXL];
static int inactivePath[MAXL], inactivePathCnt;
static uint8_t activePath[MAXL];
static double  PM[MAXL];
static uint8_t scl_uh[MAXL][MAXN];

static void scl_alloc(void)
{
    for (int lam = 0; lam <= n; lam++) {
        int len = 1 << (n - lam);
        for (int s = 0; s < Lmax; s++) {
            arrP[lam][s] = (double *)malloc(sizeof(double) * len);
            arrC[lam][s] = (uint8_t *)malloc(2 * len);
        }
    }
}

static void scl_reset(void)
{
    inactivePathCnt = 0;
    for (int l = Lmax - 1; l >= 0; l--) { inactivePath[inactivePathCnt++] = l; activePath[l] = 0; }
    for (int lam = 0; lam <= n; lam++) {
        inactiveArrCnt[lam] = 0;
        for (int s = Lmax - 1; s >= 0; s--) { inactiveArr[lam][inactiveArrCnt[lam]++] = s; arrRef[lam][s] = 0; }
    }
}

static int assign_initial_path(void)
{
    int l = inactivePath[--inactivePathCnt];
    activePath[l] = 1;
    for (int lam = 0; lam <= n; lam++) {
        int s = inactiveArr[lam][--inactiveArrCnt[lam]];
        pathToArr[lam][l] = s; arrRef[lam][s] = 1;
    }
    return l;
}

static int clone_path(int l)
{
    int c = inactivePath[--inactivePathCnt];
    activePath[c] = 1; PM[c] = PM[l];
    for (int lam = 0; lam <= n; lam++) {
        int s = pathToArr[lam][l];
        pathToArr[lam][c] = s; arrRef[lam][s]++;
    }
    return c;
}

static void kill_path(int l)
{
    activePath[l] = 0; inactivePath[inactivePathCnt++] = l;
    for (int lam = 0; lam <= n; lam++) {
        int s = pathToArr[lam][l];
        if (--arrRef[lam][s] == 0) inactiveArr[lam][inactiveArrCnt[lam]++] = s;
    }
}

/* (layer, path)의 배열을 이 경로 전용으로 만들고 배열 번호 반환 (필요하면 복사) */
static int own_arr(int lam, int l)
{
    int s = pathToArr[lam][l];
    if (arrRef[lam][s] == 1) return s;
    int s2 = inactiveArr[lam][--inactiveArrCnt[lam]];
    int len = 1 << (n - lam);
    memcpy(arrP[lam][s2], arrP[lam][s], sizeof(double) * len);
    memcpy(arrC[lam][s2], arrC[lam][s], 2 * len);
    arrRef[lam][s]--; arrRef[lam][s2] = 1;
    pathToArr[lam][l] = s2;
    return s2;
}

static void calc_P(int lam, int phi)
{
    if (lam == 0) return;
    int psi = phi >> 1;
    if ((phi & 1) == 0) calc_P(lam - 1, psi);
    int len = 1 << (n - lam);
    for (int l = 0; l < Lmax; l++) {
        if (!activePath[l]) continue;
        int sl = own_arr(lam, l), sm = own_arr(lam - 1, l);
        double *Pl = arrP[lam][sl], *Pm = arrP[lam - 1][sm];
        uint8_t *Cl = arrC[lam][sl];
        if ((phi & 1) == 0) {
            for (int b = 0; b < len; b++) Pl[b] = f_exact(Pm[2 * b], Pm[2 * b + 1]);
        } else {
            for (int b = 0; b < len; b++) {
                uint8_t u = Cl[2 * b];
                Pl[b] = (u ? -Pm[2 * b] : Pm[2 * b]) + Pm[2 * b + 1];
            }
        }
    }
}

static void update_C(int lam, int phi)
{
    int psi = phi >> 1;
    int len = 1 << (n - lam);
    for (int l = 0; l < Lmax; l++) {
        if (!activePath[l]) continue;
        int sl = own_arr(lam, l), sm = own_arr(lam - 1, l);
        uint8_t *Cl = arrC[lam][sl], *Cm = arrC[lam - 1][sm];
        for (int b = 0; b < len; b++) {
            Cm[2 * (2 * b)     + (psi & 1)] = Cl[2 * b] ^ Cl[2 * b + 1];
            Cm[2 * (2 * b + 1) + (psi & 1)] = Cl[2 * b + 1];
        }
    }
    if (psi & 1) update_C(lam - 1, psi);
}

static inline double pm_pen(double L, int u)      /* ln(1 + exp(-(1-2u) L)) */
{
    double x = -(1 - 2 * u) * L;
    return x > 0 ? x + log1p(exp(-x)) : log1p(exp(x));
}

static void continue_paths_unfrozen(int phi)
{
    double fm[MAXL][2]; int cont[MAXL][2];
    int cnt = 0;
    for (int l = 0; l < Lmax; l++) {
        if (!activePath[l]) continue;
        int s = own_arr(n, l);
        double L0 = arrP[n][s][0];
        fm[l][0] = PM[l] + pm_pen(L0, 0);
        fm[l][1] = PM[l] + pm_pen(L0, 1);
        cont[l][0] = cont[l][1] = 0;
        cnt++;
    }
    int rho = 2 * cnt < Lmax ? 2 * cnt : Lmax;
    /* 가장 작은 metric rho개 선택 */
    for (int t = 0; t < rho; t++) {
        double best = 1e300; int bl = -1, bb = 0;
        for (int l = 0; l < Lmax; l++) {
            if (!activePath[l]) continue;
            for (int b = 0; b < 2; b++)
                if (!cont[l][b] && fm[l][b] < best) { best = fm[l][b]; bl = l; bb = b; }
        }
        cont[bl][bb] = 1;
    }
    for (int l = 0; l < Lmax; l++)
        if (activePath[l] && !cont[l][0] && !cont[l][1]) kill_path(l);
    uint8_t was[MAXL];                      /* 복제로 새로 생긴 경로는 이번 루프에서 건너뛴다 */
    for (int l = 0; l < Lmax; l++) was[l] = activePath[l];
    for (int l = 0; l < Lmax; l++) {
        if (!was[l]) continue;
        if (cont[l][0] && cont[l][1]) {
            int s = own_arr(n, l);
            arrC[n][s][phi & 1] = 0;
            int c = clone_path(l);
            int s2 = own_arr(n, c);
            arrC[n][s2][phi & 1] = 1;
            memcpy(scl_uh[c], scl_uh[l], phi);
            scl_uh[l][phi] = 0; scl_uh[c][phi] = 1;
            PM[l] = fm[l][0]; PM[c] = fm[l][1];
        } else {
            int b = cont[l][1] ? 1 : 0;
            int s = own_arr(n, l);
            arrC[n][s][phi & 1] = (uint8_t)b;
            scl_uh[l][phi] = (uint8_t)b;
            PM[l] = fm[l][b];
        }
    }
}

static uint8_t scl_out[MAXN];      /* 선택된 경로의 u 추정 */
static int bitrev_tab[MAXN];

/* 반환: 1이면 CRC 통과한 경로 선택, 0이면 CRC 통과 경로 없음(PM 최소 경로 선택) */
static int scl_decode(const double *llr)
{
    scl_reset();
    int l0 = assign_initial_path();
    PM[l0] = 0.0;
    int s0 = own_arr(0, l0);
    for (int j = 0; j < N; j++) arrP[0][s0][j] = llr[bitrev_tab[j]];
    for (int phi = 0; phi < N; phi++) {
        calc_P(n, phi);
        if (frozen[phi]) {
            for (int l = 0; l < Lmax; l++) {
                if (!activePath[l]) continue;
                int s = own_arr(n, l);
                arrC[n][s][phi & 1] = 0;
                PM[l] += pm_pen(arrP[n][s][0], 0);
                scl_uh[l][phi] = 0;
            }
        } else {
            continue_paths_unfrozen(phi);
        }
        if (phi & 1) update_C(n, phi);
    }
    /* 경로 선택: PM 오름차순으로 CRC 통과하는 첫 경로 */
    int order[MAXL], na = 0;
    for (int l = 0; l < Lmax; l++) if (activePath[l]) order[na++] = l;
    for (int i = 1; i < na; i++) {          /* 삽입 정렬 */
        int v = order[i], j = i - 1;
        while (j >= 0 && PM[order[j]] > PM[v]) { order[j + 1] = order[j]; j--; }
        order[j + 1] = v;
    }
    static uint8_t tmp[MAXN];
    int chosen = order[0], ok = 0;
    for (int i = 0; i < na; i++) {
        for (int k = 0; k < K; k++) tmp[k] = scl_uh[order[i]][infoIdx[k]];
        if (crc_rem(tmp, K) == 0) { chosen = order[i]; ok = 1; break; }
    }
    memcpy(scl_out, scl_uh[chosen], N);
    return ok;
}

/* ===================== d_min 패턴 생성 ===================== */
typedef struct {
    int A, nw;
    uint64_t *bits;       /* A x nw : 정보 비트 순서(K비트)로 압축한 u-도메인 오류 패턴 */
    uint32_t *syn;        /* 패턴의 CRC 신드롬 */
    int *start;           /* 패턴이 시작하는 u 인덱스 */
    int dmin_log;         /* log2(d_min) */
    int Nr;               /* 시작 인덱스 개수 */
} PatSet;
static PatSet pat;

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static int cmp_u32(const void *a, const void *b)
{
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return (x > y) - (x < y);
}

/* 반환: 자체 검증 오류 개수 (0이면 정상) */
static int build_patterns(int verbose)
{
    int errors = 0;
    int wmin = n + 1;
    for (int k = 0; k < K; k++) { int w = __builtin_popcount(infoIdx[k]); if (w < wmin) wmin = w; }
    pat.dmin_log = wmin;

    unsigned long long Aexp = 0; int Nr = 0;
    for (int k = 0; k < K; k++) {
        int i = infoIdx[k];
        if (__builtin_popcount(i) != wmin) continue;
        int inv = 0;
        for (int l = 0; l < n; l++) if (!((i >> l) & 1)) inv += __builtin_popcount(i & ((1 << l) - 1));
        Aexp += 1ULL << ((n - wmin) + inv); Nr++;
    }
    if (Aexp > 4000000ULL) { fprintf(stderr, "패턴 수(%llu)가 너무 많아 생성하지 않습니다\n", Aexp); pat.A = 0; return 1; }
    pat.A = (int)Aexp; pat.Nr = Nr; pat.nw = (K + 63) / 64;
    pat.bits  = (uint64_t *)calloc((size_t)pat.A * pat.nw, sizeof(uint64_t));
    pat.syn   = (uint32_t *)malloc(sizeof(uint32_t) * pat.A);
    pat.start = (int *)malloc(sizeof(int) * pat.A);

    static uint8_t w[MAXN], kb[MAXN];
    int a = 0;
    for (int kk = 0; kk < K; kk++) {
        int i = infoIdx[kk];
        if (__builtin_popcount(i) != wmin) continue;
        int S[MAXLOG], Z[MAXLOG], ns = 0, nz = 0;
        for (int l = 0; l < n; l++) { if ((i >> l) & 1) S[ns++] = l; else Z[nz++] = l; }
        /* 매개변수: (k in S, l in Z, l > k) 쌍마다 1비트, 그리고 Z 위치마다 평행이동 1비트 */
        int pk[MAXLOG * MAXLOG], pl[MAXLOG * MAXLOG], np = 0;
        for (int s = 0; s < ns; s++) for (int z = 0; z < nz; z++) if (Z[z] > S[s]) { pk[np] = s; pl[np] = Z[z]; np++; }
        uint64_t total = 1ULL << (np + nz);
        for (uint64_t mask = 0; mask < total; mask++) {
            int v[MAXLOG];
            for (int s = 0; s < ns; s++) v[s] = 1 << S[s];
            for (int p = 0; p < np; p++) if ((mask >> p) & 1) v[pk[p]] |= 1 << pl[p];
            int b = 0;
            for (int z = 0; z < nz; z++) if ((mask >> (np + z)) & 1) b |= 1 << Z[z];
            memset(w, 0, N);
            for (int t = 0; t < (1 << ns); t++) {          /* flat의 모든 점 */
                int p = b;
                for (int s = 0; s < ns; s++) if ((t >> s) & 1) p ^= v[s];
                w[p] = 1;
            }
            polar_transform(w);                              /* 코드워드 -> u-도메인 패턴 */
            for (int j = 0; j < N; j++) if (frozen[j] && w[j]) errors++;   /* frozen 위치에 1이 있으면 안 됨 */
            for (int k = 0; k < K; k++) kb[k] = w[infoIdx[k]];
            int first = -1;
            for (int k = 0; k < K; k++) if (kb[k]) { first = infoIdx[k]; break; }
            if (first != i) errors++;
            uint64_t *dst = pat.bits + (size_t)a * pat.nw;
            for (int k = 0; k < K; k++) if (kb[k]) dst[k >> 6] |= 1ULL << (k & 63);
            pat.syn[a] = crc_rem(kb, K);
            pat.start[a] = i;
            a++;
        }
    }
    if (a != pat.A) errors++;
    if (verbose) {
        /* 중복 검사 (해시) */
        uint64_t *h = (uint64_t *)malloc(sizeof(uint64_t) * pat.A);
        for (int q = 0; q < pat.A; q++) {
            uint64_t x = 1469598103934665603ULL;
            for (int t = 0; t < pat.nw; t++) { x ^= pat.bits[(size_t)q * pat.nw + t]; x *= 1099511628211ULL; x ^= x >> 29; }
            h[q] = x;
        }
        qsort(h, pat.A, sizeof(uint64_t), cmp_u64);
        int dup = 0; for (int q = 1; q < pat.A; q++) if (h[q] == h[q - 1]) dup++;
        free(h);
        if (dup) errors++;
        /* 신드롬 충돌(서로 다른 패턴이 같은 신드롬) 개수 */
        uint32_t *sy = (uint32_t *)malloc(sizeof(uint32_t) * pat.A);
        memcpy(sy, pat.syn, sizeof(uint32_t) * pat.A);
        qsort(sy, pat.A, sizeof(uint32_t), cmp_u32);
        int coll = 0; for (int q = 1; q < pat.A; q++) if (sy[q] == sy[q - 1]) coll++;
        free(sy);
        printf("[패턴] d_min=%d, A_dmin=%d (공식 %llu), Nr=%d, 중복=%d, 신드롬 충돌=%d, 자체검증 오류=%d\n",
               1 << wmin, pat.A, Aexp, Nr, dup, coll, errors);
        printf("[메모리] 패턴 저장: %.1f KB (패턴 %d개 x %d바이트 + 신드롬 4바이트)\n",
               ((double)pat.A * (pat.nw * 8 + 4 + 4)) / 1024.0, pat.A, pat.nw * 8);
    }
    return errors;
}

/* CRC 실패 시 패턴으로 보정. 반환: 0 = CRC 이미 통과, 1 = 보정 적용, -1 = 매칭 패턴 없음 */
static int flip_decode(uint8_t *uinfo, const double *Ldec)
{
    uint32_t s = crc_rem(uinfo, K);
    if (s == 0) return 0;
    int best = -1; double bm = 1e300;
    for (int a = 0; a < pat.A; a++) {
        if (pat.syn[a] != s) continue;
        const uint64_t *pb = pat.bits + (size_t)a * pat.nw;
        double m = 0.0;
        for (int k = 0; k < K; k++) if ((pb[k >> 6] >> (k & 63)) & 1) m += fabs(Ldec[infoIdx[k]]);
        if (m < bm) { bm = m; best = a; }
    }
    if (best < 0) return -1;
    const uint64_t *pb = pat.bits + (size_t)best * pat.nw;
    for (int k = 0; k < K; k++) if ((pb[k >> 6] >> (k & 63)) & 1) uinfo[k] ^= 1;
    return 1;
}

/* ===================== 자체 검증 ===================== */
static int self_check(void)
{
    int bad = 0;
    printf("== 자체 검증 (N=%d, K=%d) ==\n", N, K);
    bad += build_patterns(1);

    /* SC 와 SCL(L=1) 이 같은 결정을 내리는지 (랜덤 LLR) */
    int saveL = Lmax; Lmax = 1;
    for (int d = 0; d <= n; d++) {
        for (int s = 0; s < MAXL; s++) { if (arrP[d][s]) { free(arrP[d][s]); free(arrC[d][s]); arrP[d][s] = NULL; } }
    }
    scl_alloc();
    int mism = 0, trials = 300;
    double *llr = (double *)calloc(N, sizeof(double));
    for (int t = 0; t < trials; t++) {
        for (int j = 0; j < N; j++) llr[j] = 3.0 * rng_gauss() + ((t & 1) ? 1.0 : -0.5);
        sc_decode(llr);
        scl_decode(llr);
        for (int j = 0; j < N; j++) if (sc_uhat[j] != scl_out[j]) { mism++; break; }
    }
    printf("[SC vs SCL(L=1)] %d회 중 불일치 %d회\n", trials, mism);
    bad += mism;
    free(llr);

    /* 패턴 하나를 직접 오류로 넣었을 때 flip이 복구하는지 */
    int fixed = 0, tot = 0;
    uint8_t *ui = (uint8_t *)malloc(K), *ue = (uint8_t *)malloc(K);
    double *Ld = (double *)calloc(MAXN, sizeof(double));
    for (int a = 0; a < pat.A && tot < 2000; a += (pat.A / 1000) + 1) {
        for (int k = 0; k < Kmsg; k++) ui[k] = rng_next() & 1;
        uint32_t c = crc_rem(ui, Kmsg);
        for (int j = 0; j < CRC_LEN; j++) ui[Kmsg + j] = (c >> (CRC_LEN - 1 - j)) & 1;
        memcpy(ue, ui, K);
        const uint64_t *pb = pat.bits + (size_t)a * pat.nw;
        for (int k = 0; k < K; k++) if ((pb[k >> 6] >> (k & 63)) & 1) ue[k] ^= 1;
        for (int j = 0; j < N; j++) Ld[j] = 1.0;
        flip_decode(ue, Ld);
        tot++; if (memcmp(ue, ui, K) == 0) fixed++;
    }
    printf("[패턴 주입 복구] %d개 중 %d개 복구 (신드롬 충돌로 일부는 다른 패턴이 선택될 수 있음)\n", tot, fixed);
    if (fixed < tot * 0.9) bad++;
    free(ui); free(ue); free(Ld);
    Lmax = saveL;
    printf("== 결과: %s ==\n", bad ? "문제 있음" : "정상");
    return bad;
}

/* ===================== 메인 ===================== */
int main(int argc, char **argv)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001);      /* 콘솔 출력을 UTF-8로 */
#endif
    int argi = 1, check = 0;
    if (argc > 1 && strcmp(argv[1], "check") == 0) { check = 1; argi = 2; }
    if (argc - argi < 3) {
        fprintf(stderr, "사용법: %s [check] N num den [L] [EbN0_start] [EbN0_stop] [step] [max_frames] [target_err] [seed]\n", argv[0]);
        return 1;
    }
    N = atoi(argv[argi]); int num = atoi(argv[argi + 1]), den = atoi(argv[argi + 2]);
    int Lval = (argc > argi + 3) ? atoi(argv[argi + 3]) : 8;
    double e0 = (argc > argi + 4) ? atof(argv[argi + 4]) : 1.0;
    double e1 = (argc > argi + 5) ? atof(argv[argi + 5]) : 3.0;
    double es = (argc > argi + 6) ? atof(argv[argi + 6]) : 0.5;
    long maxf = (argc > argi + 7) ? atol(argv[argi + 7]) : 100000;
    int  targ = (argc > argi + 8) ? atoi(argv[argi + 8]) : 100;
    uint64_t seed = (argc > argi + 9) ? (uint64_t)atoll(argv[argi + 9]) : 1;
    rng_s ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 20; i++) rng_next();

    n = 0; while ((1 << n) < N) n++;
    if ((1 << n) != N || N > MAXN || Lval < 1 || Lval > MAXL) { fprintf(stderr, "N은 2의 거듭제곱(<=%d), 1<=L<=%d\n", MAXN, MAXL); return 1; }
    Lmax = Lval;
    Kmsg = (int)floor((double)N * num / den + 0.5);
    K = Kmsg + CRC_LEN;
    construct_code();
    for (int j = 0; j < N; j++) { int r = 0, x = j; for (int b = 0; b < n; b++) { r = (r << 1) | (x & 1); x >>= 1; } bitrev_tab[j] = r; }
    sc_alloc(); scl_alloc();

    if (check) return self_check() ? 2 : 0;

    printf("N=%d R=%d/%d K=%d (메시지 %d + CRC %d), SCL L=%d\n", N, num, den, K, Kmsg, CRC_LEN, Lmax);
    if (build_patterns(1)) { fprintf(stderr, "패턴 생성 오류\n"); return 2; }

    char fname[128];
    snprintf(fname, sizeof(fname), "../results/perf_N%d_R%d_%d.csv", N, num, den);
    FILE *fp = fopen(fname, "w");
    if (!fp) { snprintf(fname, sizeof(fname), "perf_N%d_R%d_%d.csv", N, num, den); fp = fopen(fname, "w"); }
    if (!fp) { perror("fopen"); return 1; }
    fprintf(fp, "ebn0_db,frames,sc_fer,sc_ber,flip_fer,flip_ber,scl_fer,scl_ber,"
                "sc_err,sc_err_crc_detected,flip_fixed,flip_miscorrected,flip_nomatch,flip_err,scl_err\n");

    double Reff = (double)Kmsg / N;
    uint8_t *uinfo = (uint8_t *)malloc(K), *uh = (uint8_t *)malloc(K);
    uint8_t *u = (uint8_t *)malloc(N), *x = (uint8_t *)malloc(N);
    double *llr = (double *)malloc(sizeof(double) * N);

    printf("%-6s %-9s | %-9s %-9s | %-9s %-9s | %-9s %-9s | (SC오류 CRC검출 복구 오정정 매칭없음)\n",
           "Eb/N0", "frames", "SC FER", "SC BER", "FLIP FER", "FLIP BER", "SCL FER", "SCL BER");
    for (double eb = e0; eb <= e1 + 1e-9; eb += es) {
        double sigma2 = 1.0 / (2.0 * Reff * pow(10.0, eb / 10.0)), sigma = sqrt(sigma2);
        long frames = 0;
        long sc_err = 0, flip_err = 0, scl_err = 0, sc_bit = 0, flip_bit = 0, scl_bit = 0;
        long det = 0, fixed = 0, miscorr = 0, nomatch = 0;
        while (frames < maxf && flip_err < targ) {
            frames++;
            for (int k = 0; k < Kmsg; k++) uinfo[k] = rng_next() & 1;
            uint32_t c = crc_rem(uinfo, Kmsg);
            for (int j = 0; j < CRC_LEN; j++) uinfo[Kmsg + j] = (c >> (CRC_LEN - 1 - j)) & 1;
            memset(u, 0, N);
            for (int k = 0; k < K; k++) u[infoIdx[k]] = uinfo[k];
            memcpy(x, u, N); polar_transform(x);
            for (int j = 0; j < N; j++) llr[j] = 2.0 * ((1.0 - 2.0 * x[j]) + sigma * rng_gauss()) / sigma2;

            /* --- SC --- */
            sc_decode(llr);
            for (int k = 0; k < K; k++) uh[k] = sc_uhat[infoIdx[k]];
            int e_sc = 0; for (int k = 0; k < Kmsg; k++) e_sc += (uh[k] != uinfo[k]);
            sc_bit += e_sc; if (e_sc) sc_err++;
            int crc_fail = (crc_rem(uh, K) != 0);

            /* --- SC + 뒤집기 --- */
            int rc = flip_decode(uh, sc_Ldec);
            int e_fl = 0; for (int k = 0; k < Kmsg; k++) e_fl += (uh[k] != uinfo[k]);
            flip_bit += e_fl; if (e_fl) flip_err++;
            if (e_sc && crc_fail) det++;
            if (e_sc && !e_fl) fixed++;
            if (rc == 1 && e_fl) miscorr++;
            if (rc == -1) nomatch++;

            /* --- SCL --- */
            scl_decode(llr);
            int e_sl = 0; for (int k = 0; k < Kmsg; k++) e_sl += (scl_out[infoIdx[k]] != uinfo[k]);
            scl_bit += e_sl; if (e_sl) scl_err++;
        }
        double fn = (double)frames, bn = (double)frames * Kmsg;
        printf("%-6.2f %-9ld | %-9.2e %-9.2e | %-9.2e %-9.2e | %-9.2e %-9.2e | %ld %ld %ld %ld %ld\n",
               eb, frames, sc_err / fn, sc_bit / bn, flip_err / fn, flip_bit / bn, scl_err / fn, scl_bit / bn,
               sc_err, det, fixed, miscorr, nomatch);
        fprintf(fp, "%.2f,%ld,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n",
                eb, frames, sc_err / fn, sc_bit / bn, flip_err / fn, flip_bit / bn, scl_err / fn, scl_bit / bn,
                sc_err, det, fixed, miscorr, nomatch, flip_err, scl_err);
        fflush(fp); fflush(stdout);
    }
    fclose(fp);
    printf("결과 저장: %s\n", fname);
    return 0;
}
