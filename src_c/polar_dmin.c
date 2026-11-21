/*
 * polar_dmin.c
 * Polar code의 d_min, A_dmin(최소 weight 코드워드 개수)와
 * "정보 비트 인덱스별 최소 weight 코드워드 개수"를 계산해서 CSV로 저장한다.
 * (그래프는 matlab/plot_dmin.m 으로 MATLAB에서 그린다)
 *
 * 컴파일:  gcc -O2 -o polar_dmin polar_dmin.c -lm
 * 사용법:
 *   ./polar_dmin                          -> N=256,512,1024,2048 x R=1/2,1/3 전부 실행
 *   ./polar_dmin N num den [CRC] [beta|ga] [EbN0dB(ga일 때만)]
 *       예) ./polar_dmin 256 1 3            (CRC=24, beta-expansion)
 *           ./polar_dmin 256 1 3 24 ga 1.0  (CRC=24, Gaussian Approximation, 설계 Eb/N0=1.0dB)
 *
 * 가정 / 약속
 *  - G = F^{(x)n}, F=[1 0;1 1] (natural order). 행 i의 weight = 2^popcount(i)
 *  - K = round(N*R) + CRC (기본 CRC=24, CRC 비트도 정보 집합에 포함)
 *  - frozen set(채널 신뢰도 순서)
 *      beta : 5G NR 방식 beta-expansion.  PW(i) = sum_j b_j(i) * 2^(j/4)  (b_j: i의 j번째 비트)
 *             -> 선배님 그래프 8가지 케이스(Nr, A_dmin, d_min)와 모두 일치 (기본값)
 *      ga   : Gaussian Approximation (설계 Eb/N0에 따라 결과가 바뀜, 선배님 결과와는 다름)
 *  - 출력 인덱스는 1-based (MATLAB 그래프와 맞추기 위해)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI_CONST 3.14159265358979323846

enum { METHOD_BETA = 0, METHOD_GA = 1 };

typedef struct { int idx; double z; } Rel;

/* ---------- 신뢰도 정렬 공통 ---------- */
static int cmp_rel(const void *a, const void *b)   /* z 내림차순, 동점이면 큰 인덱스 먼저 */
{
    const Rel *x = (const Rel *)a, *y = (const Rel *)b;
    if (x->z > y->z) return -1;
    if (x->z < y->z) return 1;
    return y->idx - x->idx;
}

/* ---------- 방법 1: beta-expansion (5G NR) ---------- */
static void reliability_beta(int N, Rel *r)
{
    int n = 0; while ((1 << n) < N) n++;
    double beta = pow(2.0, 0.25);
    for (int i = 0; i < N; i++) {
        double pw = 0.0;
        for (int j = 0; j < n; j++)
            if ((i >> j) & 1) pw += pow(beta, j);
        r[i].idx = i; r[i].z = pw;
    }
}

/* ---------- 방법 2: Gaussian Approximation ---------- */
static double phi(double x)
{
    if (x < 10.0) return exp(-0.4527 * pow(x, 0.86) + 0.0218);
    return sqrt(PI_CONST / x) * exp(-x / 4.0) * (1.0 - 10.0 / (7.0 * x));
}

static double phi_inv(double y)          /* phi는 감소함수 -> 이분법 */
{
    double lo = 1e-10, hi = 1e5;
    if (y >= phi(lo)) return lo;
    if (y <= phi(hi)) return hi;
    for (int it = 0; it < 100; it++) {
        double mid = 0.5 * (lo + hi);
        if (phi(mid) > y) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

static int bitrev(int x, int n)
{
    int r = 0;
    for (int i = 0; i < n; i++) { r = (r << 1) | (x & 1); x >>= 1; }
    return r;
}

static void reliability_ga(int N, double R, double ebn0_db, Rel *r)
{
    int n = 0; while ((1 << n) < N) n++;
    double *z = (double *)malloc(sizeof(double) * N);
    double sigma2 = 1.0 / (2.0 * R * pow(10.0, ebn0_db / 10.0));
    z[0] = 2.0 / sigma2;                                   /* 채널 LLR 평균 */
    for (int lev = 0; lev < n; lev++) {
        int B = 1 << lev;
        for (int j = 0; j < B; j++) {
            double T = z[j];
            z[j]     = phi_inv(1.0 - pow(1.0 - phi(T), 2.0));  /* 나쁜 채널 */
            z[j + B] = 2.0 * T;                                /* 좋은 채널 */
        }
    }
    for (int i = 0; i < N; i++) { r[i].idx = i; r[i].z = z[bitrev(i, n)]; }  /* natural order */
    free(z);
}

/* 신뢰도 상위 K개를 정보 비트(1), 나머지를 frozen(0)으로 */
static void construct(int N, int K, int method, double ebn0_db, int *inA)
{
    Rel *r = (Rel *)malloc(sizeof(Rel) * N);
    if (method == METHOD_BETA) reliability_beta(N, r);
    else                       reliability_ga(N, (double)K / N, ebn0_db, r);
    qsort(r, N, sizeof(Rel), cmp_rel);
    for (int i = 0; i < N; i++) inA[i] = 0;
    for (int k = 0; k < K; k++) inA[r[k].idx] = 1;
    free(r);
}

/* 선배님 그래프의 기준값 (N, 1/den, Nr, A_dmin, d_min) : 비교용 */
typedef struct { int N, den, Nr; unsigned long long A; unsigned long long dmin; } Ref;
static const Ref REF[] = {
    { 256, 2,  8,  1632ULL,  8}, { 256, 3,  1,    32ULL,  8},
    { 512, 2,  2,   192ULL,  8}, { 512, 3, 10,  3168ULL, 16},
    {1024, 2,  1,   128ULL,  8}, {1024, 3,  3,   448ULL, 16},
    {2048, 2, 18, 61824ULL, 16}, {2048, 3,  1,   128ULL, 16},
};

/* ---------- d_min, A_dmin, 인덱스별 개수 ---------- */
static void run(int N, int num, int den, int crc, int method, double ebn0_db, int verbose)
{
    int n = 0; while ((1 << n) < N) n++;
    int K = (int)floor((double)N * num / den + 0.5) + crc;   /* 정보 비트 + CRC */
    int *inA = (int *)calloc(N, sizeof(int));
    construct(N, K, method, ebn0_db, inA);

    /* d_min = 2^(정보 행들의 최소 popcount) */
    int wmin = n + 1;
    for (int i = 0; i < N; i++)
        if (inA[i]) { int w = __builtin_popcount(i); if (w < wmin) wmin = w; }
    unsigned long long dmin = 1ULL << wmin;

    /* 이 개수 공식이 맞으려면 정보 집합이 decreasing(상위 집합으로 닫혀 있음)이어야 한다 */
    int viol = 0;
    for (int i = 0; i < N; i++) {
        if (!inA[i]) continue;
        for (int l = 0; l < n; l++) {
            if ((i >> l) & 1) continue;
            if (!inA[i | (1 << l)]) viol++;                      /* 0 -> 1 */
            for (int m = 0; m < l; m++)
                if (((i >> m) & 1) && !inA[i ^ (1 << m) ^ (1 << l)]) viol++;  /* 1을 위로 이동 */
        }
    }

    char fname[64];
    snprintf(fname, sizeof(fname), "dmin_N%d_R%d_%d.csv", N, num, den);
    FILE *fp = fopen(fname, "w");
    if (!fp) { perror("fopen"); free(inA); return; }
    fprintf(fp, "info_bit_index,num_min_weight_codewords,d_min\n");

    unsigned long long total = 0;
    int nidx = 0;
    for (int i = 0; i < N; i++) {
        if (!inA[i] || __builtin_popcount(i) != wmin) continue;
        /* 인덱스 i에서 시작하는 최소 weight 코드워드 수
         *   = 2^( (0인 비트 수) + (1이 0보다 낮은 자리에 있는 (m<l, i_m=1, i_l=0) 쌍의 수) ) */
        int inv = 0;
        for (int l = 0; l < n; l++)
            if (!((i >> l) & 1)) inv += __builtin_popcount(i & ((1 << l) - 1));
        unsigned long long cnt = 1ULL << ((n - wmin) + inv);
        total += cnt; nidx++;
        fprintf(fp, "%d,%llu,%llu\n", i + 1, cnt, dmin);        /* 1-based 인덱스 */
        if (verbose) printf("  index %4d : %llu\n", i + 1, cnt);
    }
    fclose(fp);

    /* 선배님 기준값과 비교 (beta, CRC=24, R=1/2 또는 1/3일 때) */
    const char *cmp = "";
    if (method == METHOD_BETA && crc == 24 && num == 1) {
        for (size_t t = 0; t < sizeof(REF) / sizeof(REF[0]); t++)
            if (REF[t].N == N && REF[t].den == den) {
                cmp = (REF[t].Nr == nidx && REF[t].A == total && REF[t].dmin == dmin)
                      ? "  [선배님 값과 일치]" : "  [선배님 값과 불일치]";
            }
    }

    printf("N=%4d R=%d/%d K=%4d | d_min=%4llu A_dmin=%llu (Nr=%d) -> %s%s%s\n",
           N, num, den, K, dmin, total, nidx, fname,
           viol ? "  [경고: 정보 집합이 decreasing이 아님, 개수 공식 검증 필요]" : "", cmp);
    free(inA);
}

int main(int argc, char **argv)
{
    int crc = 24;               /* CRC 비트 수 (K에 포함) */
    int method = METHOD_BETA;   /* 기본: beta-expansion */
    double ebn0 = 1.0;          /* GA 설계 Eb/N0 [dB] (method=ga일 때만 사용) */

    if (argc >= 4) {
        if (argc >= 5) crc = atoi(argv[4]);
        if (argc >= 6) method = (strcmp(argv[5], "ga") == 0) ? METHOD_GA : METHOD_BETA;
        if (argc >= 7) ebn0 = atof(argv[6]);
        run(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), crc, method, ebn0, 1);
        return 0;
    }
    int Ns[4] = {256, 512, 1024, 2048};
    int Rs[2][2] = {{1, 2}, {1, 3}};
    for (int a = 0; a < 2; a++)
        for (int b = 0; b < 4; b++)
            run(Ns[b], Rs[a][0], Rs[a][1], crc, method, ebn0, 0);
    return 0;
}