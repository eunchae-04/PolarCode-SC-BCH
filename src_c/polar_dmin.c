/*
 * polar_dmin.c
 * Polar code의 d_min, A_dmin(최소 weight 코드워드 개수)와
 * "정보 비트 인덱스별 최소 weight 코드워드 개수"를 계산해서 CSV로 저장한다.
 * (그래프는 plot_dmin.m 으로 MATLAB에서 그린다)
 *
 * 컴파일:  gcc -O2 -o polar_dmin polar_dmin.c -lm
 * 사용법:
 *   ./polar_dmin                     -> N=256,512,1024,2048 x R=1/2,1/3 전부 실행
 *   ./polar_dmin N num den [CRC] [EbN0dB] -> 하나만 실행 (예: ./polar_dmin 256 1 3 24 1.0)
 *
 * 가정 / 약속
 *  - G = F^{(x)n}, F=[1 0;1 1] (natural order). 행 i의 weight = 2^popcount(i)
 *  - frozen set: Gaussian Approximation, 설계 Eb/N0 (기본: R=1/2는 2.0 dB, R=1/3은 1.0 dB, 설계 rate=K/N), 상위 K개 채널이 정보 비트
 *  - K = round(N*R) + CRC (기본 CRC=24). 선배님 그래프(N=256, R=1/3 -> 인덱스 225, A=32)와 맞추려면 CRC를 K에 포함해야 함
 *  - 출력 인덱스는 1-based (MATLAB 그래프와 맞추기 위해)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI_CONST 3.14159265358979323846

typedef struct { int idx; double z; } Rel;

/* ---------- Gaussian Approximation ---------- */
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

static int cmp_rel(const void *a, const void *b)   /* z 내림차순, 동점이면 큰 인덱스 먼저 */
{
    const Rel *x = (const Rel *)a, *y = (const Rel *)b;
    if (x->z > y->z) return -1;
    if (x->z < y->z) return 1;
    return y->idx - x->idx;
}

static void construct(int N, int K, double R, double ebn0_db, int *inA)
{
    int n = 0; while ((1 << n) < N) n++;
    double *z = (double *)malloc(sizeof(double) * N);
    Rel *r = (Rel *)malloc(sizeof(Rel) * N);

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
    qsort(r, N, sizeof(Rel), cmp_rel);
    for (int i = 0; i < N; i++) inA[i] = 0;
    for (int k = 0; k < K; k++) inA[r[k].idx] = 1;
    free(z); free(r);
}

/* ---------- d_min, A_dmin, 인덱스별 개수 ---------- */
static void run(int N, int num, int den, int crc, double ebn0_db, int verbose)
{
    int n = 0; while ((1 << n) < N) n++;
    int K = (int)floor((double)N * num / den + 0.5) + crc;   /* 정보 비트 + CRC */
    int *inA = (int *)calloc(N, sizeof(int));
    construct(N, K, (double)K / N, ebn0_db, inA);   /* 설계 rate = (정보+CRC)/N */

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

    printf("N=%4d R=%d/%d K=%4d | d_min=%4llu A_dmin=%llu (#index=%d) -> %s%s\n",
           N, num, den, K, dmin, total, nidx, fname,
           viol ? "  [경고: 정보 집합이 decreasing이 아님, 개수 공식 검증 필요]" : "");
    free(inA);
}

int main(int argc, char **argv)
{
    double ebn0 = -1.0;  /* 설계 Eb/N0 [dB]; 음수면 자동(R=1/2 -> 2.0, 그 외 1.0) */
    int crc = 24;        /* CRC 비트 수 (K에 포함) */
    if (argc >= 4) {
        if (argc >= 5) crc  = atoi(argv[4]);
        if (argc >= 6) ebn0 = atof(argv[5]);
        if (ebn0 < 0) ebn0 = (atoi(argv[3]) == 2) ? 2.0 : 1.0;
        run(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), crc, ebn0, 1);
        return 0;
    }
    int Ns[4] = {256, 512, 1024, 2048};
    int Rs[2][2] = {{1, 2}, {1, 3}};
    for (int a = 0; a < 2; a++)
        for (int b = 0; b < 4; b++)
            run(Ns[b], Rs[a][0], Rs[a][1], crc, (ebn0 >= 0) ? ebn0 : (Rs[a][1] == 2 ? 2.0 : 1.0), 0);
    return 0;
}
