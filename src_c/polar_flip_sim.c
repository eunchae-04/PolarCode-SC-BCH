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
static uint8_t frozen[MAXN];        /* "현재 선택된 코드"의 frozen 표시 */
static int infoIdx[MAXN];           /* "현재 선택된 코드"의 정보 비트 u 인덱스 (오름차순) */

/* 코드는 두 종류: A = 일반 polar(K=메시지+CRC), B = BCH 연접 polar(K=BCH 부호어 총 길이) */
typedef struct { int K; uint8_t frozen[MAXN]; int infoIdx[MAXN]; } CodeDef;
static CodeDef codeA, codeB;

static void select_code(const CodeDef *cd)
{
    K = cd->K;
    memcpy(frozen, cd->frozen, N);
    memcpy(infoIdx, cd->infoIdx, sizeof(int) * cd->K);
}

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
static void construct_code(int Kc, CodeDef *cd)     /* 신뢰도 상위 Kc개 채널을 정보 비트로 */
{
    Rel *r = (Rel *)malloc(sizeof(Rel) * N);
    double beta = pow(2.0, 0.25);
    for (int i = 0; i < N; i++) {
        double pw = 0.0;
        for (int j = 0; j < n; j++) if ((i >> j) & 1) pw += pow(beta, j);
        r[i].idx = i; r[i].z = pw;
    }
    qsort(r, N, sizeof(Rel), cmp_rel);
    cd->K = Kc;
    for (int i = 0; i < N; i++) cd->frozen[i] = 1;
    for (int k = 0; k < Kc; k++) cd->frozen[r[k].idx] = 0;
    int c = 0;
    for (int i = 0; i < N; i++) if (!cd->frozen[i]) cd->infoIdx[c++] = i;
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

/* ===================== BCH (t=1, 단일 오류 정정) 외부 부호 ===================== */
/* 길이 2^m-1 인 원시 BCH 부호 중 t=1 인 것 = Hamming 부호. (7,4), (15,11), (31,26), (63,57), ...
 *  - 생성다항식 g(x) = 원시다항식 (차수 m = 중복 비트 r)
 *  - 마지막 블록은 짧은 부호(shortened)로 처리해서 불필요한 중복 비트를 보내지 않는다
 *  - 인코딩: 시스템적 [메시지 | 패리티], 디코딩: 신드롬 = 수신 다항식 mod g -> 표로 오류 위치 찾기
 *  - 인터리빙(il=1): 블록 b의 t번째 비트를 polar 정보 위치에 t*D+b 순으로 넣어서,
 *    SC가 연속된 정보 비트를 틀려도 블록 하나에는 오류가 1개만 들어가게 한다. (il=0: 블록을 차례로 이어 붙임) */
#define MAXBLK  512
#define MAXBLEN 256
typedef struct {
    int on, m, r, nb, kb, il;
    int P, D, KB;                    /* 페이로드(메시지+CRC) 비트 수, 블록 수, 연접 후 polar 정보 비트 수 */
    uint32_t glow;                   /* g(x)에서 x^r 항을 뺀 하위 r비트 */
    int klen[MAXBLK], len[MAXBLK];   /* 블록별 메시지 비트 수, 전송 비트 수(메시지+패리티) */
    int s_blk[MAXN], s_t[MAXN];      /* polar 정보 위치 s -> (블록, 블록 내 위치) */
    short tab[MAXBLK][MAXBLEN];      /* 블록별 신드롬 -> 오류 위치 (-1이면 없음) */
} Bch;
static Bch bch;
static const uint32_t PRIM_LOW[9] = {0, 0, 0, 0x3, 0x3, 0x5, 0x3, 0x9, 0x1D};   /* x^3+x+1, x^4+x+1, x^5+x^2+1, x^6+x+1, x^7+x^3+1, x^8+x^4+x^3+x^2+1 */

static int bch_reg(const uint8_t *b, int len)     /* LFSR 나눗셈 (CRC와 같은 방식) */
{
    uint32_t reg = 0, mask = (1u << bch.r) - 1;
    for (int i = 0; i < len; i++) {
        uint32_t fb = ((reg >> (bch.r - 1)) & 1u) ^ b[i];
        reg = (reg << 1) & mask;
        if (fb) reg ^= bch.glow;
    }
    return (int)reg;
}

static void bch_encode_block(const uint8_t *msg, int k, uint8_t *out)
{
    memcpy(out, msg, k);
    int reg = bch_reg(msg, k);
    for (int j = 0; j < bch.r; j++) out[k + j] = (reg >> (bch.r - 1 - j)) & 1;
}

static int bch_decode_block(uint8_t *y, int b)    /* 0: 오류 없음, 1: 1비트 정정, -1: 정정 불가 */
{
    int s = bch_reg(y, bch.len[b]);
    if (s == 0) return 0;
    int t = bch.tab[b][s];
    if (t < 0) return -1;
    y[t] ^= 1;
    return 1;
}

/* 반환: 0 성공, -1 불가(m 범위 밖, 연접 후 정보 비트 수가 N 초과 등) */
static int bch_setup(int P, int m, int il)
{
    memset(&bch, 0, sizeof(bch));
    if (m < 3 || m > 8) return -1;
    bch.m = m; bch.r = m; bch.nb = (1 << m) - 1; bch.kb = bch.nb - m;
    bch.il = il; bch.P = P; bch.glow = PRIM_LOW[m];
    bch.D = (P + bch.kb - 1) / bch.kb;
    if (bch.D > MAXBLK) return -1;
    int maxlen = 0;
    for (int b = 0; b < bch.D; b++) {
        bch.klen[b] = (b < bch.D - 1) ? bch.kb : P - (bch.D - 1) * bch.kb;
        bch.len[b] = bch.klen[b] + bch.r;
        bch.KB += bch.len[b];
        if (bch.len[b] > maxlen) maxlen = bch.len[b];
    }
    if (bch.KB > N) { bch.on = 0; return -1; }
    int s = 0;
    if (il) {
        for (int t = 0; t < maxlen; t++)
            for (int b = 0; b < bch.D; b++)
                if (t < bch.len[b]) { bch.s_blk[s] = b; bch.s_t[s] = t; s++; }
    } else {
        for (int b = 0; b < bch.D; b++)
            for (int t = 0; t < bch.len[b]; t++) { bch.s_blk[s] = b; bch.s_t[s] = t; s++; }
    }
    uint8_t e[MAXBLEN];
    for (int b = 0; b < bch.D; b++) {
        for (int i = 0; i < MAXBLEN; i++) bch.tab[b][i] = -1;
        memset(e, 0, sizeof(e));
        for (int t = 0; t < bch.len[b]; t++) {         /* 단위 오류 벡터의 신드롬 = 그 위치의 신드롬 */
            e[t] = 1; bch.tab[b][bch_reg(e, bch.len[b])] = (short)t; e[t] = 0;
        }
    }
    bch.on = 1;
    return 0;
}

/* 페이로드 pl[0..P-1] -> polar 정보 비트 스트림 sB[0..KB-1] */
static void bch_pack(const uint8_t *pl, uint8_t *sB)
{
    static uint8_t cw[MAXBLK][MAXBLEN];
    for (int b = 0; b < bch.D; b++) bch_encode_block(pl + b * bch.kb, bch.klen[b], cw[b]);
    for (int s = 0; s < bch.KB; s++) sB[s] = cw[bch.s_blk[s]][bch.s_t[s]];
}

/* 수신 스트림 -> 정정된 페이로드. 반환: 정정 불가 블록 수, *ncorr: 정정한 블록 수 */
static int bch_unpack(const uint8_t *sB, uint8_t *pl, int *ncorr)
{
    static uint8_t cw[MAXBLK][MAXBLEN];
    for (int s = 0; s < bch.KB; s++) cw[bch.s_blk[s]][bch.s_t[s]] = sB[s];
    int unc = 0; *ncorr = 0;
    for (int b = 0; b < bch.D; b++) {
        int rc = bch_decode_block(cw[b], b);
        if (rc < 0) unc++; else if (rc > 0) (*ncorr)++;
    }
    for (int b = 0; b < bch.D; b++) memcpy(pl + b * bch.kb, cw[b], bch.klen[b]);
    return unc;
}

/* BCH 정정 없이 수신 스트림에서 페이로드만 꺼냄 (BCH 정정이 실제로 얼마나 도움이 되는지 비교용) */
static void bch_raw(const uint8_t *sB, uint8_t *pl)
{
    static uint8_t cw[MAXBLK][MAXBLEN];
    for (int s = 0; s < bch.KB; s++) cw[bch.s_blk[s]][bch.s_t[s]] = sB[s];
    for (int b = 0; b < bch.D; b++) memcpy(pl + b * bch.kb, cw[b], bch.klen[b]);
}

/* BCH 자체 검증. 반환: 문제 개수 */
static int bch_selftest(void)
{
    if (!bch.on) { printf("[BCH] 사용 안 함\n"); return 0; }
    int bad = 0;
    uint8_t *pl = (uint8_t *)malloc(bch.P), *pl2 = (uint8_t *)malloc(bch.P);
    uint8_t *sB = (uint8_t *)malloc(bch.KB), *sE = (uint8_t *)malloc(bch.KB);
    for (int q = 0; q < bch.P; q++) pl[q] = rng_next() & 1;
    bch_pack(pl, sB);
    int nc, unc = bch_unpack(sB, pl2, &nc);
    if (unc || nc || memcmp(pl, pl2, bch.P)) bad++;            /* 오류 없으면 그대로 */
    int fail1 = 0;
    for (int s = 0; s < bch.KB; s++) {                          /* 스트림 위치마다 1비트 오류 -> 정정되어야 함 */
        memcpy(sE, sB, bch.KB); sE[s] ^= 1;
        unc = bch_unpack(sE, pl2, &nc);
        if (unc || nc != 1 || memcmp(pl, pl2, bch.P)) fail1++;
    }
    static uint8_t seen[MAXBLK][MAXBLEN];
    memset(seen, 0, sizeof(seen));
    int dup = 0;
    for (int s = 0; s < bch.KB; s++) { if (seen[bch.s_blk[s]][bch.s_t[s]]++) dup++; }
    int cnt = 0; for (int b = 0; b < bch.D; b++) cnt += bch.len[b];
    if (cnt != bch.KB) dup++;
    printf("[BCH(%d,%d) t=1] 블록 %d개, 연접 후 polar K=%d (CRC 포함 페이로드 %d), 인터리브 %s\n",
           bch.nb, bch.kb, bch.D, bch.KB, bch.P, bch.il ? "ON" : "OFF");
    printf("[BCH] 1비트 오류 정정 실패 %d / %d, 인터리버 중복 %d\n", fail1, bch.KB, dup);
    bad += fail1 + dup;
    free(pl); free(pl2); free(sB); free(sE);
    return bad;
}

/* BCH 오류 주입 실험: 정보 비트 스트림에 오류를 e개 넣었을 때 메시지가 완전히 복구되는 비율
 *   - 연속 오류: SC 디코딩 순서로 인접한 e개 비트가 틀린 경우 (인터리브가 효과를 내는 상황)
 *   - 임의 위치 오류: 서로 다른 e개 위치가 무작위로 틀린 경우 */
static void bch_inject_test(int trials)
{
    if (!bch.on) { printf("BCH를 사용할 수 없는 설정입니다\n"); return; }
    uint8_t *pl = (uint8_t *)malloc(bch.P), *pl2 = (uint8_t *)malloc(bch.P);
    uint8_t *sB = (uint8_t *)malloc(bch.KB), *sE = (uint8_t *)malloc(bch.KB);
    uint8_t *mk = (uint8_t *)calloc(bch.KB, 1);
    const int es[] = {1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 32};
    printf("BCH(%d,%d) t=1, 블록 %d개, 인터리브 %s, polar 정보 비트 %d (메시지+CRC %d), 시행 %d회\n",
           bch.nb, bch.kb, bch.D, bch.il ? "ON" : "OFF", bch.KB, bch.P, trials);
    printf("%-8s %-24s %-24s\n", "오류 수", "연속 오류 복구율", "임의 위치 오류 복구율");
    for (size_t i = 0; i < sizeof(es) / sizeof(es[0]); i++) {
        int e = es[i];
        if (e > bch.KB) break;
        int okb = 0, okr = 0, nc;
        for (int t = 0; t < trials; t++) {
            for (int q = 0; q < bch.P; q++) pl[q] = rng_next() & 1;
            bch_pack(pl, sB);
            memcpy(sE, sB, bch.KB);                              /* 연속 e개 */
            int s0 = (int)(rng_next() % (uint64_t)(bch.KB - e + 1));
            for (int j = 0; j < e; j++) sE[s0 + j] ^= 1;
            bch_unpack(sE, pl2, &nc);
            if (memcmp(pl, pl2, bch.P) == 0) okb++;
            memcpy(sE, sB, bch.KB);                              /* 서로 다른 임의 위치 e개 */
            int chosen = 0;
            while (chosen < e) {
                int p = (int)(rng_next() % (uint64_t)bch.KB);
                if (!mk[p]) { mk[p] = 1; sE[p] ^= 1; chosen++; }
            }
            memset(mk, 0, bch.KB);
            bch_unpack(sE, pl2, &nc);
            if (memcmp(pl, pl2, bch.P) == 0) okr++;
        }
        printf("%-8d %-24.1f %-24.1f\n", e, 100.0 * okb / trials, 100.0 * okr / trials);
    }
    printf("(%% 단위. 100%%면 메시지가 항상 완전히 복구됨)\n");
    free(pl); free(pl2); free(sB); free(sE); free(mk);
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
static int     sc_fpos = -1;      /* >=0 이면 이 u 인덱스의 결정을 강제로 sc_fval 로 둔다 (SC-Flip용) */
static uint8_t sc_fval = 0;

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
        uint8_t u = (base == sc_fpos) ? sc_fval : (frozen[base] ? 0 : (l < 0));
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
        {   /* d_min 오류 하나가 u-도메인(정보 비트)에서는 몇 개의 비트가 틀린 것으로 나타나는지 */
            int wmn = 1 << 30, wmx = 0; double wsum = 0;
            for (int q = 0; q < pat.A; q++) {
                int w = 0;
                for (int t = 0; t < pat.nw; t++) w += __builtin_popcountll(pat.bits[(size_t)q * pat.nw + t]);
                if (w < wmn) wmn = w;
                if (w > wmx) wmx = w;
                wsum += w;
            }
            printf("[패턴] d_min 오류 1개 = 틀린 정보 비트 수: 최소 %d, 평균 %.1f, 최대 %d\n", wmn, wsum / pat.A, wmx);
        }
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

/* ===================== SC-Flip: 첫 번째 틀린 비트를 고치고 SC를 이어서 다시 ===================== */
/* 아이디어: SC는 앞 비트 결정이 뒤 비트에 전부 영향을 주므로, "첫 번째 틀린 비트"만 정확히 고치면
 *           그 뒤는 SC가 알아서 올바르게 복구될 수 있다. 후보 위치 k를 하나 골라
 *           그 비트만 반대로 강제하고 SC를 처음부터 다시 돌린 뒤 CRC로 확인한다.
 *   mode 0 : 후보 = 시작 인덱스 집합 (popcount(index) == log2 d_min 인 정보 비트, Nr개)  <- 교수님 그래프의 x축
 *   mode 1 : 후보 = 1차 SC에서 |LLR|이 가장 작은 정보 비트 T개 (일반적인 SC-Flip, 비교 기준)
 * 후보는 |LLR| 오름차순(믿음이 약한 것부터)으로 시도. 반환: 0=CRC 이미 통과, 1=고쳐서 CRC 통과, -1=실패
 * *runs 에는 SC 실행 횟수(첫 번째 포함)를 돌려준다. */
static int cand_k[MAXN];
static double cand_m[MAXN];
static int cmp_cand(const void *a, const void *b)
{
    double x = cand_m[*(const int *)a], y = cand_m[*(const int *)b];
    return (x > y) - (x < y);
}

static int extraW = 0;   /* 시작 인덱스 집합을 popcount <= log2(d_min)+extraW 까지 넓힘 (0이면 교수님 그래프의 인덱스만) */
static int in_start_set(int k) { return __builtin_popcount(infoIdx[k]) <= pat.dmin_log + extraW; }

/* 후보 목록 만들기: 반환 개수 */
static int build_cands(const double *Ld0, int mode, int T)
{
    static int all[MAXN];
    int na = 0;
    for (int k = 0; k < K; k++) {
        if (mode == 0 && !in_start_set(k)) continue;
        cand_m[k] = fabs(Ld0[infoIdx[k]]);
        all[na++] = k;
    }
    qsort(all, na, sizeof(int), cmp_cand);
    int nc = (T < na) ? T : na;
    for (int i = 0; i < nc; i++) cand_k[i] = all[i];
    return nc;
}

static int scflip_decode(const double *llr, const uint8_t *uh0, const double *Ld0,
                         int mode, int T, uint8_t *out, int *runs)
{
    memcpy(out, uh0, K);
    *runs = 1;
    if (crc_rem(uh0, K) == 0) return 0;
    int nc = build_cands(Ld0, mode, T);
    static uint8_t tmp[MAXN];
    for (int i = 0; i < nc; i++) {
        int k = cand_k[i];
        sc_fpos = infoIdx[k]; sc_fval = uh0[k] ^ 1;
        sc_decode(llr);
        (*runs)++;
        for (int q = 0; q < K; q++) tmp[q] = sc_uhat[infoIdx[q]];
        if (crc_rem(tmp, K) == 0) { memcpy(out, tmp, K); sc_fpos = -1; return 1; }
    }
    sc_fpos = -1;
    return -1;
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

    /* 강제 결정: 정해진 위치의 결정만 바뀌고, 그 앞 비트는 그대로인지 */
    {
        int bad_f = 0, tr = 200;
        double *l2 = (double *)calloc(N, sizeof(double));
        uint8_t *ref = (uint8_t *)malloc(N);
        for (int t = 0; t < tr; t++) {
            for (int j = 0; j < N; j++) l2[j] = 3.0 * rng_gauss() + 0.5;
            sc_decode(l2); memcpy(ref, sc_uhat, N);
            int k = rng_next() % K, pos = infoIdx[k];
            sc_fpos = pos; sc_fval = ref[pos] ^ 1;
            sc_decode(l2); sc_fpos = -1;
            if (sc_uhat[pos] != (ref[pos] ^ 1)) bad_f++;
            for (int j = 0; j < pos; j++) if (sc_uhat[j] != ref[j]) { bad_f++; break; }
        }
        printf("[강제 결정 SC] %d회 중 이상 %d회\n", tr, bad_f);
        bad += bad_f; free(l2); free(ref);
    }

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
    bad += bch_selftest();
    printf("== 결과: %s ==\n", bad ? "문제 있음" : "정상");
    return bad;
}

/* ===================== 한 케이스 설정 / 해제 ===================== */
static int setup_case(int Nval, int num, int den, int Lval, int bch_m, int il)
{
    N = Nval; n = 0; while ((1 << n) < N) n++;
    if ((1 << n) != N || N > MAXN || Lval < 1 || Lval > MAXL) {
        fprintf(stderr, "N은 2의 거듭제곱(<=%d), 1<=L<=%d\n", MAXN, MAXL);
        return 1;
    }
    Lmax = Lval;
    Kmsg = (int)floor((double)N * num / den + 0.5);
    K = Kmsg + CRC_LEN;
    construct_code(K, &codeA);
    select_code(&codeA);
    bch.on = 0;
    if (bch_m > 0) {
        if (bch_setup(Kmsg + CRC_LEN, bch_m, il) == 0) construct_code(bch.KB, &codeB);
        else fprintf(stderr, "[BCH] m=%d: 이 N, R에서는 적용할 수 없어 BCH 비교를 생략합니다 (연접 후 정보 비트가 N을 넘음)\n", bch_m);
    }
    for (int j = 0; j < N; j++) {
        int r = 0, x = j;
        for (int b = 0; b < n; b++) { r = (r << 1) | (x & 1); x >>= 1; }
        bitrev_tab[j] = r;
    }
    sc_alloc(); scl_alloc();
    return 0;
}

static void teardown_case(void)
{
    for (int d = 0; d < n; d++) {
        free(scL[d]); free(scXl[d]); free(scXr[d]);
        scL[d] = NULL; scXl[d] = NULL; scXr[d] = NULL;
    }
    for (int d = 0; d <= n; d++)
        for (int s = 0; s < MAXL; s++)
            if (arrP[d][s]) { free(arrP[d][s]); free(arrC[d][s]); arrP[d][s] = NULL; arrC[d][s] = NULL; }
    free(pat.bits); free(pat.syn); free(pat.start);
    memset(&pat, 0, sizeof(pat));
}

/* ===================== 시뮬레이션 (한 케이스, SNR 구간 전체) ===================== */
/* stop_on: 어느 디코더의 프레임 오류가 targ개가 되면 그 SNR을 끝낼지 (0: SC, 1: SC+뒤집기, 2: SCL) */
static int run_sim(int num, int den, double e0, double e1, double es, long maxf, int targ, int stop_on, int Tflip)
{
    static const char *stop_name[3] = {"SC", "SC+flip", "SCL"};
    printf("N=%d R=%d/%d K=%d (메시지 %d + CRC %d), SCL L=%d\n", N, num, den, K, Kmsg, CRC_LEN, Lmax);
    printf("종료 조건: %s 프레임 오류 %d개 도달 또는 %ld프레임\n", stop_name[stop_on], targ, maxf);
    if (bch.on)
        printf("SC+BCH: BCH(%d,%d) t=1, 블록 %d개, 인터리브 %s, polar K %d -> %d\n",
               bch.nb, bch.kb, bch.D, bch.il ? "ON" : "OFF", K, bch.KB);
    else
        printf("SC+BCH: 사용 안 함\n");
    if (build_patterns(1)) { fprintf(stderr, "패턴 생성 오류\n"); return 2; }

    char fname[128];
    snprintf(fname, sizeof(fname), "../results/perf_N%d_R%d_%d.csv", N, num, den);
    FILE *fp = fopen(fname, "w");
    if (!fp) { snprintf(fname, sizeof(fname), "perf_N%d_R%d_%d.csv", N, num, den); fp = fopen(fname, "w"); }
    if (!fp) { perror("fopen"); return 1; }
    fprintf(fp, "ebn0_db,frames,sc_fer,sc_ber,bch_fer,bch_ber,flip_fer,flip_ber,scl_fer,scl_ber,"
                "sc_err,bch_err,flip_err,scl_err,sc_err_crc_detected,flip_fixed,flip_miscorrected,flip_nomatch,"
                "wrong1,wrong2,wrong3_4,wrong5_8,wrong9p,bch_raw_err,"
                "sfs_fer,sfs_ber,sfl_fer,sfl_ber,sfs_runs,sfl_runs,first_total,first_in_start,first_in_llr,T\n");

    double Reff = (double)Kmsg / N;
    uint8_t *uinfo = (uint8_t *)malloc(K), *uh = (uint8_t *)malloc(K), *uh0 = (uint8_t *)malloc(K), *uh1 = (uint8_t *)malloc(K);
    double *Ld0 = (double *)malloc(sizeof(double) * N);
    uint8_t *u = (uint8_t *)malloc(N), *x = (uint8_t *)malloc(N);
    double *llr = (double *)malloc(sizeof(double) * N);
    uint8_t *sB = (uint8_t *)malloc(N), *sBh = (uint8_t *)malloc(N), *plh = (uint8_t *)malloc(K + 8);

    printf("%-6s %-9s | %-9s | %-9s | %-9s | %-9s | (프레임 오류 수: SC BCH FLIP SCL)\n",
           "Eb/N0", "frames", "SC FER", "BCH FER", "FLIP FER", "SCL FER");
    for (double eb = e0; eb <= e1 + 1e-9; eb += es) {
        double sigma2 = 1.0 / (2.0 * Reff * pow(10.0, eb / 10.0)), sigma = sqrt(sigma2);
        long frames = 0;
        long sc_err = 0, flip_err = 0, scl_err = 0, sc_bit = 0, flip_bit = 0, scl_bit = 0;
        long bch_err = 0, bch_bit = 0, bch_raw_err = 0;
        long det = 0, fixed = 0, miscorr = 0, nomatch = 0;
        long sfs_err = 0, sfs_bit = 0, sfl_err = 0, sfl_bit = 0, sfs_runs = 0, sfl_runs = 0;
        long first_tot = 0, first_start = 0, first_llr = 0;
        long hw[5] = {0, 0, 0, 0, 0};      /* SC가 틀린 프레임에서 틀린 정보 비트 수 분포: 1, 2, 3-4, 5-8, 9+ */
        while (frames < maxf &&
               (stop_on == 0 ? sc_err : stop_on == 1 ? flip_err : scl_err) < targ) {
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
            {   /* 틀린 정보 비트(CRC 포함) 개수 분포: BCH t=1 정정이 가능한 오류 모양인지 보기 위함 */
                int nw = 0; for (int k = 0; k < K; k++) nw += (uh[k] != uinfo[k]);
                if (nw) hw[nw == 1 ? 0 : nw == 2 ? 1 : nw <= 4 ? 2 : nw <= 8 ? 3 : 4]++;
            }

            /* 1차 SC 결과 보존 (이후 디코더들이 덮어쓰기 때문) */
            memcpy(uh0, uh, K);
            for (int j = 0; j < N; j++) Ld0[j] = sc_Ldec[j];
            {   /* 첫 번째로 틀린 정보 비트가 (a) 시작 인덱스 집합 (b) |LLR| 하위 T개 안에 있는가 */
                int fk = -1;
                for (int k = 0; k < K; k++) if (uh0[k] != uinfo[k]) { fk = k; break; }
                if (fk >= 0) {
                    for (int k = fk + 1; k < K; k++)       /* infoIdx 가 오름차순이 아닐 경우 대비: u 인덱스가 가장 작은 것 */
                        if (uh0[k] != uinfo[k] && infoIdx[k] < infoIdx[fk]) fk = k;
                    first_tot++;
                    if (in_start_set(fk)) first_start++;
                    int nc = build_cands(Ld0, 1, Tflip);
                    for (int i = 0; i < nc; i++) if (cand_k[i] == fk) { first_llr++; break; }
                }
            }

            /* --- SC-Flip: 시작 인덱스 후보 / |LLR| 하위 T개 후보 (첫 오류 비트만 고치고 SC 재실행) --- */
            {
                int rn, e;
                scflip_decode(llr, uh0, Ld0, 0, Tflip, uh1, &rn); sfs_runs += rn;
                e = 0; for (int k = 0; k < Kmsg; k++) e += (uh1[k] != uinfo[k]);
                sfs_bit += e; if (e) sfs_err++;
                scflip_decode(llr, uh0, Ld0, 1, Tflip, uh1, &rn); sfl_runs += rn;
                e = 0; for (int k = 0; k < Kmsg; k++) e += (uh1[k] != uinfo[k]);
                sfl_bit += e; if (e) sfl_err++;
            }

            /* --- SC + 패턴 뒤집기(신드롬 매칭) --- */
            int rc = flip_decode(uh, Ld0);
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

            /* --- SC + BCH (별도의 BCH 연접 polar 코드 B, 같은 메시지, 새 잡음) --- */
            if (bch.on) {
                bch_pack(uinfo, sB);                          /* 메시지+CRC -> BCH 블록 -> 인터리브 */
                select_code(&codeB);
                memset(u, 0, N);
                for (int s = 0; s < bch.KB; s++) u[infoIdx[s]] = sB[s];
                memcpy(x, u, N); polar_transform(x);
                for (int j = 0; j < N; j++) llr[j] = 2.0 * ((1.0 - 2.0 * x[j]) + sigma * rng_gauss()) / sigma2;
                sc_decode(llr);
                for (int s = 0; s < bch.KB; s++) sBh[s] = sc_uhat[infoIdx[s]];
                int ncorr;
                bch_raw(sBh, plh);                            /* 정정 전 (같은 SC 결과) */
                { int e_r = 0; for (int q = 0; q < Kmsg; q++) e_r += (plh[q] != uinfo[q]); if (e_r) bch_raw_err++; }
                bch_unpack(sBh, plh, &ncorr);                 /* BCH 디코딩(블록당 1비트 정정) */
                int e_b = 0; for (int q = 0; q < Kmsg; q++) e_b += (plh[q] != uinfo[q]);
                bch_bit += e_b; if (e_b) bch_err++;
                select_code(&codeA);
            }
        }
        double fn = (double)frames, bn = (double)frames * Kmsg;
        double bfer = bch.on ? bch_err / fn : -1.0, bber = bch.on ? bch_bit / bn : -1.0;
        printf("%-6.2f %-9ld | %-9.2e | %-9.2e | %-9.2e | %-9.2e | %ld %ld %ld %ld\n",
               eb, frames, sc_err / fn, bfer, flip_err / fn, scl_err / fn,
               sc_err, bch.on ? bch_err : -1L, flip_err, scl_err);
        printf("       SC 오류 프레임의 틀린 정보 비트 수: 1비트 %ld, 2비트 %ld, 3~4비트 %ld, 5~8비트 %ld, 9비트 이상 %ld\n",
               hw[0], hw[1], hw[2], hw[3], hw[4]);
        if (bch.on)
            printf("       BCH: 정정 전 프레임 오류 %ld -> 정정 후 %ld (정정이 줄인 오류 %ld)\n",
                   bch_raw_err, bch_err, bch_raw_err - bch_err);
        printf("       SC-Flip(시작 인덱스 %d개): FER %.2e, 평균 SC 실행 %.2f회 | SC-Flip(|LLR| 하위 T=%d): FER %.2e, 평균 %.2f회\n",
               pat.Nr, sfs_err / fn, sfs_runs / fn, Tflip, sfl_err / fn, sfl_runs / fn);
        printf("       SC의 첫 오류 비트 위치: 전체 %ld개 중 시작 인덱스 집합 %ld개(%.1f%%), |LLR| 하위 T 안 %ld개(%.1f%%)\n",
               first_tot, first_start, first_tot ? 100.0 * first_start / first_tot : 0.0,
               first_llr, first_tot ? 100.0 * first_llr / first_tot : 0.0);
        printf("       뒤집기: SC 오류 중 CRC 검출 %ld, 복구 %ld, 오정정 %ld, 매칭 패턴 없음 %ld\n", det, fixed, miscorr, nomatch);
        fprintf(fp, "%.2f,%ld,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%.6e,%.6e,%.6e,%.6e,%.3f,%.3f,%ld,%ld,%ld,%d\n",
                eb, frames, sc_err / fn, sc_bit / bn, bfer, bber, flip_err / fn, flip_bit / bn, scl_err / fn, scl_bit / bn,
                sc_err, bch.on ? bch_err : -1L, flip_err, scl_err, det, fixed, miscorr, nomatch,
                hw[0], hw[1], hw[2], hw[3], hw[4], bch.on ? bch_raw_err : -1L,
                sfs_err / fn, sfs_bit / bn, sfl_err / fn, sfl_bit / bn, sfs_runs / fn, sfl_runs / fn,
                first_tot, first_start, first_llr, Tflip);
        fflush(fp); fflush(stdout);
    }
    fclose(fp);
    free(uinfo); free(uh); free(uh0); free(uh1); free(Ld0); free(u); free(x); free(llr); free(sB); free(sBh); free(plh);
    printf("결과 저장: %s\n\n", fname);
    return 0;
}

/* ===================== 메인 ===================== */
static void usage(const char *prog)
{
    fprintf(stderr,
        "사용법:\n"
        "  %s N num den [L] [EbN0_start] [EbN0_stop] [step] [max_frames] [target_err] [seed] [stop_on] [bch_m] [il] [T] [W]   (한 케이스)\n"
        "  %s all [L] [EbN0_start] [EbN0_stop] [step] [max_frames] [target_err] [seed] [stop_on] [bch_m] [il] [T]         (N=256,512,1024,2048 x R=1/2,1/3 전부)\n"
        "  %s check N num den [L] ... [bch_m] [il]                                                                    (자체 검증)\n"
        "  %s bchtest N num den [bch_m] [il] [trials] [seed]            (BCH 오류 주입 실험: 오류 1개, 연속/임의 e개 복구율)\n"
        "  stop_on: scl(기본) | flip | sc  -- 이 디코더의 프레임 오류가 target_err개가 되면 그 SNR 종료\n"
        "  bch_m  : BCH(2^m-1, 2^m-1-m) t=1 외부 부호. 3~8, 0이면 SC+BCH 생략 (기본 0 = 생략, 4 = BCH(15,11))\n"
        "  W      : 시작 인덱스 후보를 popcount <= log2(d_min)+W 인 정보 비트까지 넓힘 (기본 0)\n"
        "  T      : 일반 SC-Flip 후보 수(|LLR| 가장 작은 정보 비트 T개, 기본 10)\n"
        "  il     : 1이면 BCH 블록을 인터리브(기본), 0이면 블록을 차례로 이어 붙임\n",
        prog, prog, prog, prog);
}

int main(int argc, char **argv)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001);      /* 콘솔 출력을 UTF-8로 */
#endif
    if (argc < 2) { usage(argv[0]); return 1; }

    int mode = 0;                   /* 0: 한 케이스, 1: all, 2: check */
    int argi = 1;
    if (strcmp(argv[1], "all") == 0) { mode = 1; argi = 2; }
    else if (strcmp(argv[1], "check") == 0) { mode = 2; argi = 2; }
    else if (strcmp(argv[1], "bchtest") == 0) {
        /* ./polar_flip_sim bchtest N num den [bch_m=4] [il=1] [trials=5000] [seed=1] : BCH 오류 주입 실험 (polar 디코딩 없음) */
        if (argc < 5) { usage(argv[0]); return 1; }
        int m = (argc > 5) ? atoi(argv[5]) : 4, ilv = (argc > 6) ? atoi(argv[6]) : 1;
        int trials = (argc > 7) ? atoi(argv[7]) : 5000;
        uint64_t sd = (argc > 8) ? (uint64_t)atoll(argv[8]) : 1;
        rng_s ^= sd * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 20; i++) rng_next();
        if (setup_case(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), 8, m, ilv)) return 1;
        bch_inject_test(trials);
        teardown_case();
        return 0;
    }

    int Nval = 0, num = 1, den = 2;
    if (mode != 1) {
        if (argc - argi < 3) { usage(argv[0]); return 1; }
        Nval = atoi(argv[argi]); num = atoi(argv[argi + 1]); den = atoi(argv[argi + 2]);
        argi += 3;
    }
    /* 나머지 인자 (mode 2는 사용하지 않음) */
    int Lval  = (argc > argi)     ? atoi(argv[argi])     : 8;
    double e0 = (argc > argi + 1) ? atof(argv[argi + 1]) : (mode == 1 ? 1.0 : 1.0);
    double e1 = (argc > argi + 2) ? atof(argv[argi + 2]) : (mode == 1 ? 5.0 : 3.0);
    double es = (argc > argi + 3) ? atof(argv[argi + 3]) : 0.5;
    long maxf = (argc > argi + 4) ? atol(argv[argi + 4]) : (mode == 1 ? 20000 : 100000);
    int  targ = (argc > argi + 5) ? atoi(argv[argi + 5]) : 100;
    uint64_t seed = (argc > argi + 6) ? (uint64_t)atoll(argv[argi + 6]) : 1;
    int stop_on = 2;                /* 기본: SCL (오류가 가장 드물어서 SCL 점까지 확보) */
    if (argc > argi + 7) {
        if (strcmp(argv[argi + 7], "sc") == 0) stop_on = 0;
        else if (strcmp(argv[argi + 7], "flip") == 0) stop_on = 1;
        else if (strcmp(argv[argi + 7], "scl") == 0) stop_on = 2;
        else { usage(argv[0]); return 1; }
    }
    int bch_m = (argc > argi + 8) ? atoi(argv[argi + 8]) : 0;
    int il    = (argc > argi + 9) ? atoi(argv[argi + 9]) : 1;
    int Tflip = (argc > argi + 10) ? atoi(argv[argi + 10]) : 10;
    extraW = (argc > argi + 11) ? atoi(argv[argi + 11]) : 0;
    rng_s ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 20; i++) rng_next();

    if (mode == 1) {
        int Ns[4] = {256, 512, 1024, 2048};
        int Rs[2][2] = {{1, 2}, {1, 3}};
        for (int a = 0; a < 2; a++) {
            for (int b = 0; b < 4; b++) {
                printf("==== 케이스 %d/8 ====\n", a * 4 + b + 1);
                if (setup_case(Ns[b], Rs[a][0], Rs[a][1], Lval, bch_m, il)) return 1;
                int rc = run_sim(Rs[a][0], Rs[a][1], e0, e1, es, maxf, targ, stop_on, Tflip);
                teardown_case();
                if (rc) return rc;
            }
        }
        return 0;
    }

    if (setup_case(Nval, num, den, Lval, bch_m, il)) return 1;
    if (mode == 2) return self_check() ? 2 : 0;
    int rc = run_sim(num, den, e0, e1, es, maxf, targ, stop_on, Tflip);
    teardown_case();
    return rc;
}
