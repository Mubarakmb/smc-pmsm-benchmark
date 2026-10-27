#define _POSIX_C_SOURCE 199309L
/* Per-sample execution time of the six speed controllers (same laws as sim.py).
 * Each step: read e, w, dwr -> compute iq* (incl. saturation + anti-windup + state update).
 * Compile: gcc -O2 -std=c11 ctrl_bench.c -lm -o ctrl_bench
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define LM 2000
static inline double sgn(double x){ return (x>0)-(x<0); }
static inline double clip(double u,double m){ return u>m?m:(u<-m?-m:u); }

typedef struct { double chi, ah, Ts, Imax; } Com;
typedef struct { double I1, e0, K, z, v; double ebuf[LM]; double wa[LM], wb[LM]; int bidx, nb; } St;

static double csmc(const Com*c, St*s, const double*G, double e, double w, double dwr){
    double N = G[0]*sgn(e) + G[1]*e;
    return clip((dwr + c->ah*w + N)/c->chi, c->Imax); }

static double ismc(const Com*c, St*s, const double*G, double e, double w, double dwr){
    double sv = e - s->e0 + G[0]*s->I1;
    double N = G[0]*e + G[1]*sgn(sv) + G[3]*sv;
    double u = (dwr + c->ah*w + N)/c->chi, uc = clip(u,c->Imax);
    if (!(uc!=u && sgn(u)==sgn(e))) s->I1 += e*c->Ts;
    return uc; }

static double tsmc(const Com*c, St*s, const double*G, double e, double w, double dwr){
    double ft = G[0]*pow(fabs(e),G[1])*sgn(e);
    double sv = e + s->I1;
    double N = ft + G[3]*sgn(sv) + G[4]*sv;
    double u = (dwr + c->ah*w + N)/c->chi, uc = clip(u,c->Imax);
    if (!(uc!=u && sgn(u)==sgn(e))) s->I1 += ft*c->Ts;
    return uc; }

static double asmc(const Com*c, St*s, const double*G, double e, double w, double dwr){
    double sv = e, sw = sgn(sv);
    double N = s->K*sw + G[1]*sv;
    double u = (dwr + c->ah*w + N)/c->chi, uc = clip(u,c->Imax);
    s->z += c->Ts/G[3]*(sw - s->z);
    double dl = fabs(s->z) - G[4];
    double dK = (G[5]+G[6]*fabs(e))*s->K*sgn(dl) - G[7]*fmax(s->K-G[9],0.0) + G[7]*fmax(G[8]-s->K,0.0);
    s->K += dK*c->Ts;
    return uc; }

static double fosmc(const Com*c, St*s, const double*G, double e, double w, double dwr){
    s->ebuf[s->bidx] = e; if (s->nb < LM) s->nb++;
    double A=0.0, B=0.0;
    for (int j=0;j<s->nb;j++){ int k = s->bidx-j; if (k<0) k+=LM; A += s->wa[j]*s->ebuf[k]; B += s->wb[j]*s->ebuf[k]; }
    double al = G[1];
    double Dm = A*pow(c->Ts,al), D1 = B*pow(c->Ts,al-1.0);
    s->bidx = (s->bidx+1)%LM;
    double sv = e + G[0]*Dm;
    double N = G[0]*D1 + G[2]*sgn(sv) + G[3]*sv;
    return clip((dwr + c->ah*w + N)/c->chi, c->Imax); }

static double stsmc(const Com*c, St*s, const double*G, double e, double w, double dwr){
    double sv = e;
    double N = G[0]*sqrt(fabs(sv))*sgn(sv) + s->v;
    double u = (dwr + c->ah*w + N)/c->chi, uc = clip(u,c->Imax);
    if (!(uc!=u && sgn(u)==sgn(e))) s->v += G[1]*sgn(sv)*c->Ts;
    return uc; }

typedef double (*Fn)(const Com*, St*, const double*, double, double, double);

static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec + 1e-9*t.tv_nsec; }

int main(int argc, char**argv){
    long N = argc>1 ? atol(argv[1]) : 10000000L;
    Com c = { 1.5*4*0.1706/15.1e-4, 0.0, 1e-4, 39.6 };
    double G[6][12] = {
        {5800, 400}, {100, 5800, 0, 400}, {100, 0.7, 0, 5800, 400},
        {0, 400, 58, 2e-3, 0.8, 50, 0, 1e3, 58, 14500}, {100, 0.5, 5800, 400}, {1500, 1e6} };
    const char* nm[6] = {"CSMC","ISMC","TSMC","ASMC","FOSMC","STSMC"};
    Fn f[6] = {csmc, ismc, tsmc, asmc, fosmc, stsmc};
    St *s = calloc(1, sizeof(St));
    /* GL weights for FOSMC */
    double al = 0.5; s->wa[0]=1; s->wb[0]=1;
    for (int j=1;j<LM;j++){ s->wa[j]=s->wa[j-1]*(1-(-al+1)/j); s->wb[j]=s->wb[j-1]*(1-((1-al)+1)/j); }
    /* pseudo-random but deterministic inputs, precomputed so input generation is not timed */
    int NI = 4096; double *E = malloc(NI*sizeof(double)), *W = malloc(NI*sizeof(double));
    srand(1); for (int i=0;i<NI;i++){ E[i] = ((double)rand()/RAND_MAX-0.5)*2.0; W[i] = 300.0 + E[i]; }
    printf("controller,ns_per_step,N\n");
    for (int k=0;k<6;k++){
        St st = *s; volatile double sink = 0;
        long n = (k==4) ? N/50 : N;               /* FOSMC is much slower */
        for (long i=0;i<10000;i++) sink += f[k](&c,&st,G[k],E[i&(NI-1)],W[i&(NI-1)],0.0);  /* warm-up */
        double t0 = now();
        for (long i=0;i<n;i++) sink += f[k](&c,&st,G[k],E[i&(NI-1)],W[i&(NI-1)],0.0);
        double t1 = now();
        printf("%s,%.2f,%ld\n", nm[k], 1e9*(t1-t0)/n, n);
    }
    return 0;
}
