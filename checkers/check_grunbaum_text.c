/*
 * Exact checker for the four-hyperplane tensor certificate.
 *
 * Usage:
 *   cc -O3 -std=c11 check_grunbaum_text.c -o check_grunbaum_text
 *   ./check_grunbaum_text grunbaum_certificate.txt
 *
 * The checker:
 *   1. derives the five six-variable chart polynomials from the displayed
 *      symmetric tensors A and B and the quaternionic Spin(4) formulas;
 *   2. converts them exactly to degree-4 tensor-product Bernstein form;
 *   3. reads and replays the human-readable dyadic subdivision tree; and
 *   4. verifies every leaf separator using signed integer arithmetic.
 *
 * No floating-point arithmetic is used.  The certificate producer and its
 * search heuristics are not part of the trusted computation.  Compilation
 * requires GCC or Clang support for signed __int128; no third-party library
 * is required.
 */

#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef __SIZEOF_INT128__
#error "This checker requires a compiler with signed __int128 support (GCC or Clang)."
#endif

_Static_assert(sizeof(int64_t)==8, "This checker requires 64-bit int64_t.");

#define NV 6
#define DEG 4
#define SIDE 5
#define NCTRL 15625
#define NCOMP 5
#define LEAF_TOKEN 6
#define FIXED_SCALE ((int64_t)1 << 52)

static const int STRIDE[NV] = {3125, 625, 125, 25, 5, 1};
static const int64_t T6[SIDE][SIDE] = {
    { 6, -6,  6, -6,  6},
    { 6, -3,  0,  3, -6},
    { 6,  0, -2,  0,  6},
    { 6,  3,  0, -3, -6},
    { 6,  6,  6,  6,  6}
};

/* A sparse polynomial used only for the 16 bilinear frame entries. */
typedef struct {
    int n;
    uint16_t code[8];
    int8_t coeff[8];
} SPoly;

typedef struct {
    int64_t *lo; /* NCOMP contiguous control nets */
    int64_t *hi;
} Net;

typedef struct {
    uint8_t *tokens;
    uint64_t ntokens;
    int64_t *witnesses; /* 5 integers per leaf */
    uint64_t nleaves;
    uint64_t token_pos;
    uint64_t leaf_pos;
    uint64_t checked_nodes;
    uint64_t checked_leaves;
    int max_depth;
    __int128 min_leaf_sum;
} Certificate;

static void die(const char *message) {
    fprintf(stderr, "ERROR: %s\n", message);
    exit(1);
}

static void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (!p) die("out of memory");
    return p;
}

static uint16_t mono_var(int v) {
    return (uint16_t)STRIDE[v];
}

static uint16_t mono_prod(int v, int w) {
    return (uint16_t)(STRIDE[v] + STRIDE[w]);
}

static void sp_add(SPoly *p, int coeff, uint16_t code) {
    if (p->n >= 8) die("internal sparse-polynomial capacity exceeded");
    p->coeff[p->n] = (int8_t)coeff;
    p->code[p->n] = code;
    p->n++;
}

/*
 * Unnormalized columns q_j = p e_j conjugate(r), with
 * p=(1,x1,x2,x3), r=(1,y1,y2,y3).
 * q[row][column] is a sparse polynomial in the six chart variables.
 */
static void build_frame(SPoly q[4][4]) {
    memset(q, 0, 16 * sizeof(SPoly));
#define C(R,COL,A) sp_add(&q[(R)][(COL)], (A), 0)
#define V(R,COL,A,I) sp_add(&q[(R)][(COL)], (A), mono_var((I)))
#define P(R,COL,A,I,J) sp_add(&q[(R)][(COL)], (A), mono_prod((I),(J)))
    /* column 0 */
    C(0,0, 1); P(0,0, 1,0,3); P(0,0, 1,1,4); P(0,0, 1,2,5);
    V(1,0, 1,0); P(1,0,-1,1,5); P(1,0, 1,2,4); V(1,0,-1,3);
    P(2,0, 1,0,5); V(2,0, 1,1); P(2,0,-1,2,3); V(2,0,-1,4);
    P(3,0,-1,0,4); P(3,0, 1,1,3); V(3,0, 1,2); V(3,0,-1,5);
    /* column 1 */
    V(0,1,-1,0); P(0,1,-1,1,5); P(0,1, 1,2,4); V(0,1, 1,3);
    C(1,1, 1); P(1,1, 1,0,3); P(1,1,-1,1,4); P(1,1,-1,2,5);
    P(2,1, 1,0,4); P(2,1, 1,1,3); V(2,1, 1,2); V(2,1, 1,5);
    P(3,1, 1,0,5); V(3,1,-1,1); P(3,1, 1,2,3); V(3,1,-1,4);
    /* column 2 */
    P(0,2, 1,0,5); V(0,2,-1,1); P(0,2,-1,2,3); V(0,2, 1,4);
    P(1,2, 1,0,4); P(1,2, 1,1,3); V(1,2,-1,2); V(1,2,-1,5);
    C(2,2, 1); P(2,2,-1,0,3); P(2,2, 1,1,4); P(2,2,-1,2,5);
    V(3,2, 1,0); P(3,2, 1,1,5); P(3,2, 1,2,4); V(3,2, 1,3);
    /* column 3 */
    P(0,3,-1,0,4); P(0,3, 1,1,3); V(0,3,-1,2); V(0,3, 1,5);
    P(1,3, 1,0,5); V(1,3, 1,1); P(1,3, 1,2,3); V(1,3, 1,4);
    V(2,3,-1,0); P(2,3, 1,1,5); P(2,3, 1,2,4); V(2,3,-1,3);
    C(3,3, 1); P(3,3,-1,0,3); P(3,3,-1,1,4); P(3,3, 1,2,5);
#undef C
#undef V
#undef P
}

static void sort3(int *a, int *b, int *c) {
    int t;
    if (*a > *b) { t=*a; *a=*b; *b=t; }
    if (*b > *c) { t=*b; *b=*c; *c=t; }
    if (*a > *b) { t=*a; *a=*b; *b=t; }
}

static void sort4(int *a, int *b, int *c, int *d) {
    int x[4] = {*a,*b,*c,*d};
    for (int i=1;i<4;i++) {
        int v=x[i],j=i-1;
        while (j>=0 && x[j]>v) { x[j+1]=x[j]; j--; }
        x[j+1]=v;
    }
    *a=x[0]; *b=x[1]; *c=x[2]; *d=x[3];
}

/* Five times the cubic tensor A from the manuscript. */
static int64_t A5(int a, int b, int c) {
    sort3(&a,&b,&c);
    const int key = 16*a + 4*b + c;
    switch (key) {
        case  0: return -65; /* 000 */
        case  1: return  60; /* 001 */
        case  2: return -15; /* 002 */
        case  3: return -25; /* 003 */
        case  5: return  25; /* 011 */
        case 10: return -75; /* 022 */
        case 15: return  70; /* 033 */
        case 21: return -40; /* 111 */
        case 22: return -50; /* 112 */
        case 26: return  30; /* 122 */
        case 31: return  30; /* 133 */
        case 42: return -15; /* 222 */
        case 47: return -50; /* 233 */
        case 63: return -65; /* 333 */
        case 27: return  -1; /* 123 */
        default: return 0;
    }
}

/* The sparse integer quartic tensor B from the manuscript. */
static int64_t B4(int a, int b, int c, int d) {
    sort4(&a,&b,&c,&d);
    const int key = 64*a + 16*b + 4*c + d;
    switch (key) {
        case   7: return -2; /* 0013 */
        case  15: return -4; /* 0033 */
        case  31: return  3; /* 0133 */
        case  43: return -1; /* 0223 */
        case  91: return -1; /* 1123 */
        case  95: return -1; /* 1133 */
        case 107: return -1; /* 1223 */
        case 111: return  1; /* 1233 */
        case 191: return -6; /* 2333 */
        default: return 0;
    }
}

static void add_product3(int64_t *out, int64_t tc,
                         const SPoly *a, const SPoly *b, const SPoly *c) {
    for (int i=0;i<a->n;i++) for (int j=0;j<b->n;j++) for (int k=0;k<c->n;k++) {
        const int code = a->code[i] + b->code[j] + c->code[k];
        out[code] += tc * (int64_t)a->coeff[i] * b->coeff[j] * c->coeff[k];
    }
}

static void add_product4(int64_t *out, int64_t tc,
                         const SPoly *a, const SPoly *b,
                         const SPoly *c, const SPoly *d) {
    for (int i=0;i<a->n;i++) for (int j=0;j<b->n;j++)
    for (int k=0;k<c->n;k++) for (int l=0;l<d->n;l++) {
        const int code = a->code[i] + b->code[j] + c->code[k] + d->code[l];
        out[code] += tc * (int64_t)a->coeff[i] * b->coeff[j]
                            * c->coeff[k] * d->coeff[l];
    }
}

static int64_t *build_power_polynomials(void) {
    SPoly q[4][4];
    build_frame(q);
    int64_t *power = (int64_t *)calloc((size_t)NCOMP*NCTRL, sizeof(int64_t));
    if (!power) die("out of memory");

    for (int omit=0; omit<4; omit++) {
        int cols[3], m=0;
        for (int j=0;j<4;j++) if (j!=omit) cols[m++]=j;
        int64_t *out = power + (size_t)omit*NCTRL;
        for (int a=0;a<4;a++) for (int b=0;b<4;b++) for (int c=0;c<4;c++) {
            int64_t tc=A5(a,b,c);
            if (tc) add_product3(out,tc,&q[a][cols[0]],&q[b][cols[1]],&q[c][cols[2]]);
        }
    }
    int64_t *out = power + (size_t)4*NCTRL;
    for (int a=0;a<4;a++) for (int b=0;b<4;b++)
    for (int c=0;c<4;c++) for (int d=0;d<4;d++) {
        int64_t tc=B4(a,b,c,d);
        if (tc) add_product4(out,tc,&q[a][0],&q[b][1],&q[c][2],&q[d][3]);
    }
    return power;
}

static int count_nonzero(const int64_t *a) {
    int n=0;
    for (int i=0;i<NCTRL;i++) if (a[i]) n++;
    return n;
}

static int64_t checked_i128_to_i64(__int128 x) {
    if (x > INT64_MAX || x < INT64_MIN) die("integer overflow");
    return (int64_t)x;
}

static void transform_axis(const int64_t *src, int64_t *dst, int axis) {
    const int stride=STRIDE[axis];
    const int block=stride*SIDE;
    for (int base=0;base<NCTRL;base+=block) {
        for (int off=0;off<stride;off++) {
            int64_t v[SIDE];
            for (int e=0;e<SIDE;e++) v[e]=src[base+e*stride+off];
            for (int i=0;i<SIDE;i++) {
                __int128 s=0;
                for (int e=0;e<SIDE;e++) s+=(__int128)T6[i][e]*v[e];
                dst[base+i*stride+off]=checked_i128_to_i64(s);
            }
        }
    }
}

static int64_t floor_ratio(__int128 n, int64_t d) {
    __int128 q=n/d, r=n%d;
    if (r && n<0) q--;
    return checked_i128_to_i64(q);
}

static int64_t ceil_ratio(__int128 n, int64_t d) {
    __int128 q=n/d, r=n%d;
    if (r && n>0) q++;
    return checked_i128_to_i64(q);
}

static Net build_root_net(void) {
    int64_t *power=build_power_polynomials();
    const int expected[NCOMP]={316,316,316,316,832};
    for (int c=0;c<NCOMP;c++) {
        int n=count_nonzero(power+(size_t)c*NCTRL);
        if (n!=expected[c]) {
            fprintf(stderr,"component %d has %d power terms, expected %d\n",c,n,expected[c]);
            die("tensor-to-polynomial expansion mismatch");
        }
    }

    int64_t *tmp1=(int64_t *)xmalloc((size_t)NCTRL*sizeof(int64_t));
    int64_t *tmp2=(int64_t *)xmalloc((size_t)NCTRL*sizeof(int64_t));
    Net net;
    net.lo=(int64_t *)xmalloc((size_t)NCOMP*NCTRL*sizeof(int64_t));
    net.hi=(int64_t *)xmalloc((size_t)NCOMP*NCTRL*sizeof(int64_t));

    for (int c=0;c<NCOMP;c++) {
        memcpy(tmp1,power+(size_t)c*NCTRL,(size_t)NCTRL*sizeof(int64_t));
        int64_t *src=tmp1,*dst=tmp2,*swap;
        for (int axis=0;axis<NV;axis++) {
            transform_axis(src,dst,axis);
            swap=src;src=dst;dst=swap;
        }
        int64_t den=0;
        for (int i=0;i<NCTRL;i++) {
            int64_t a=src[i] < 0 ? -src[i] : src[i];
            if (a>den) den=a;
        }
        if (!den) die("zero polynomial component");
        for (int i=0;i<NCTRL;i++) {
            __int128 n=(__int128)src[i]*FIXED_SCALE;
            net.lo[(size_t)c*NCTRL+i]=floor_ratio(n,den);
            net.hi[(size_t)c*NCTRL+i]=ceil_ratio(n,den);
        }
    }
    free(power); free(tmp1); free(tmp2);
    return net;
}

static void free_net(Net *n) {
    free(n->lo); free(n->hi); n->lo=n->hi=NULL;
}

static int64_t floor16(__int128 n) { return floor_ratio(n,16); }
static int64_t ceil16(__int128 n) { return ceil_ratio(n,16); }

static void split_line(const int64_t b[5], int64_t l[5], int64_t r[5], int upper) {
#define ROUND(X) (upper ? ceil16((X)) : floor16((X)))
    l[0]=b[0];
    l[1]=ROUND((__int128)8*b[0]+8*b[1]);
    l[2]=ROUND((__int128)4*b[0]+8*b[1]+4*b[2]);
    l[3]=ROUND((__int128)2*b[0]+6*b[1]+6*b[2]+2*b[3]);
    l[4]=ROUND((__int128)b[0]+4*b[1]+6*b[2]+4*b[3]+b[4]);
    r[0]=l[4];
    r[1]=ROUND((__int128)2*b[1]+6*b[2]+6*b[3]+2*b[4]);
    r[2]=ROUND((__int128)4*b[2]+8*b[3]+4*b[4]);
    r[3]=ROUND((__int128)8*b[3]+8*b[4]);
    r[4]=b[4];
#undef ROUND
}

static void split_array(const int64_t *src, int64_t *left, int64_t *right,
                        int axis, int upper) {
    const int stride=STRIDE[axis];
    const int block=stride*SIDE;
    for (int base=0;base<NCTRL;base+=block) {
        for (int off=0;off<stride;off++) {
            int64_t b[5],l[5],r[5];
            for (int e=0;e<5;e++) b[e]=src[base+e*stride+off];
            split_line(b,l,r,upper);
            for (int e=0;e<5;e++) {
                left[base+e*stride+off]=l[e];
                right[base+e*stride+off]=r[e];
            }
        }
    }
}

static void split_net(Net *parent, int axis, Net *left, Net *right) {
    const size_t bytes=(size_t)NCOMP*NCTRL*sizeof(int64_t);
    left->lo=(int64_t *)xmalloc(bytes); left->hi=(int64_t *)xmalloc(bytes);
    right->lo=(int64_t *)xmalloc(bytes); right->hi=(int64_t *)xmalloc(bytes);
    for (int c=0;c<NCOMP;c++) {
        split_array(parent->lo+(size_t)c*NCTRL,
                    left->lo+(size_t)c*NCTRL,right->lo+(size_t)c*NCTRL,axis,0);
        split_array(parent->hi+(size_t)c*NCTRL,
                    left->hi+(size_t)c*NCTRL,right->hi+(size_t)c*NCTRL,axis,1);
    }
}

static int axis_from_name(const char *name) {
    static const char *axis_names[NV] = {"x1","x2","x3","y1","y2","y3"};
    for (int i=0;i<NV;i++) if (strcmp(name,axis_names[i])==0) return i;
    return -1;
}

static Certificate load_certificate(const char *path) {
    FILE *f=fopen(path,"r");
    if (!f) { perror(path); exit(1); }

    char line[512];
    if (!fgets(line,sizeof(line),f)) die("empty text certificate");
    uint64_t nodes=0, leaves=0, scale=0;
    if (sscanf(line,
               "G4BCERT-TEXT-1 nodes=%" SCNu64 " leaves=%" SCNu64
               " witness_scale=%" SCNu64,
               &nodes,&leaves,&scale)!=3)
        die("bad text certificate header");
    if (nodes != 2*leaves-1) die("certificate tree count is inconsistent");
    if (scale != (UINT64_C(1)<<52)) die("unexpected witness scale");

    Certificate c;
    memset(&c,0,sizeof(c));
    c.ntokens=nodes;
    c.nleaves=leaves;
    c.tokens=(uint8_t *)xmalloc((size_t)nodes);
    c.witnesses=(int64_t *)xmalloc((size_t)leaves*NCOMP*sizeof(int64_t));

    uint64_t node_pos=0, leaf_pos=0;
    while (fgets(line,sizeof(line),f)) {
        char kind[16];
        if (sscanf(line,"%15s",kind)!=1) continue;
        if (node_pos>=nodes) die("too many node records in text certificate");
        if (strcmp(kind,"split")==0) {
            char axis_name[16], extra[16];
            int n=sscanf(line,"split %15s %15s",axis_name,extra);
            if (n!=1) die("malformed split record in text certificate");
            int axis=axis_from_name(axis_name);
            if (axis<0) die("unknown split axis in text certificate");
            c.tokens[node_pos++]=(uint8_t)axis;
        } else if (strcmp(kind,"leaf")==0) {
            if (leaf_pos>=leaves) die("too many leaf records in text certificate");
            int64_t *w=c.witnesses+(size_t)leaf_pos*NCOMP;
            char extra[16];
            int n=sscanf(line,
                         "leaf %" SCNd64 " %" SCNd64 " %" SCNd64
                         " %" SCNd64 " %" SCNd64 " %15s",
                         &w[0],&w[1],&w[2],&w[3],&w[4],extra);
            if (n!=5) die("malformed leaf record in text certificate");
            int nonzero=0;
            for (int i=0;i<NCOMP;i++) if (w[i]) nonzero=1;
            if (!nonzero) die("zero leaf witness in text certificate");
            c.tokens[node_pos++]=LEAF_TOKEN;
            leaf_pos++;
        } else {
            die("unknown record in text certificate");
        }
    }
    fclose(f);

    if (node_pos!=nodes) die("text certificate node count does not match header");
    if (leaf_pos!=leaves) die("text certificate leaf count does not match header");
    c.min_leaf_sum=((__int128)1)<<126;
    return c;
}

static int verify_leaf(const Net *net, const int64_t w[5], __int128 *minimum) {
    int nonzero=0;
    for (int c=0;c<NCOMP;c++) if (w[c]) nonzero=1;
    if (!nonzero) return 0;
    __int128 m=((__int128)1)<<126;
    for (int i=0;i<NCTRL;i++) {
        __int128 s=0;
        for (int c=0;c<NCOMP;c++) {
            const int64_t a = (w[c]>=0)
                ? net->lo[(size_t)c*NCTRL+i]
                : net->hi[(size_t)c*NCTRL+i];
            s += (__int128)w[c]*a;
        }
        if (s<=0) return 0;
        if (s<m) m=s;
    }
    *minimum=m;
    return 1;
}

static void verify_node(Certificate *c, Net net, int depth) {
    c->checked_nodes++;
    if (depth>c->max_depth) c->max_depth=depth;
    if (c->token_pos>=c->ntokens) die("tree ended prematurely");
    uint8_t token=c->tokens[c->token_pos++];
    if (token==LEAF_TOKEN) {
        if (c->leaf_pos>=c->nleaves) die("too many leaves in tree");
        const int64_t *w=c->witnesses+(size_t)c->leaf_pos*NCOMP;
        __int128 m;
        if (!verify_leaf(&net,w,&m)) {
            fprintf(stderr,"leaf %" PRIu64 " at depth %d failed\n",c->leaf_pos,depth);
            die("invalid separating witness");
        }
        if (m<c->min_leaf_sum) c->min_leaf_sum=m;
        c->leaf_pos++; c->checked_leaves++;
        if (c->checked_leaves%5000==0) {
            printf("checked %" PRIu64 " / %" PRIu64 " leaves\n",
                   c->checked_leaves,c->nleaves);
            fflush(stdout);
        }
        free_net(&net);
        return;
    }
    if (token>=NV) die("invalid split-axis token");
    Net left,right;
    split_net(&net,token,&left,&right);
    free_net(&net);
    verify_node(c,left,depth+1);
    verify_node(c,right,depth+1);
}

static void print_i128(__int128 x) {
    if (x==0) { putchar('0'); return; }
    if (x<0) { putchar('-'); x=-x; }
    char s[64]; int n=0;
    while (x) { s[n++]=(char)('0'+x%10); x/=10; }
    while (n) putchar(s[--n]);
}

int main(int argc, char **argv) {
    if (argc!=2) {
        fprintf(stderr,"usage: %s grunbaum_certificate.txt\n",argv[0]);
        return 2;
    }
    printf("Deriving the five chart polynomials from A and B...\n");
    Net root=build_root_net();
    printf("Loading and checking the human-readable certificate...\n");
    Certificate c=load_certificate(argv[1]);
    verify_node(&c,root,0);
    if (c.token_pos!=c.ntokens || c.leaf_pos!=c.nleaves)
        die("certificate was not consumed exactly");
    if (c.checked_nodes!=c.ntokens || c.checked_leaves!=c.nleaves)
        die("internal count mismatch");
    printf("CERTIFIED\n");
    printf("  nodes:      %" PRIu64 "\n",c.checked_nodes);
    printf("  leaves:     %" PRIu64 "\n",c.checked_leaves);
    printf("  max depth:  %d\n",c.max_depth);
    printf("  minimum exact integer leaf margin: ");
    print_i128(c.min_leaf_sum); putchar('\n');
    printf("The five tensor contractions have no common zero on [-1,1]^6.\n");
    free(c.tokens); free(c.witnesses);
    return 0;
}
