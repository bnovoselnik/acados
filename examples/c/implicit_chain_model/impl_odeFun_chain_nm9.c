/* This function was automatically generated by CasADi */
#ifdef __cplusplus
extern "C" {
#endif

#ifdef CODEGEN_PREFIX
  #define NAMESPACE_CONCAT(NS, ID) _NAMESPACE_CONCAT(NS, ID)
  #define _NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else /* CODEGEN_PREFIX */
  #define CASADI_PREFIX(ID) impl_odeFun_chain_nm9_ ## ID
#endif /* CODEGEN_PREFIX */

#include <math.h>

#ifndef real_t
#define real_t double
#endif /* real_t */

#define to_double(x) (double) x
#define to_int(x) (int) x
/* Pre-c99 compatibility */
#if __STDC_VERSION__ < 199901L
real_t CASADI_PREFIX(fmin)(real_t x, real_t y) { return x<y ? x : y;}
#define fmin(x,y) CASADI_PREFIX(fmin)(x,y)
real_t CASADI_PREFIX(fmax)(real_t x, real_t y) { return x>y ? x : y;}
#define fmax(x,y) CASADI_PREFIX(fmax)(x,y)
#endif

#define PRINTF printf
real_t CASADI_PREFIX(sq)(real_t x) { return x*x;}
#define sq(x) CASADI_PREFIX(sq)(x)

real_t CASADI_PREFIX(sign)(real_t x) { return x<0 ? -1 : x>0 ? 1 : x;}
#define sign(x) CASADI_PREFIX(sign)(x)

static const int CASADI_PREFIX(s0)[52] = {48, 1, 0, 48, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47};
#define s0 CASADI_PREFIX(s0)
static const int CASADI_PREFIX(s1)[7] = {3, 1, 0, 3, 0, 1, 2};
#define s1 CASADI_PREFIX(s1)
/* impl_odeFun_chain_nm9 */
int impl_odeFun_chain_nm9(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  real_t a0=arg[1] ? arg[1][0] : 0;
  real_t a1=arg[0] ? arg[0][3] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][0]=a0;
  a0=arg[1] ? arg[1][1] : 0;
  a1=arg[0] ? arg[0][4] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][1]=a0;
  a0=arg[1] ? arg[1][2] : 0;
  a1=arg[0] ? arg[0][5] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][2]=a0;
  a0=arg[0] ? arg[0][6] : 0;
  a1=arg[0] ? arg[0][0] : 0;
  real_t a2=(a0-a1);
  real_t a3=sq(a2);
  real_t a4=arg[0] ? arg[0][7] : 0;
  real_t a5=arg[0] ? arg[0][1] : 0;
  real_t a6=(a4-a5);
  real_t a7=sq(a6);
  a3=(a3+a7);
  a7=arg[0] ? arg[0][8] : 0;
  real_t a8=arg[0] ? arg[0][2] : 0;
  real_t a9=(a7-a8);
  real_t a10=sq(a9);
  a3=(a3+a10);
  a3=sqrt(a3);
  a10=3.3000000000000002e-02;
  a3=(a10/a3);
  real_t a11=1.;
  a3=(a11-a3);
  a2=(a3*a2);
  real_t a12=sq(a1);
  real_t a13=sq(a5);
  a12=(a12+a13);
  a13=sq(a8);
  a12=(a12+a13);
  a12=sqrt(a12);
  a12=(a10/a12);
  a12=(a11-a12);
  a1=(a12*a1);
  a1=(a2-a1);
  a13=3.3333333333333336e+01;
  a1=(a13*a1);
  real_t a14=arg[1] ? arg[1][3] : 0;
  a14=(a14-a1);
  if (res[0]!=0) res[0][3]=a14;
  a6=(a3*a6);
  a5=(a12*a5);
  a5=(a6-a5);
  a5=(a13*a5);
  a14=arg[1] ? arg[1][4] : 0;
  a14=(a14-a5);
  if (res[0]!=0) res[0][4]=a14;
  a3=(a3*a9);
  a12=(a12*a8);
  a12=(a3-a12);
  a12=(a13*a12);
  a8=9.8100000000000005e+00;
  a12=(a12-a8);
  a9=arg[1] ? arg[1][5] : 0;
  a9=(a9-a12);
  if (res[0]!=0) res[0][5]=a9;
  a9=arg[1] ? arg[1][6] : 0;
  a12=arg[0] ? arg[0][9] : 0;
  a9=(a9-a12);
  if (res[0]!=0) res[0][6]=a9;
  a9=arg[1] ? arg[1][7] : 0;
  a12=arg[0] ? arg[0][10] : 0;
  a9=(a9-a12);
  if (res[0]!=0) res[0][7]=a9;
  a9=arg[1] ? arg[1][8] : 0;
  a12=arg[0] ? arg[0][11] : 0;
  a9=(a9-a12);
  if (res[0]!=0) res[0][8]=a9;
  a9=arg[0] ? arg[0][12] : 0;
  a0=(a9-a0);
  a12=sq(a0);
  a14=arg[0] ? arg[0][13] : 0;
  a4=(a14-a4);
  a5=sq(a4);
  a12=(a12+a5);
  a5=arg[0] ? arg[0][14] : 0;
  a7=(a5-a7);
  a1=sq(a7);
  a12=(a12+a1);
  a12=sqrt(a12);
  a12=(a10/a12);
  a12=(a11-a12);
  a0=(a12*a0);
  a2=(a0-a2);
  a2=(a13*a2);
  a1=arg[1] ? arg[1][9] : 0;
  a1=(a1-a2);
  if (res[0]!=0) res[0][9]=a1;
  a4=(a12*a4);
  a6=(a4-a6);
  a6=(a13*a6);
  a1=arg[1] ? arg[1][10] : 0;
  a1=(a1-a6);
  if (res[0]!=0) res[0][10]=a1;
  a12=(a12*a7);
  a3=(a12-a3);
  a3=(a13*a3);
  a3=(a3-a8);
  a7=arg[1] ? arg[1][11] : 0;
  a7=(a7-a3);
  if (res[0]!=0) res[0][11]=a7;
  a7=arg[1] ? arg[1][12] : 0;
  a3=arg[0] ? arg[0][15] : 0;
  a7=(a7-a3);
  if (res[0]!=0) res[0][12]=a7;
  a7=arg[1] ? arg[1][13] : 0;
  a3=arg[0] ? arg[0][16] : 0;
  a7=(a7-a3);
  if (res[0]!=0) res[0][13]=a7;
  a7=arg[1] ? arg[1][14] : 0;
  a3=arg[0] ? arg[0][17] : 0;
  a7=(a7-a3);
  if (res[0]!=0) res[0][14]=a7;
  a7=arg[0] ? arg[0][18] : 0;
  a9=(a7-a9);
  a3=sq(a9);
  a1=arg[0] ? arg[0][19] : 0;
  a14=(a1-a14);
  a6=sq(a14);
  a3=(a3+a6);
  a6=arg[0] ? arg[0][20] : 0;
  a5=(a6-a5);
  a2=sq(a5);
  a3=(a3+a2);
  a3=sqrt(a3);
  a3=(a10/a3);
  a3=(a11-a3);
  a9=(a3*a9);
  a0=(a9-a0);
  a0=(a13*a0);
  a2=arg[1] ? arg[1][15] : 0;
  a2=(a2-a0);
  if (res[0]!=0) res[0][15]=a2;
  a14=(a3*a14);
  a4=(a14-a4);
  a4=(a13*a4);
  a2=arg[1] ? arg[1][16] : 0;
  a2=(a2-a4);
  if (res[0]!=0) res[0][16]=a2;
  a3=(a3*a5);
  a12=(a3-a12);
  a12=(a13*a12);
  a12=(a12-a8);
  a5=arg[1] ? arg[1][17] : 0;
  a5=(a5-a12);
  if (res[0]!=0) res[0][17]=a5;
  a5=arg[1] ? arg[1][18] : 0;
  a12=arg[0] ? arg[0][21] : 0;
  a5=(a5-a12);
  if (res[0]!=0) res[0][18]=a5;
  a5=arg[1] ? arg[1][19] : 0;
  a12=arg[0] ? arg[0][22] : 0;
  a5=(a5-a12);
  if (res[0]!=0) res[0][19]=a5;
  a5=arg[1] ? arg[1][20] : 0;
  a12=arg[0] ? arg[0][23] : 0;
  a5=(a5-a12);
  if (res[0]!=0) res[0][20]=a5;
  a5=arg[0] ? arg[0][24] : 0;
  a7=(a5-a7);
  a12=sq(a7);
  a2=arg[0] ? arg[0][25] : 0;
  a1=(a2-a1);
  a4=sq(a1);
  a12=(a12+a4);
  a4=arg[0] ? arg[0][26] : 0;
  a6=(a4-a6);
  a0=sq(a6);
  a12=(a12+a0);
  a12=sqrt(a12);
  a12=(a10/a12);
  a12=(a11-a12);
  a7=(a12*a7);
  a9=(a7-a9);
  a9=(a13*a9);
  a0=arg[1] ? arg[1][21] : 0;
  a0=(a0-a9);
  if (res[0]!=0) res[0][21]=a0;
  a1=(a12*a1);
  a14=(a1-a14);
  a14=(a13*a14);
  a0=arg[1] ? arg[1][22] : 0;
  a0=(a0-a14);
  if (res[0]!=0) res[0][22]=a0;
  a12=(a12*a6);
  a3=(a12-a3);
  a3=(a13*a3);
  a3=(a3-a8);
  a6=arg[1] ? arg[1][23] : 0;
  a6=(a6-a3);
  if (res[0]!=0) res[0][23]=a6;
  a6=arg[1] ? arg[1][24] : 0;
  a3=arg[0] ? arg[0][27] : 0;
  a6=(a6-a3);
  if (res[0]!=0) res[0][24]=a6;
  a6=arg[1] ? arg[1][25] : 0;
  a3=arg[0] ? arg[0][28] : 0;
  a6=(a6-a3);
  if (res[0]!=0) res[0][25]=a6;
  a6=arg[1] ? arg[1][26] : 0;
  a3=arg[0] ? arg[0][29] : 0;
  a6=(a6-a3);
  if (res[0]!=0) res[0][26]=a6;
  a6=arg[0] ? arg[0][30] : 0;
  a5=(a6-a5);
  a3=sq(a5);
  a0=arg[0] ? arg[0][31] : 0;
  a2=(a0-a2);
  a14=sq(a2);
  a3=(a3+a14);
  a14=arg[0] ? arg[0][32] : 0;
  a4=(a14-a4);
  a9=sq(a4);
  a3=(a3+a9);
  a3=sqrt(a3);
  a3=(a10/a3);
  a3=(a11-a3);
  a5=(a3*a5);
  a7=(a5-a7);
  a7=(a13*a7);
  a9=arg[1] ? arg[1][27] : 0;
  a9=(a9-a7);
  if (res[0]!=0) res[0][27]=a9;
  a2=(a3*a2);
  a1=(a2-a1);
  a1=(a13*a1);
  a9=arg[1] ? arg[1][28] : 0;
  a9=(a9-a1);
  if (res[0]!=0) res[0][28]=a9;
  a3=(a3*a4);
  a12=(a3-a12);
  a12=(a13*a12);
  a12=(a12-a8);
  a4=arg[1] ? arg[1][29] : 0;
  a4=(a4-a12);
  if (res[0]!=0) res[0][29]=a4;
  a4=arg[1] ? arg[1][30] : 0;
  a12=arg[0] ? arg[0][33] : 0;
  a4=(a4-a12);
  if (res[0]!=0) res[0][30]=a4;
  a4=arg[1] ? arg[1][31] : 0;
  a12=arg[0] ? arg[0][34] : 0;
  a4=(a4-a12);
  if (res[0]!=0) res[0][31]=a4;
  a4=arg[1] ? arg[1][32] : 0;
  a12=arg[0] ? arg[0][35] : 0;
  a4=(a4-a12);
  if (res[0]!=0) res[0][32]=a4;
  a4=arg[0] ? arg[0][36] : 0;
  a6=(a4-a6);
  a12=sq(a6);
  a9=arg[0] ? arg[0][37] : 0;
  a0=(a9-a0);
  a1=sq(a0);
  a12=(a12+a1);
  a1=arg[0] ? arg[0][38] : 0;
  a14=(a1-a14);
  a7=sq(a14);
  a12=(a12+a7);
  a12=sqrt(a12);
  a12=(a10/a12);
  a12=(a11-a12);
  a6=(a12*a6);
  a5=(a6-a5);
  a5=(a13*a5);
  a7=arg[1] ? arg[1][33] : 0;
  a7=(a7-a5);
  if (res[0]!=0) res[0][33]=a7;
  a0=(a12*a0);
  a2=(a0-a2);
  a2=(a13*a2);
  a7=arg[1] ? arg[1][34] : 0;
  a7=(a7-a2);
  if (res[0]!=0) res[0][34]=a7;
  a12=(a12*a14);
  a3=(a12-a3);
  a3=(a13*a3);
  a3=(a3-a8);
  a14=arg[1] ? arg[1][35] : 0;
  a14=(a14-a3);
  if (res[0]!=0) res[0][35]=a14;
  a14=arg[1] ? arg[1][36] : 0;
  a3=arg[0] ? arg[0][39] : 0;
  a14=(a14-a3);
  if (res[0]!=0) res[0][36]=a14;
  a14=arg[1] ? arg[1][37] : 0;
  a3=arg[0] ? arg[0][40] : 0;
  a14=(a14-a3);
  if (res[0]!=0) res[0][37]=a14;
  a14=arg[1] ? arg[1][38] : 0;
  a3=arg[0] ? arg[0][41] : 0;
  a14=(a14-a3);
  if (res[0]!=0) res[0][38]=a14;
  a14=arg[0] ? arg[0][42] : 0;
  a14=(a14-a4);
  a4=sq(a14);
  a3=arg[0] ? arg[0][43] : 0;
  a3=(a3-a9);
  a9=sq(a3);
  a4=(a4+a9);
  a9=arg[0] ? arg[0][44] : 0;
  a9=(a9-a1);
  a1=sq(a9);
  a4=(a4+a1);
  a4=sqrt(a4);
  a10=(a10/a4);
  a11=(a11-a10);
  a14=(a11*a14);
  a14=(a14-a6);
  a14=(a13*a14);
  a6=arg[1] ? arg[1][39] : 0;
  a6=(a6-a14);
  if (res[0]!=0) res[0][39]=a6;
  a3=(a11*a3);
  a3=(a3-a0);
  a3=(a13*a3);
  a0=arg[1] ? arg[1][40] : 0;
  a0=(a0-a3);
  if (res[0]!=0) res[0][40]=a0;
  a11=(a11*a9);
  a11=(a11-a12);
  a13=(a13*a11);
  a13=(a13-a8);
  a8=arg[1] ? arg[1][41] : 0;
  a8=(a8-a13);
  if (res[0]!=0) res[0][41]=a8;
  a8=arg[1] ? arg[1][42] : 0;
  a13=arg[0] ? arg[0][45] : 0;
  a8=(a8-a13);
  if (res[0]!=0) res[0][42]=a8;
  a8=arg[1] ? arg[1][43] : 0;
  a13=arg[0] ? arg[0][46] : 0;
  a8=(a8-a13);
  if (res[0]!=0) res[0][43]=a8;
  a8=arg[1] ? arg[1][44] : 0;
  a13=arg[0] ? arg[0][47] : 0;
  a8=(a8-a13);
  if (res[0]!=0) res[0][44]=a8;
  a8=arg[1] ? arg[1][45] : 0;
  a13=arg[2] ? arg[2][0] : 0;
  a8=(a8-a13);
  if (res[0]!=0) res[0][45]=a8;
  a8=arg[1] ? arg[1][46] : 0;
  a13=arg[2] ? arg[2][1] : 0;
  a8=(a8-a13);
  if (res[0]!=0) res[0][46]=a8;
  a8=arg[1] ? arg[1][47] : 0;
  a13=arg[2] ? arg[2][2] : 0;
  a8=(a8-a13);
  if (res[0]!=0) res[0][47]=a8;
  return 0;
}

void impl_odeFun_chain_nm9_incref(void) {
}

void impl_odeFun_chain_nm9_decref(void) {
}

int impl_odeFun_chain_nm9_n_in(void) { return 3;}

int impl_odeFun_chain_nm9_n_out(void) { return 1;}

const char* impl_odeFun_chain_nm9_name_in(int i){
  switch (i) {
  case 0: return "i0";
  case 1: return "i1";
  case 2: return "i2";
  default: return 0;
  }
}

const char* impl_odeFun_chain_nm9_name_out(int i){
  switch (i) {
  case 0: return "o0";
  default: return 0;
  }
}

const int* impl_odeFun_chain_nm9_sparsity_in(int i) {
  switch (i) {
  case 0: return s0;
  case 1: return s0;
  case 2: return s1;
  default: return 0;
  }
}

const int* impl_odeFun_chain_nm9_sparsity_out(int i) {
  switch (i) {
  case 0: return s0;
  default: return 0;
  }
}

int impl_odeFun_chain_nm9_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 3;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 15;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
