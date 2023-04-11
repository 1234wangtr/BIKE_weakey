#pragma once
#include <omp.h>

const int n = 12323* 2;
const int r = n / 2;
const int w = 142;// (int)sqrt(n);
const int v = w / 2;
const int t = 134;

/*
const int n = 15 * 2;
const int r = n / 2;
const int w = 4;// (int)sqrt(n);
const int v = w / 2;
const int t = 2;
*/
const int dbg_level = 0;
const int max_iter = 5;

extern long long errFlipNum[r][max_iter+2];
extern long long corFlipNum[r][max_iter + 2];
extern long long errRemNum[r][max_iter + 2];
extern long long corRemNum[r][max_iter + 2];

extern char poly1[r];
#pragma omp threadprivate(poly1)
extern char poly2[r];
#pragma omp threadprivate(poly2)

extern char init_err[n];
//#pragma omp threadprivate(init_err)

extern char init_synd[r];
#pragma omp threadprivate(init_synd)
//#pragma omp threadprivate(poly1,poly2,init_synd,init_err)

extern char poly[n];

extern char all_zero[r];

extern int err_weight[n];

void vectorAdd(char a[], char b[], char c[], int length);

void getRandVec(char a[],int len, int weight);

void printVec(char a[], int len);

void matrixMult(char p1[], char p2[], char err[], char res[], int r);

int checkEqual(char a[], char b[], int len);

void checkErr(char candidate[]);

void getBerVec(char a[], int n, int t);

int getErrNum(char candidate[]);

int getWeiNum(char candidate[]);

void printTwoNum();
void printOneNum();
int getTwoNum(int p1[], int p2[], char err[]);
int getOneNum(int p1[], int p2[], char err[]);

extern int two_num_succ[n];
extern int two_num_fail[n];

extern int one_num_succ[n];
extern int one_num_fail[n];
