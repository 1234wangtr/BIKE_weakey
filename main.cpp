#include <iostream>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "intersections.h"


#include "upcStatistical.h"
#include "util.h"
#include "fastbbf.h"
#include <omp.h>

#include <algorithm>

using namespace std;

int small_r = 4100;
int out_of = 1;
int only_key = 0;	// set only_key=1 iff key is special and error is rand
int tot_rand = 0;	// set tot_rand=1 iff key is rand and error is special (the DFR is almost average DFR)
int rand_err = 0;	// set rand_err=1 iff err is rand



char poly1[r];
char poly2[r];
char poly[n];
char init_err[n];
char init_synd[r];
//#pragma omp threadprivate(poly1,poly2,init_synd,init_err)

char all_zero[r];
const int tau = 3; // T T-tau for threshold

long long errFlipNum[r][max_iter + 2];
long long corFlipNum[r][max_iter + 2];
long long errRemNum[r][max_iter + 2];
long long corRemNum[r][max_iter + 2];

int err_weight[n];

int two_num_succ[n];		// 记录synd=2的数量 
int two_num_fail[n];

int one_num_succ[n];		// 记录synd=1的数量 
int one_num_fail[n];

/*
what is matrix be like?

p1[0]	p1[r-1]	...	...	p1[1]	p2[0]	p2[r-1]	...	...	p2[1]
p1[1]	p1[0]					p2[1]
p1[2]	p1[1]					p2[2]
...
...
...
p1[r-1]	p1[r-2]	...	...	p1[0]	p2[r-1]

*/

int bitFlippingParallel();
int bitFlippingSteps();

void BFMaskedIter(char tmp_s[], char y[], char mask[], int T);
void BFIter(char tmp_s[], char y[], char black[], char gray[], int T);
int blackGrayFlip();
int blackGrayFlip2();
int sortedFlip();

void printParameters() {
	printf("n=%d r=%d w=%d t=%d\n", n, r, w, t);
}

void initAll() {
	printParameters();
	
	//for (int i = 0; i < 10; i++)poly1[i] = 1;
	//getRandVec(poly1+10, r-10, v-10);
	getRandVec(poly1, r, v);
	getRandVec(poly2, r, v);
	getRandVec(init_err, n, t);
	getRandVec(init_synd, r, 0);

	getRandVec(all_zero, r, 0);

	matrixMult(poly1, poly2, init_err, init_synd, r);

	//dbg
	//printVec(poly1, r);
	//printVec(poly2, r);
}

int test_bfp() {
	//printParameters();

	getRandVec(poly1, r, v);
	getRandVec(poly2, r, v);
	getRandVec(init_err, n, t);
	getRandVec(init_synd, r, 0);

	//printVec(poly1, r);
	//printVec(poly2, r);
	//printVec(init_err, n);

	matrixMult(poly1, poly2, init_err, init_synd, r);
	//printVec(init_synd, r);

	return bitFlippingParallel()*10+bitFlippingSteps();
}

int test_bgf() {
	getRandVec(poly1, r, v);
	getRandVec(poly2, r, v);
	getRandVec(init_err, n, t);
	getRandVec(init_synd, r, 0);

	//printVec(poly1, r);
	//printVec(poly2, r);
	//printVec(init_err, n);

	matrixMult(poly1, poly2, init_err, init_synd, r);
	//printVec(init_synd, r);
	//return blackGrayFlip2();
	return blackGrayFlip2Fast()*10+ blackGrayFlip2();
}

int test_bbfast() {
	getRandVec(poly1, r, v);
	getRandVec(poly2, r, v);
	//getRandVec(init_err, n, t);
	// change to beinouli init
	getBerVec(init_err,n, t);

	getRandVec(init_synd, r, 0);

	int idx1[r];
	int idx2[r];
	int i1 = 0, i2 = 0;
	for (int i = 0; i < r; i++) {
		if (poly1[i])idx1[i1++] = i;
		if (poly2[i])idx2[i2++] = i;
	}
	idx1[i1] = idx2[i2] = -1;
	//matrixMult(poly1, poly2, init_err, init_synd, r);
	matrixMultFast(idx1, idx2, init_err, init_synd, r);
	//return blackGrayFlipFast()*10+ blackGrayFlip2Fast();
	return blackGrayFlipFast();
}

int test_bbfast_weak_key() {
	
	char init_err[n];
	getRandVec(init_err, n, 0);

	

	
	getRandVec(init_synd, r, 0);

	// rand again
	getRandVec(poly1, r, 0);
	getRandVec(poly2, r, 0);
	//poly1[0] = 1;
	
	poly1[0] = 1;
	getRandVec(&poly1[1], small_r-1, v-out_of-1);
	if(small_r<r)getRandVec(&poly1[small_r], r - small_r,out_of);
	
	
	getRandVec(poly2, r, v);
	if(tot_rand){
		getRandVec(poly1,r,v);
	}

	// dbg input

	
	int idx1[r];
	int idx2[r];
	int i1 = 0, i2 = 0;
	for (int i = 0; i < r; i++) {
		if (poly1[i])idx1[i1++] = i;
		if (poly2[i])idx2[i2++] = i;
	}
	idx1[i1] = idx2[i2] = -1;
	
	
	
	if (only_key) {
		getRandVec(init_err, n, t);
	}
	else {
		init_err[0] = 1;
		getRandVec(&init_err[1], small_r - 1, t/2 - out_of - 1);
		getRandVec(&init_err[small_r], r - small_r, out_of);
		getRandVec(&init_err[r], r, t / 2);
	}
	if(tot_rand){
		init_err[0] = 1;
		getRandVec(&init_err[1], small_r - 1, t/2 - out_of - 1);
		getRandVec(&init_err[small_r], r - small_r, out_of);
		getRandVec(&init_err[r], r, t / 2);
	}
	if(rand_err){
		getRandVec(init_err,n,t);
	}

	

	int err_idx[2000];
	int err_bg = 0;
	for (int i = 0; i < r; i++) {
		if (init_err[i] == 1)err_idx[err_bg++] = i;
	}
	int err_mid = err_bg;
	for (int i = r; i < n; i++) {
		if (init_err[i] == 1)err_idx[err_mid++] = i;
	}
	



	matrixMultFast(idx1, idx2, init_err, init_synd, r);
	int res1= blackGrayFlipFast();
	int res2 = res1;
	if (res1 == 0) {
		cout << "poly1:" << endl;
		int xxx=-1000000;
		for (int i = 0; i < v; i++) {
			cout << idx1[i] << ",";
			if(xxx<= ( (r-idx1[i] + idx1[(i+1)%v])%r  )           )xxx=(r-idx1[i] + idx1[(i+1)%v])%r;
		}
		
		cout<<endl<<"width of poly1:"<<r-xxx<<endl;
		cout <<endl<< "poly2:" << endl;
		for (int i = 0; i < v; i++) {
			cout << idx2[i] << ",";
		}
		cout << endl<<"err:" << endl;
		for (int i = 0; i < t; i++) {
			cout << err_idx[i] << ",";
		}cout << endl;

		//eps=1
		int d = err_idx[err_bg - 1];
		int a = err_idx[1];
		int c = err_idx[err_bg - 2];
		int b = err_idx[err_bg - 3];
		if (d - a < small_r || (b + r - d)<small_r) {
			cout << "bad one error" << endl;
			res2 = 1;
		}
		else {
			cout << "good one error" << endl;
		}

		d = idx1[v - 1];
		a = idx1[1];
		c = idx1[v - 2];
		b = idx1[v - 3];
		if (d - a < small_r || (b + r - d)<small_r) {
			cout << "bad one key" << endl;
			res2 = 1;
		}
		else {
			cout << "good one key" << endl;
		}
	}


	// change
	return res1 * 10 + res2;
	
}

int test_bb() {
	getRandVec(poly1, r, v);
	getRandVec(poly2, r, v);
	getRandVec(init_err, n, t);
	getRandVec(init_synd, r, 0);

	int idx1[r];
	int idx2[r];
	int i1 = 0, i2 = 0;
	for (int i = 0; i < r; i++) {
		if (poly1[i])idx1[i1++] = i;
		if (poly2[i])idx2[i2++] = i;
	}
	idx1[i1] = idx2[i2] = -1;
	//matrixMult(poly1, poly2, init_err, init_synd, r);
	matrixMultFast(idx1, idx2, init_err, init_synd, r);
	//return blackGrayFlipFast()*10+ blackGrayFlip2Fast();
	return blackGrayFlip();
}

int test_sortedf() {
	getRandVec(poly1, r, v);
	getRandVec(poly2, r, v);
	getRandVec(init_err, n, t);
	getRandVec(init_synd, r, 0);

	//printVec(poly1, r);
	//printVec(poly2, r);
	//printVec(init_err, n);

	matrixMult(poly1, poly2, init_err, init_synd, r);
	//printVec(init_synd, r);

	return sortedFlip();
}

int bitFlippingParallel() {
	// given H s, output e s.t. He=s
	char y[n];
	getRandVec(y, n, 0);
	int max_iter = 5;
	char tmp_s[r];
	for (int i = 0; i < r; i++)tmp_s[i] = init_synd[i];

	int iter = 0;
	while (!checkEqual(tmp_s,all_zero,r)) {
		iter++;
		if (iter >= max_iter)break;
		
		for (int j = 0; j < r; j++) { //deal with poly1
			// i行j列 poly1   r+i-j
			int cnt = 0;
			int not_cnt = 0;
			for (int i = 0; i < r; i++) {
				char mat_ele = poly1[(r + i - j) % r];
				if (mat_ele == 1) {
					if (tmp_s[i] == 1)cnt++;
					else not_cnt++;
				}
			}
			assert(not_cnt + cnt == v);
			if (cnt > not_cnt) {
				y[j] = 1 - y[j];
			}
		}
		for (int j = r; j < n; j++) {
			int cnt = 0;
			int not_cnt = 0;
			for (int i = 0; i < r; i++) {
				char mat_ele = poly2[(n + i - j) % r];
				if (mat_ele == 1) {
					if (tmp_s[i] == 1)cnt++;
					else not_cnt++;
				}
			}
			assert(not_cnt + cnt == v);
			if (cnt > not_cnt) {
				y[j] = 1 - y[j];
			}
		}

		//refresh tmp_s
		matrixMult(poly1, poly2, y, tmp_s, r);
		vectorAdd(tmp_s, init_synd, tmp_s,r);
	}
	 
	if (checkEqual(tmp_s, all_zero, r)) {
		int sum = 0;
		for (int i = 0; i < n; i++)sum += y[i];
		//printVec(y, n);
		cout << "decode succ with wt(e)=" <<sum<< endl;
		return 1;
	}
	else {
		cout << "decode fail" << endl;
		return 0;
	}
}

int bitFlippingSteps() {
	// given H s, output e s.t. He=s
	char y[n];
	getRandVec(y, n, 0);
	int max_iter = 5;
	char tmp_s[r];
	for (int i = 0; i < r; i++)tmp_s[i] = init_synd[i];

	int iter = 0;
	while (!checkEqual(tmp_s, all_zero, r)) {
		iter++;
		if (iter >= max_iter)break;

		for (int j = 0; j < r; j++) { //deal with poly1
									  // i行j列 poly1   r+i-j
			int cnt = 0;
			int not_cnt = 0;
			for (int i = 0; i < r; i++) {
				char mat_ele = poly1[(r + i - j) % r];
				if (mat_ele == 1) {
					if (tmp_s[i] == 1)cnt++;
					else not_cnt++;
				}
			}
			assert(not_cnt + cnt == v);
			if (cnt > not_cnt) {
				y[j] = 1 - y[j];
				// refresh it
				// column j change -> tmp_s += column[j]
				for (int i = 0; i < r; i++) {
					tmp_s[i] = (tmp_s[i]+poly1[(r-j+i)%r]) % 2;
				}
			}
		}
		for (int j = r; j < n; j++) {
			int cnt = 0;
			int not_cnt = 0;
			for (int i = 0; i < r; i++) {
				char mat_ele = poly2[(n + i - j) % r];
				if (mat_ele == 1) {
					if (tmp_s[i] == 1)cnt++;
					else not_cnt++;
				}
			}
			assert(not_cnt + cnt == v);
			if (cnt > not_cnt) {
				y[j] = 1 - y[j];
				// refresh it
				// column j change -> tmp_s += column[j]
				for (int i = 0; i < r; i++) {
					tmp_s[i] = (tmp_s[i] + poly2[(n - j + i) % r]) % 2;
				}
			}
		}

	}

	if (checkEqual(tmp_s, all_zero, r)) {
		int sum = 0;
		for (int i = 0; i < n; i++)sum += y[i];
		//printVec(y, n);
		cout << "decode succ with wt(e)=" << sum << endl;
		return 1;
	}
	else {
		cout << "decode fail" << endl;
		return 0;
	}
}

// threshold(S)=max(0.0069722S+13.53,36)

int threshold(char s[]) {
	//change
	//return v / 2;
	int sum = 0;
	for (int i = 0; i < r; i++) {
		sum += s[i];
	}
	//int res=int(0.009*sum + 11.5);
	int res = int(0.0069722*sum + 13.530);
	//if (res < v/2) res = v/2;
	// change to 36 !!!
	if (res < 36)res = 36;
	
	//int res = w / 2 + 3;
	//cout << "***dbg:S=" << sum << "	res=" << res << endl;
	return res;
}



int blackGrayFlip() {
	// given H s, output e s.t. He=s
	char y[n];
	getRandVec(y, n, 0);
	int max_iter = 6;
	char tmp_s[r];
	for (int i = 0; i < r; i++)tmp_s[i] = init_synd[i];

	int iter = 0;
	while (!checkEqual(tmp_s, all_zero, r)) {
		iter++;
		if (iter >= max_iter)break;

		int T = threshold(tmp_s);
		
		char black[n];
		char gray[n];
		getRandVec(black, n, 0);
		getRandVec(gray, n, 0);

		BFIter(tmp_s, y, black, gray, T);
		matrixMult(poly1, poly2, y, tmp_s, r);
		vectorAdd(tmp_s, init_synd, tmp_s, r);
		//dbg
		checkErr(y);

		if (iter == 1) {
			BFMaskedIter(tmp_s, y, black, (v + 1) / 2 + 1);
			matrixMult(poly1, poly2, y, tmp_s, r);
			vectorAdd(tmp_s, init_synd, tmp_s, r);
			checkErr(y);

			BFMaskedIter(tmp_s, y, gray, (v + 1) / 2 + 1);
			matrixMult(poly1, poly2, y, tmp_s, r);
			vectorAdd(tmp_s, init_synd, tmp_s, r);
			checkErr(y);
		}
	}

	if (checkEqual(tmp_s, all_zero, r)) {
		int sum = 0;
		for (int i = 0; i < n; i++)sum += y[i];
		//printVec(y, n);
		cout << "decode succ with wt(e)=" << sum << endl;
		return 1;
	}
	else {
		cout << "decode fail" << endl;
		return 0;
	}
}

int blackGrayFlip2() {
	// given H s, output e s.t. He=s
	char y[n];
	getRandVec(y, n, 0);
	int max_iter = 6;
	char tmp_s[r];
	for (int i = 0; i < r; i++)tmp_s[i] = init_synd[i];

	if (dbg_level) {
		cout << "dbg p1 in normal" << endl;
		for (int i = 0; i < r; i++) {
			cout << (int)poly1[i] << " ";
		}cout << endl;
	}

	int iter = 0;
	while (!checkEqual(tmp_s, all_zero, r)) {
		iter++;
		if (iter >= max_iter)break;

		int T = threshold(tmp_s);
		char black[n];
		char gray[n];
		getRandVec(black, n, 0);
		getRandVec(gray, n, 0);

		if (dbg_level)cout << "dbg Thre=" << T << endl;
		BFIter(tmp_s, y, black, gray, T);
		matrixMult(poly1, poly2, y, tmp_s, r);
		vectorAdd(tmp_s, init_synd, tmp_s, r);
		if (iter == 1) {
			BFMaskedIter(tmp_s, y, black, (v + 1) / 2 + 1);
			matrixMult(poly1, poly2, y, tmp_s, r);
			vectorAdd(tmp_s, init_synd, tmp_s, r);
			BFMaskedIter(tmp_s, y, black, (v + 1) / 2 + 1);
			matrixMult(poly1, poly2, y, tmp_s, r);
			vectorAdd(tmp_s, init_synd, tmp_s, r);
		}
	}

	if (checkEqual(tmp_s, all_zero, r)) {
		int sum = 0;
		for (int i = 0; i < n; i++)sum += y[i];
		//printVec(y, n);
		cout << "decode succ with wt(e)=" << sum << endl;
		return 1;
	}
	else {
		cout << "decode fail" << endl;
		return 0;
	}
}

void BFIter(char tmp_s[],char y[],char black[],char gray[],int T) {
	for (int j = 0; j < r; j++) { //deal with poly1
								  // i行j列 poly1   r+i-j
		int cnt = 0;
		for (int i = 0; i < r; i++) {
			char mat_ele = poly1[(r + i - j) % r];
			if (mat_ele == 1) {
				if (tmp_s[i] == 1) { cnt++; }
			}
		}
		if (cnt >= T) {
			y[j] = 1 - y[j];
			black[j] = 1;
			
		}
		else if (cnt >= T - tau) {
			gray[j] = 1;
		}
	}

	for (int j = r; j < n; j++) {
		int cnt = 0;
		for (int i = 0; i < r; i++) {
			char mat_ele = poly2[(n + i - j) % r];
			if (mat_ele == 1) {
				if (tmp_s[i] == 1)cnt++;
			}
		}
		
		if (cnt >= T) {
			y[j] = 1 - y[j];
			black[j] = 1;
			
		}
		else if (cnt >= T - tau) {
			gray[j] = 1;
		}
	}

	
}

void BFMaskedIter(char tmp_s[], char y[], char mask[], int T) {
	for (int j = 0; j < r; j++) { //deal with poly1
								  // i行j列 poly1   r+i-j
		int cnt = 0;
		for (int i = 0; i < r; i++) {
			char mat_ele = poly1[(r + i - j) % r];
			if (mat_ele == 1) {
				if (tmp_s[i] == 1)cnt++;
			}
		}
		if (cnt >= T) {
			if(mask[j])y[j] = 1-y[j] ;
		}
	}

	for (int j = r; j < n; j++) {
		int cnt = 0;
		for (int i = 0; i < r; i++) {
			char mat_ele = poly2[(n + i - j) % r];
			if (mat_ele == 1) {
				if (tmp_s[i] == 1)cnt++;
			}
		}

		if (cnt >= T) {
			if (mask[j])y[j] = 1 - y[j];
		}
	}
}

int sortedFlip() {
	char y[n];
	getRandVec(y, n, 0);
	int max_iter = 6;
	char tmp_s[r];
	for (int i = 0; i < r; i++)tmp_s[i] = init_synd[i];

	

	int iter = 0;
	while (!checkEqual(tmp_s, all_zero, r)) {
		iter++;
		if (iter >= max_iter)break;

		int T[n];
		for (int i = 0; i < n; i++) {
			T[i] = i;
		}
		// first of all get Cnts and sort  T[i]/n=cnt[i]  T[i]%n=i
		for (int j = 0; j < r; j++) { //deal with poly1
									  // i行j列 poly1   r+i-j
			int cnt = 0;
			int not_cnt = 0;
			for (int i = 0; i < r; i++) {
				char mat_ele = poly1[(r + i - j) % r];
				if (mat_ele == 1) {
					if (tmp_s[i] == 1)cnt++;
					else not_cnt++;
				}
			}
			
			T[j] += cnt*n;
		}
		for (int j = r; j < n; j++) {
			int cnt = 0;
			int not_cnt = 0;
			for (int i = 0; i < r; i++) {
				char mat_ele = poly2[(n + i - j) % r];
				if (mat_ele == 1) {
					if (tmp_s[i] == 1)cnt++;
					else not_cnt++;
				}
			}
			T[j] += cnt*n;
		}
		sort(T, T + n);
		for (int ii = n-1; ii >= 0; ii--) {
			int idx = T[ii] % n;
			if (idx < r) {
				int cnt = 0;
				int not_cnt = 0;
				for (int i = 0; i < r; i++) {
					char mat_ele = poly1[(r + i - idx) % r];
					if (mat_ele == 1) {
						if (tmp_s[i] == 1)cnt++;
						else not_cnt++;
					}
				}
				if (cnt > not_cnt) {
					y[idx] = 1 - y[idx];
					// refresh it
					// column j change -> tmp_s += column[j]
					for (int i = 0; i < r; i++) {
						tmp_s[i] = (tmp_s[i] + poly1[(r - idx + i) % r]) % 2;
					}
				}
			}
			else {
				int cnt = 0;
				int not_cnt = 0;
				for (int i = 0; i < r; i++) {
					char mat_ele = poly2[(n + i - idx) % r];
					if (mat_ele == 1) {
						if (tmp_s[i] == 1)cnt++;
						else not_cnt++;
					}
				}
				if (cnt > not_cnt) {
					y[idx] = 1 - y[idx];
					// refresh it
					// column j change -> tmp_s += column[j]
					for (int i = 0; i < r; i++) {
						tmp_s[i] = (tmp_s[i] + poly2[(n - idx + i) % r]) % 2;
					}
				}
				
			}
		}
	}

	if (checkEqual(tmp_s, all_zero, r)) {
		int sum = 0;
		for (int i = 0; i < n; i++)sum += y[i];
		//printVec(y, n);
		cout << "decode succ with wt(e)=" << sum << endl;
		return 1;
	}
	else {
		cout << "decode fail" << endl;
		return 0;
	}

}


int testInters() {
	initAll();
	getIntersections(poly1, poly2, r,t);
	return 0;
}

int testIntersUnion() {
	initAll();
	getUnionIntersections(poly1, poly2, r, t);
	return 0;
}



void printFlipNum() {
	cout << "print Flip Num" << endl;
	for (int i = 0; i < r; i++) {
		long long sum = 0;
		for (int j = 0; j < max_iter + 2; j++) {
			sum += errFlipNum[i][j]+errRemNum[i][j]+corFlipNum[i][j]+corRemNum[i][j];
		}
		if (sum) {
			cout << "err num=" << i << " sum="<<sum<<endl;
			for (int j = 0; j < max_iter + 2; j++) {
				printf("iter=%d errFlip=%lld errRem=%lld corFlip=%lld corRem=%lld\n",j,errFlipNum[i][j],errRemNum[i][j],corFlipNum[i][j],corRemNum[i][j]);
			}
		}
	}
}

void printAvgThres() {
	double avg[r];
	for (int i = 0; i < r; i++)avg[i] = 0;
	for (int tmpT = 0; tmpT <= t; tmpT++) {
		double sum = 0.0;
		const int max_iter = 100;
		for (int i = 0; i < max_iter; i++) {
			getRandVec(poly1, r, v);
			getRandVec(poly2, r, v);
			getRandVec(init_err, n, tmpT);
			getRandVec(init_synd, r, 0);

			int idx1[r];
			int idx2[r];
			int i1 = 0, i2 = 0;
			for (int i = 0; i < r; i++) {
				if (poly1[i])idx1[i1++] = i;
				if (poly2[i])idx2[i2++] = i;
			}
			idx1[i1] = idx2[i2] = -1;
			matrixMultFast(idx1, idx2, init_err, init_synd, r);
			sum += threshold(init_synd);
		}
		avg[tmpT] = sum / max_iter;
	}
	for (int i = 0; i <= t; i++) {
		cout << i << ":" << avg[i] << "," << endl;
	}
}

void printErrWei() {
	cout << "printErrWei:" << endl;
	for (int i = 0; i < n; i++) {
		if (err_weight[i])cout << "i=" << i << "\t sum=" << err_weight[i] << endl;
	}
}


int main() {
	srand((unsigned int)time(NULL));
	

	

	printParameters();
	cout<<"small_r="<<small_r<<endl;

	//printAvgThres();
	for (int i = 0; i < n; i++) {
		err_weight[i] = 0;
	}
	
	
	
	
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < max_iter + 2; j++) {
			errFlipNum[i][j] = errRemNum[i][j] = corFlipNum[i][j] = corRemNum[i][j] = 0;
		}
	}
	
	
	
	for (int i = 0; i < n; i++) {
		two_num_fail[i] = two_num_succ[i] = one_num_fail[i] = one_num_succ[i] = 0;
	}


	int succ1 = 0;
	int succ2 = 0;

	int tot=0;
	int numProcs = omp_get_num_procs();
	//numProcs /= 4;
	cout<<"out of="<<out_of<<endl;
	cout<<"using "<<numProcs<<" cores"<<endl;

	volatile bool flag= false;
	#pragma omp parallel for  shared(flag) num_threads(numProcs)
	for (long long i = 1; i <= 50000000000; i++) {
		if(flag)continue;
		int res = test_bbfast_weak_key();
		#pragma omp critical 
		{
		tot += 1;
		//succ1 += res;
		succ1 += res / 10;
		succ2 += res % 10; 
		if(tot-succ1 >= 16*1000)flag=true;;
		if (tot % 100000 == 0 || tot==100 || tot==1000||tot==2000||tot==5000) {
			cout<<"????--------------------------"<<endl;
			cout << "tot=" << tot << " succ1=" << succ1 << " succ2=" << succ2 << endl;
			//cout << "i=" << i << "succ=" << succ1 << endl;
		}

		}
	}
	//printFlipNum();

	//printErrWei();
	

	
	
	/*int succ = 0;
	for (int i = 0; i < 100; i++){
		succ+=testProbFlip();
		cout << "succ=" << succ << " tot=" << i + 1 << endl;
	}*/
	/*
	for (int i = 0; i < 1; i++) {
		initAll();
		upcStatistical(poly1, poly2, t, r);
	}*/
	cout << "-------------------------\ntot=" << tot << " succ1=" << succ1 << " succ2=" << succ2 << endl;
			
	system("pause");
	return 0;
}
