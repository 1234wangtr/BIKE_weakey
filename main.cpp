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

/*
const int n = 12323*2;
const int r = n / 2;
const int w = 142;// (int)sqrt(n);
const int v = w/2;
const int t = 134;
*/





/*
// just for test
const int n = 6 * 2;
const int r = n / 2;
const int w = 4;// (int)sqrt(n);
const int v = w / 2;
const int t = 5;
*/

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
int small_r = 4000;
int out_of = 0;
int test_bbfast_weak_key() {
	
	char init_err[n];
	getRandVec(init_err, n, 0);

	

	
	getRandVec(init_synd, r, 0);

	// rand again
	getRandVec(poly1, r, 0);
	getRandVec(poly2, r, 0);
	//poly1[0] = 1;
	/*
	int wei_p1 = 0;
	while (wei_p1 < v) {
		int rd = rand() % (5000);
		if (poly1[(3 * rd ) % r] == 0) {
			poly1[(3 * rd ) % r] = 1;
			wei_p1++;
		}
	}*/
	getRandVec(&poly1[0], small_r, v-out_of);
	getRandVec(&poly1[small_r], r - small_r,out_of);
	//getRandVec(&poly2[0], 5000, v);
	//getRandVec(&poly1[3000], 2500, v/2+1);
	
	getRandVec(poly2, r, v);


	// dbg input

	
	int idx1[r];
	int idx2[r];
	int i1 = 0, i2 = 0;
	for (int i = 0; i < r; i++) {
		if (poly1[i])idx1[i1++] = i;
		if (poly2[i])idx2[i2++] = i;
	}
	idx1[i1] = idx2[i2] = -1;
	
	/*
	int idx1[] = { 10, 312, 351, 420, 542, 633, 1045, 1118, 1143, 1147, 1162, 1169, 1230, 1255, 1301, 1358, 1446, 1480, 1483, 1487, 1544, 1548, 1684, 1735, 1750, 1785, 1808, 1843, 1930, 1956, 1988, 2003, 2004, 2024, 2090, 2154, 2168, 2169, 2362, 2365, 2370, 2371, 2445, 2587, 2623, 2626, 2927, 3020, 3238, 3278, 3362, 3525, 3530, 3547, 3557, 3582, 3672, 3760, 3787, 3815, 3965, 3966, 4056, 4313, 4316, 4501, 4516, 4576, 4585, 4657, 4776,-1 };
	int idx2[] = { 264, 491, 508, 836, 879, 887, 1342, 1372, 1691, 2056, 2065, 2139, 2186, 2287, 2336, 2355, 2367, 2568, 2690, 2784, 2886, 2920, 3033, 3047, 3258, 3292, 3331, 3461, 3519, 3694, 3907, 4005, 4040, 4090, 4113, 4244, 4550, 4651, 4812, 4817, 4934, 4984, 5118, 5137, 5165, 5198, 5280, 5696, 5776, 6156, 6321, 6770, 6931, 7007, 7054, 7236, 7738, 8148, 8690, 8711, 8979, 9061, 9313, 9828, 10190, 10403, 11148, 11285, 11430, 11452, 11776 ,-1 };
	getRandVec(poly1, r, 0);
	getRandVec(poly2, r, 0);
	for (int i = 0; i < v; i++) {
		poly1[idx1[i]] = poly2[idx2[i]] = 1;
	}*/

	getRandVec(init_err, n, 0);
	getRandVec(&init_err[0], small_r, t/2-out_of);
	getRandVec(&init_err[small_r], r - small_r, out_of);
	//getRandVec(&init_err[3000], 2500, t / 2 / 2 + 1);
	/*
	int wei_er = 0;
	while (wei_er < t/2) {
		int rd = rand() % (5000);
		if (init_err[(3 * rd + 1)%r] == 0) {
			init_err[(3 * rd + 1) % r] = 1;
			wei_er++;
		}
	}*/


	getRandVec(&init_err[r], r, t / 2);

	// dbg input
	/*getRandVec(init_err, n, 0);
	int err_idx[] = { 4, 72, 153, 168, 211, 253, 318, 327, 382, 387, 431, 440, 499, 566, 592, 681, 849, 856, 936, 1083, 1213, 1261, 1346, 1380, 1467, 1470, 1489, 1634, 1686, 1727, 1781, 1805, 1807, 1901, 2031, 2081, 2197, 2483, 2523, 2586, 2679, 2735, 2745, 2754, 2847, 2915, 3328, 3331, 3426, 3498, 3568, 3578, 3706, 3731, 3877, 3964, 4052, 4132, 4306, 4521, 4689, 4719, 4778, 4873, 4886, 4917, 4996, 12452, 12487, 12518, 12657, 12906, 13302, 13547, 13663, 13695, 13861, 14045, 14185, 14325, 14391, 14637, 14689, 14878, 15235, 15392, 16113, 16116, 16407, 16647, 16650, 16662, 16759, 17050, 17150, 17315, 17392, 17481, 17626, 17710, 18501, 18546, 18861, 18954, 19001, 19733, 19872, 20012, 20129, 20187, 20210, 20245, 20273, 20300, 20927, 20941, 21133, 21585, 21751, 21815, 21903, 22222, 22704, 22736, 22787, 22874, 23018, 23077, 23137, 23269, 23448, 23627, 23799, 24145 };
	for (int i = 0; i < t; i++)init_err[err_idx[i]] = 1;*/

	matrixMultFast(idx1, idx2, init_err, init_synd, r);
	int res1= blackGrayFlipFast();

	//getRandVec(init_err, n, 0);
	//init_err[0] = 1;
	//getRandVec(&init_err[1], small_r - 1, t-1);
	//matrixMultFast(idx1, idx2, init_err, init_synd, r);
	//int res2 = blackGrayFlipFast();
	
	// change
	return res1 * 10 + res1;
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
	cout<<"out of="<<out_of<<endl;
	cout<<"using "<<numProcs<<" cores"<<endl;
	#pragma omp parallel for  num_threads(numProcs)
	for (int i = 1; i <= 10000000; i++) {
		int res = test_bbfast_weak_key();
		#pragma omp critical 
		{
		tot += 1;
		//succ1 += res;
		succ1 += res / 10;
		succ2 += res % 10; 
		}
		if (tot % 10000 == 0 || tot==100 || tot==1000||tot==2000||tot==5000) {
			cout << "tot=" << tot << " succ=" << succ1  << endl;
			//cout << "i=" << i << "succ=" << succ1 << endl;
		}
		if (0 && i % 1000 == 0) {
			printTwoNum();
			printOneNum();
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
	system("pause");
	return 0;
}