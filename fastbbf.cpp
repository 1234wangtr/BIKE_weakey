#include "fastbbf.h"
#include "util.h"

#include <iostream>
#include <fstream>

#include <omp.h>

using namespace std;

int tmp_iter = 0;

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

const int dbg_err_num = 0;

void matrixMultFast(int p1[], int p2[], char err[], char res[], int r) {
	// mult the matrix by p1,p2  and  err   -> store in res
	//int nThreads = 4;
        	//omp_set_num_threads( nThreads );
	//#pragma omp parallel for
	//for(int i=0;i<100;i++){
	//	//res[i]=0;
	//}
	
	//#pragma omp parallel for // shared(p1,p2,err,res) firstprivate(r) private(i)
	for (int i = 0; i < r; i++) {
		int sum = 0;
		int bg1 = 0, bg2 = 0;
		
		while (p1[bg1]!=-1) {
			if (err[(r-p1[bg1] + i) % r])sum++;
			bg1++;
		}
		while (p2[bg2] != -1) {
			if (err[r+(r-p2[bg2] + i) % r])sum++;
			bg2++;
		}
		res[i] = sum % 2;
	}
}

void BFIterFast(char tmp_s[], char y[], int black[], int gray[],int idx1[],int idx2[], int T) {
	/*
		idx1[0] = first element a0  s.t. poly1[a0]=1
		idx1[1] -> poly1[idx[1]]=1
		...
		idx1[v-1] -> poly1[idx1[v-1]]=1
		idx1[v] = -1
	*/

	//dbg
	int eNum = 0;
	for (int i = 0; i < n; i++) {
		if (y[i] != init_err[i])eNum++;
	}

	int black_iter = 0;
	int gray_iter = 0;
	const int tau = 3;

	for (int j = 0; j < r; j++) { //deal with poly1
								  // i行j列 poly1   r+i-j
		int cnt = 0;
		/*for (int i = 0; i < r; i++) {
			char mat_ele = poly1[(r + i - j) % r];
			if (mat_ele == 1) {
				if (tmp_s[i] == 1)cnt++;
			}
		}*/
		int bg1 = 0;
		while (idx1[bg1] != -1) {
			// idx1[bg1]=(i-j)%r
			if (tmp_s[(j + idx1[bg1]) % r] == 1) { cnt++; if (dbg_level)cout << "dbg upc:\t" << (j + idx1[bg1]) % r << " bad" << endl; }
			bg1++;
		}

		if (cnt >= T) {
			y[j] = 1 - y[j];
			black[black_iter++] = j;
			if (dbg_level)cout << "dbg flip " << j << endl;
		}
		else if (cnt >= T - tau) {
			gray[gray_iter++] = j;
		}

		//dbg
		if (dbg_err_num) {
			if (cnt >= T) {
				if (y[j] == init_err[j])errFlipNum[eNum][tmp_iter]++;
				else corFlipNum[eNum][tmp_iter]++;
			}
			else {
				if (y[j] == init_err[j])corRemNum[eNum][tmp_iter]++;
				else errRemNum[eNum][tmp_iter]++;
			}
		}
	}

	for (int j = r; j < n; j++) { //deal with poly2
								  // i行j列 poly2   r+i-j
		int cnt = 0;

		int bg2 = 0;
		while (idx2[bg2] != -1) {
			// idx1[bg1]=(i-j)%r
			if (tmp_s[(j + idx2[bg2]) % r] == 1)cnt++;
			bg2++;
		}

		if (cnt >= T) {
			y[j] = 1 - y[j];
			black[black_iter++] = j;
			if (dbg_level)cout << "dbg flip " << j << endl;
		}
		else if (cnt >= T - tau) {
			gray[gray_iter++] = j;
		}

		//dbg
		if (dbg_err_num) {
			if (cnt >= T) {
				if (y[j] == init_err[j])errFlipNum[eNum][tmp_iter]++;
				else corFlipNum[eNum][tmp_iter]++;
			}
			else {
				if (y[j] == init_err[j])corRemNum[eNum][tmp_iter]++;
				else errRemNum[eNum][tmp_iter]++;
			}
		}
	}
	black[black_iter] = -1;
	gray[gray_iter] = -1;
	if (dbg_level) {
		cout << "dbg black:" << endl;
		for (int i = 0; i <= black_iter; i++) {
			cout << black[i] << " ";
		}cout << endl;
	}
}

void BFMaskedIterFast(char tmp_s[], char y[], int mask[], int idx1[], int idx2[], int T) {
	//dbg
	int eNum = 0;
	for (int i = 0; i < n; i++) {
		if (y[i] != init_err[i])eNum++;
	}

	int n = 2 * r;
	int mask_iter = 0;
	while (mask[mask_iter] != -1) {
		int refresh_idx = mask[mask_iter];
		int cnt = 0;
		if (refresh_idx < r) {
			// in poly1
			int bg1 = 0;
			while (idx1[bg1] != -1) {
				// idx1[bg1]=(i-j)%r
				if (tmp_s[(refresh_idx + idx1[bg1]) % r] == 1)cnt++;
				bg1++;
			}
			if (cnt >= T) {
				y[refresh_idx] = 1 - y[refresh_idx];
			}

		}
		else {
			int bg2 = 0;
			while (idx2[bg2] != -1) {
				// idx1[bg1]=(i-j)%r
				if (tmp_s[(refresh_idx + idx2[bg2]) % r] == 1)cnt++;
				bg2++;
			}
			if (cnt >= T) {
				y[refresh_idx] = 1 - y[refresh_idx];
			}
		}
		//dbg
		if (dbg_err_num) {
			if (cnt >= T) {
				if (y[refresh_idx] == init_err[refresh_idx])errFlipNum[eNum][tmp_iter]++;
				else corFlipNum[eNum][tmp_iter]++;
			}
			else {
				if (y[refresh_idx] == init_err[refresh_idx])corRemNum[eNum][tmp_iter]++;
				else errRemNum[eNum][tmp_iter]++;
			}
		}

		mask_iter++;
	}

	
}

int thresholdFast(char s[]) {
	//return v / 2;
	int sum = 0;
	for (int i = 0; i < r; i++) {
		sum += s[i];
	}
	//int res=int(0.009*sum + 11.5);
	int res = int(0.0069722*sum + 13.530);
	if (res < v / 2+1) res = v / 2+1;
	//int res = w / 2 + 3;
	return res;
}


int blackGrayFlip2Fast() {
	// given H s, output e s.t. He=s
	
	char y[n];
	getRandVec(y, n, 0);
	int max_iter = 6;
	char tmp_s[r];
	for (int i = 0; i < r; i++)tmp_s[i] = init_synd[i];

	int idx1[r];
	int idx2[r];
	int i1 = 0, i2 = 0;
	for (int i = 0; i < r; i++) {
		if (poly1[i])idx1[i1++] = i;
		if (poly2[i])idx2[i2++] = i;
	}
	idx1[i1] = idx2[i2] = -1;
	if (dbg_level) {
		cout << "dbg:p1=" << endl;
		for (int i = 0; i < r; i++) {
			cout << (int)poly1[i] << " ";
		}cout << endl;

		cout << "dbg:i1=" << i1 << endl;
		for (int i = 0; i < i1; i++) {
			cout << idx1[i] << " ";
		}cout << endl;
	}

	int iter = 0;
	while (!checkEqual(tmp_s, all_zero, r)) {
		iter++;
		if (iter >= max_iter)break;

		int T = thresholdFast(tmp_s);
		
		int black_idx[n];
		int gray_idx[n];

		if (dbg_level)cout << "dbg Threshold:" << T << endl;
		BFIterFast(tmp_s, y, black_idx, gray_idx,idx1,idx2,T);
		//matrixMult(poly1, poly2, y, tmp_s, r);
		matrixMultFast(idx1, idx2, y, tmp_s, r);
		vectorAdd(tmp_s, init_synd, tmp_s, r);
		if (iter == 1) {
			BFMaskedIterFast(tmp_s, y, black_idx,idx1,idx2, (v + 1) / 2 + 1);
			//matrixMult(poly1, poly2, y, tmp_s, r);
			matrixMultFast(idx1, idx2, y, tmp_s, r);
			vectorAdd(tmp_s, init_synd, tmp_s, r);
			BFMaskedIterFast(tmp_s, y, black_idx, idx1, idx2, (v + 1) / 2 + 1);
			//matrixMult(poly1, poly2, y, tmp_s, r);
			matrixMultFast(idx1, idx2, y, tmp_s, r);
			vectorAdd(tmp_s, init_synd, tmp_s, r);
		}
	}

	if (checkEqual(tmp_s, all_zero, r)) {
		int sum = 0;
		for (int i = 0; i < n; i++)sum += y[i];
		//printVec(y, n);
		//cout << "decode succ with wt(e)=" << sum << endl;

		return 1;
	}
	else {
		cout << "decode fail black-black" << endl;
		return 0;
	}

	
}

int blackGrayFlipFast() {
	// given H s, output e s.t. He=s

	char y[n];
	getRandVec(y, n, 0);
	int max_iter = 6;
	char tmp_s[r];
	for (int i = 0; i < r; i++)tmp_s[i] = init_synd[i];

	int idx1[r];
	int idx2[r];
	int i1 = 0, i2 = 0;
	for (int i = 0; i < r; i++) {
		if (poly1[i])idx1[i1++] = i;
		if (poly2[i])idx2[i2++] = i;
	}
	idx1[i1] = idx2[i2] = -1;
	if (dbg_level) {
		cout << "dbg:p1=" << endl;
		for (int i = 0; i < r; i++) {
			cout << (int)poly1[i] << " ";
		}cout << endl;

		cout << "dbg:i1=" << i1 << endl;
		for (int i = 0; i < i1; i++) {
			cout << idx1[i] << " ";
		}cout << endl;
	}

	// get Two Num
	int two_num = getTwoNum(idx1,idx2,init_err);
	int one_num = getOneNum(idx1, idx2, init_err);

	int iter = 0;
	tmp_iter = 0;

	// check err num
	int err_num_path[100];
	int path_iter = 0;

	int wei_num_path[100];
	int wei_iter = 0;


	while (!checkEqual(tmp_s, all_zero, r)) {
		iter++;
		if (iter >= max_iter)break;

		int T = thresholdFast(tmp_s);
		//cout << "dbg get threshold=" << T << endl;

		int black_idx[n];
		int gray_idx[n];
		//dbg 
		//checkErr(y);
		err_num_path[path_iter++] = getErrNum(y);
		wei_num_path[wei_iter++] = getWeiNum(y);

		if (dbg_level)cout << "dbg Threshold:" << T << endl;
		BFIterFast(tmp_s, y, black_idx, gray_idx, idx1, idx2, T);


		//matrixMult(poly1, poly2, y, tmp_s, r);
		matrixMultFast(idx1, idx2, y, tmp_s, r);
		vectorAdd(tmp_s, init_synd, tmp_s, r);
		tmp_iter++;
		if (iter == 1) {
			//dbg 
			//checkErr(y);
			err_num_path[path_iter++] = getErrNum(y);
			wei_num_path[wei_iter++] = getWeiNum(y);

			BFMaskedIterFast(tmp_s, y, black_idx, idx1, idx2, (v + 1) / 2 + 1);
			//matrixMult(poly1, poly2, y, tmp_s, r);
			matrixMultFast(idx1, idx2, y, tmp_s, r);
			vectorAdd(tmp_s, init_synd, tmp_s, r);
			tmp_iter++;
			//dbg 
			//checkErr(y);
			err_num_path[path_iter++] = getErrNum(y);
			wei_num_path[wei_iter++] = getWeiNum(y);

			BFMaskedIterFast(tmp_s, y, gray_idx, idx1, idx2, (v + 1) / 2 + 1);
			//matrixMult(poly1, poly2, y, tmp_s, r);
			matrixMultFast(idx1, idx2, y, tmp_s, r);
			vectorAdd(tmp_s, init_synd, tmp_s, r);
			tmp_iter++;
		}
	}

	if (checkEqual(tmp_s, all_zero, r)) {
		int sum = 0;
		for (int i = 0; i < n; i++)sum += y[i];
		//printVec(y, n);
		//cout << "decode succ with wt(e)=" << sum << endl;
		two_num_succ[two_num]++;
		one_num_succ[one_num]++;
		return 1;
	}
	else {
		cout << "decode fail black-gray" << endl;
		two_num_fail[two_num]++;
		one_num_fail[one_num]++;
		return 0;
		// output to file
		ofstream ofs;
		ofs.open("data/output_err_4k5_00out_example.txt",ios::app);
		ofs << "------decode fail" << endl;
		int idx1[r];
		int idx2[r];
		int i1 = 0, i2 = 0;
		for (int i = 0; i < r; i++) {
			if (poly1[i])idx1[i1++] = i;
			if (poly2[i])idx2[i2++] = i;
		}
		idx1[i1] = idx2[i2] = -1;
		int erridx[n];
		int ee = 0;
		for (int i = 0; i < n; i++) {
			if (init_err[i])erridx[ee++] = i;
		}
		erridx[ee] = -1;
		ofs << "err num" << endl;
		for (int i = 0; i < 7; i++) {
			ofs << err_num_path[i]<<'\t';
		}ofs << endl;

		ofs << "wei num" << endl;
		for (int i = 0; i < 7; i++) {
			ofs <<wei_num_path[i] << '\t';
		}ofs << endl;

		ofs << "idx1:" << endl;
		int iter = 0;
		while (idx1[iter] != -1) {
			ofs << idx1[iter++] << '\t';
		}
		ofs << endl<<"idx2:" << endl;
		iter = 0;
		while (idx2[iter] != -1) {
			ofs << idx2[iter++] << '\t';
		}
		ofs << endl<<"err:" << endl;
		iter = 0;
		while (erridx[iter] != -1) {
			ofs << erridx[iter++] << '\t';
		}ofs << endl << "................" << endl;



		ofs << "poly1:" << endl;
		for (int i = 0; i < r; i++) {
			ofs << (int)(poly1[i]);
			if (i % 500 == 0)ofs << endl;
		}
		ofs << endl;
		ofs << "poly2:" << endl;
		for (int i = 0; i < r; i++) {
			ofs << (int)(poly2[i]);
			if (i % 500 == 0)ofs << endl;
		}
		ofs << endl;
		ofs << "error:" << endl;
		for (int i = 0; i < n; i++) {
			ofs << (int)(init_err[i]) ;
			if (i % 500 == 0)ofs << endl;
		}ofs << endl;

		



		return 0;
	}


}