#include "util.h"
#include <cstdlib>
#include <iostream>

using namespace std;

int getErrNum(char candidate[]) {
	int still_error = 0;
	int false_positive = 0;
	for (int i = 0; i < n; i++) {
		if (candidate[i] && !init_err[i]) {
			false_positive++;
		}
		if (!candidate[i] && init_err[i]) {
			still_error++;
		}
	}
	return still_error + false_positive ;
}

int getWeiNum(char candidate[]) {
	int error = 0;
	
	for (int i = 0; i < n; i++) {
		if (candidate[i]) {
			error++;
		}
		
	}
	return error;
}

void checkErr(char candidate[]) {
	int still_error = 0;
	int false_positive = 0;
	for (int i = 0; i < n; i++) {
		if (candidate[i] && !init_err[i]) {
			false_positive++;
		}
		if (!candidate[i] && init_err[i]) {
			still_error++;
		}
	}
	cout << "dbg tote=" << still_error + false_positive << "\treme=" << still_error << "\tnewe=" << false_positive << endl;
}

void vectorAdd(char a[], char b[], char c[], int length) {
	for (int i = 0; i < length; i++) {
		c[i] = (a[i] + b[i]) % 2;
	}
}


void getRandVec(char a[],int len, int weight) {
	for (int i = 0; i < len; i++) {
		a[i] = 0;
	}
	int sum = 0;
	while (sum < weight) {
		int rd = rand()%len;
		if (a[rd] == 0) {
			a[rd] = 1;
			sum++;
		}
	}
	return;
}

void getBerVec(char a[], int n,int t) {
	int sum = 0;
	for (int i = 0; i < n; i++) {
		int rd = (rand()*999) % n;
		if (rd < t)a[i] = 1;
		else a[i] = 0;
		if(rd<t)sum += 1;
	}
	while (sum > t) {
		sum = 0;
		for (int i = 0; i < n; i++) {
			int rd = (rand() * 1) % n;
			if (rd < t)a[i] = 1;
			else a[i] = 0;
			if (rd<t)sum += 1;
		}
	}
	err_weight[sum]++;
	//cout << "sum=" << sum <<" t="<<t<<" n="<<n<< endl;
	return;
}

void printVec(char a[], int len) {
	cout << "dbg vec:" << endl;
	for (int i = 0; i < len; i++) {
		cout << (int)(a[i])<<" ";
	}
	cout << endl;
}

void matrixMult(char p1[], char p2[], char err[], char res[], int r) {
	// mult the matrix by p1,p2  and  err   -> store in res
	for (int i = 0; i < r; i++) {
		int sum = 0;
		for (int j = 0; j < r; j++) {
			sum += err[j] * p1[(r + i - j) % r];
		}
		for (int j = 0; j < r; j++) {
			sum += err[j + r] * p2[(r + i - j) % r];
		}
		res[i] = sum % 2;
	}
}

int checkEqual(char a[], char b[], int len) {
	for (int i = 0; i < len; i++) {
		if (a[i] != b[i])return 0;
	}
	return 1;
}

void printTwoNum() {
	cout<<endl << "dbg Two Num succ" << endl;
	for (int i = 0; i < n; i++) {
		if (two_num_succ[i])cout << i << ":" << two_num_succ[i] << ",\t";
	}
	cout << endl<<"dbg Two Num fail" << endl;
	for (int i = 0; i < n; i++) {
		if (two_num_fail[i])cout << i << ":" << two_num_fail[i] << ",\t";
	}cout << endl<<"----" << endl;
}

void printOneNum() {
	cout << endl << "dbg One Num succ" << endl;
	for (int i = 0; i < n; i++) {
		if (one_num_succ[i])cout << i << ":" << one_num_succ[i] << ",\t";
	}
	cout << endl << "dbg One Num fail" << endl;
	for (int i = 0; i < n; i++) {
		if (one_num_fail[i])cout << i << ":" << one_num_fail[i] << ",\t";
	}cout << endl << "----" << endl;
}

int getTwoNum(int p1[], int p2[], char err[]) {
	int ret = 0;
	// mult the matrix by p1,p2  and  err   -> store in res
	for (int i = 0; i < r; i++) {
		int sum = 0;
		int bg1 = 0, bg2 = 0;
		while (p1[bg1] != -1) {
			//idx = p1[bg1]+i
			if (err[(r - p1[bg1] + i) % r])sum++;
			bg1++;
		}
		while (p2[bg2] != -1) {
			//idx = p1[bg1]+i
			if (err[r + (r - p2[bg2] + i) % r])sum++;
			bg2++;
		}
		//res[i] = sum % 2;
		if (sum==4)ret++;
	}
	return ret;
}

int getOneNum(int p1[], int p2[], char err[]) {
	int ret = 0;
	// mult the matrix by p1,p2  and  err   -> store in res
	for (int i = 0; i < r; i++) {
		int sum = 0;
		int bg1 = 0, bg2 = 0;
		while (p1[bg1] != -1) {
			//idx = p1[bg1]+i
			if (err[(r - p1[bg1] + i) % r])sum++;
			bg1++;
		}
		while (p2[bg2] != -1) {
			//idx = p1[bg1]+i
			if (err[r + (r - p2[bg2] + i) % r])sum++;
			bg2++;
		}
		//res[i] = sum % 2;
		if (sum==3)ret++;
	}
	return ret;
}