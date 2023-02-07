#include "upcStatistical.h"
#include <iostream>
#include "fastbbf.h"

using namespace std;

int upcStatistical(char poly1[], char poly2[], int t, int r) {
	int n = 2 * r;
	int* errDis = new int[n];
	int* correctDis = new int[n];
	for (int i = 0; i < n; i++) {
		errDis[i] = correctDis[i] = 0;
	}

	int w = 71;

	int totNum = 500;
	for (int iter = 0; iter < totNum; iter++) {
		//rand polys
		getRandVec(poly1, r, w);
		getRandVec(poly2, r, w);
		int idx1[20000];
		int idx2[20000];
		int i1 = 0, i2 = 0;
		for (int i = 0; i < r; i++) {
			if (poly1[i])idx1[i1++] = i;
			if (poly2[i])idx2[i2++] = i;
		}
		idx1[i1] = idx2[i2] = -1;

		char* tmpErr = new char[n];
		char* tmpSynd = new char[r];
		getRandVec(tmpSynd, r, 0);
		getRandVec(tmpErr, n, t);
		//matrixMult(poly1, poly2, tmpErr, tmpSynd, r);
		matrixMultFast(idx1, idx2, tmpErr, tmpSynd, r);

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

		for (int i = 0; i < r; i++) {//left column
			int tmpUpc = 0;
			int bg1 = 0;
			while (idx1[bg1] != -1) {
				// idx1[bg1]=(i-j)%r
				if (tmpSynd[(i + idx1[bg1]) % r] == 1) { tmpUpc++; }
				bg1++;
			}
			if (tmpErr[i])errDis[tmpUpc]++;
			else correctDis[tmpUpc]++;
			
		}
		for (int i = r; i < n; i++) {//right column
			int tmpUpc = 0;
			int bg2 = 0;
			while (idx2[bg2] != -1) {
				// idx1[bg1]=(i-j)%r
				if (tmpSynd[(i + idx2[bg2]) % r] == 1) { tmpUpc++; }
				bg2++;
			}
			if (tmpErr[i])errDis[tmpUpc]++;
			else correctDis[tmpUpc]++;
		}


		delete[] tmpErr;
		delete[] tmpSynd;
	}



	int sumErr = 0, sumCorrect = 0;
	for (int i = 0; i < n; i++) {
		sumErr += errDis[i];
		sumCorrect += correctDis[i];
	}
	cout << "errDistribute:" << endl;
	for (int i = 0; i < n; i++) {
		if (errDis[i])cout << i << ":" << (double)(errDis[i]+0.0)/sumErr<<","<<endl;
	}

	cout << "correctDistribute:" << endl;
	for (int i = 0; i < n; i++) {
		if (correctDis[i])cout << i << ":" << (double)(correctDis[i] + 0.0) / sumCorrect<<"," << endl;
	}

	return 0;
}