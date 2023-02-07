#include "intersections.h"
#include <iostream>
#include <algorithm>
#include "util.h"

using namespace std;

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

int dbg = 0;

int printNum(int* x, int len) {
	cout << "dbg num:" << endl;
	int* num;
	num = new int[len];
	for (int i = 0; i < len; i++)num[i] = 0;

	for (int i = 0; i < len; i++) {
		num[x[i]]++;
	}
	for (int i = 0; i < len; i++) {
		if (num[i]) {
			cout << i << " " << num[i] << endl;
		}
	}

	delete[] num;
	return 0;
}


int getIntersections(char poly1[], char poly2[], int r,int t) {
	int n = 2 * r;
	int** inters;
	inters = new int*[n];
	for (int i = 0; i < n; i++) {
		inters[i] = new int[n];
		for (int j = 0; j < n; j++)inters[i][j] = 0;
	}

	int* ones;
	ones = new int[r];

	int iter = 0;
	for (int i = 0; i < r; i++) {
		if (poly1[i] == 1)ones[iter++] = i;
	}
	for (int i = r; i < n; i++) {
		if (poly2[i-r] == 1)ones[iter++] = i;
	}

	for (int row = 0; row < r; row++) {
		for (int i = 0; i < iter; i++) {
			for (int j = i+1; j < iter; j++) {
				int x, y;
				x = ones[i];
				y = ones[j];
				if (x < r)x = (x + row) % r;
				else x = (x + row) % r + r;
				if (y < r)y = (y + row) % r;
				else y = (y + row) % r + r;
				
				inters[x][y]++;
			}
		}
	}

	if (dbg) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				cout << inters[i][j] << " ";
			}
			cout << endl;
		}
	}

	cout << "dbg inters num" << endl;
	printNum(&inters[0][0], n);

	//统计列与列交集情况
	int *nums;
	nums = new int[n];
	for (int i = 0; i < n; i++)nums[i] = 0;

	int *maxi_nums;
	maxi_nums = new int[n];
	for (int i = 0; i < n; i++)maxi_nums[i] = 0;

	int zeros=0;

	// simulate real t errors
	char* tvec;
	tvec = new char[n];
	getRandVec(tvec, n, t);

	int* errNum;
	errNum = new int[n];
	for (int i = 0; i < n; i++)errNum[i] = 0;




	int max_tmp = 0;
	for (int i = 0; i < n; i++) {
		int tmp = 0;	//tmp统计最大值改为tmp统计前k大的和
		for (int j = 0; j < n; j++) {
			if(j==i)continue;
			if (tmp < inters[i][j])tmp = inters[i][j];
			if (inters[i][j] == 0)zeros++;
			if (tvec[j] == 1)errNum[i] += inters[i][j];
		}
		
		maxi_nums[tmp]++;

		tmp = 0;

		sort(inters[i], inters[i]+n);
		int t = 60;
		for (int j = n - 1; j >= n - t; j--) {
			tmp += inters[i][j];
		}
		nums[tmp]++;
		if (max_tmp < tmp)max_tmp = tmp;
	}
	for (int i = 0; i <= max_tmp; i++)
		if(nums[i])cout << i << " " << nums[i] << endl;
	cout << "max intersect" << endl;
	for (int i = 0; i < n; i++)
		if (maxi_nums[i])cout << i << " " << maxi_nums[i] << endl;

	cout << "zeros=" << zeros << endl;

	printNum(errNum, n);


	delete[] errNum;
	delete[] tvec;
	delete[] nums;
	delete[] ones;



	for (int i = 0; i < n; i++) {
		delete [] inters[i];
	}
	delete[] inters;

	return 0;
}

int getUnionIntersections(char poly1[], char poly2[], int r, int t) {
	int n = 2 * r;
	char* tVec=new char[n];
	getRandVec(tVec, n, t);
	
	int *tIdx = new int[n];
	for (int i = 0; i < n; i++)tIdx[i] = 0;
	int idx = 0;
	for (int i = 0; i < n; i++) {
		if (tVec[i])tIdx[idx++] = i;
	}

	int joinIdx = rand() % n;

	int tot = 0;
	int joinNum = 0;
	for (int i = 0; i < r; i++) {
		int found = 0;
		for (int j = 0; j < idx; j++) {
			int getIdx = tIdx[j];
			if (getIdx < r) {
				// in poly1
				if (poly1[(i - getIdx + r) % r] == 1) {
					found = 1; break;
				}
			}
			else {
				if (poly2[(i - getIdx + n) % r] == 1) {
					found = 1; break;
				}
			}
		}
		tot += found;
		if (joinIdx < r) {
			if (poly1[(i + joinIdx + r) % r])joinNum += found;
		}
		else {
			if (poly2[(i + joinIdx) % r])joinNum += found;
		}
	}
	cout << "tot=" << tot << endl;
	cout << "joinNum=" << joinNum << endl;


	delete[] tIdx;
	delete[] tVec;
	return 0;
}