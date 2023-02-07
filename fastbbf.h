#pragma once

int blackGrayFlip2Fast();

int blackGrayFlipFast();


void BFIterFast(char tmp_s[], char y[], int black[], int gray[], int idx1[], int idx2[], int T);

void BFMaskedIterFast(char tmp_s[], char y[], int mask[], int idx1[], int idx2[], int T);

int thresholdFast(char s[]);

void matrixMultFast(int p1[], int p2[], char err[], char res[], int r);