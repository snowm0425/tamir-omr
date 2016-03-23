#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "pnmfile.h"


int nRows, nCols;
int lThickness, lSpacing;
int** M; //input matrix
int** W; //working matrix
int** I; //output matrix

int** newMatrix(int dim1,int dim2)
{
	int** M = (int**)malloc((size_t)sizeof(int*)*dim1);
	for (int i = 0; i < dim1; i++)
	{
		M[i] = (int*) malloc((size_t)(sizeof(int)*dim2) );
		for (int j = 0; j < dim2; j++)
			M[i][j] = 0;
	}
	return M;
}

void delMatrix(int** M, int dim1,int dim2)
{
	for (int i = 0; i < dim1; i++)
	{
		free(M[i]);
	}
	free(M);
}

int OtsuThreshold()
{
	int hist[256];
	for (int i = 0; i < 256; i++)
		hist[i] = 0;
	// Calculate histogram
	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			hist[M[i][j]]++;
		}
	}
	
	// Total number of pixels
	int total = nRows*nCols;
	
	float sum = 0;
	for (int t=0 ; t<256 ; t++) sum += t * hist[t];
	
	float sumB = 0;
	int wB = 0;
	int wF = 0;
	
	float varMax = 0;
	int threshold = 0;

	for (int t=0 ; t<256 ; t++)
	{
		wB += hist[t];               // Weight Background
		if (wB == 0) continue;
		wF = total - wB;                 // Weight Foreground
		if (wF == 0) break;
		sumB += (float) (t * hist[t]);
		float mB = sumB / wB;            // Mean Background
		float mF = (sum - sumB) / wF;    // Mean Foreground
		// Calculate Between Class Variance
		float varBetween = (float)wB * (float)wF * (mB - mF) * (mB - mF);
		// Check if new maximum found
		if (varBetween > varMax)
		{
			varMax = varBetween;
			threshold = t;
		}
	}
	return threshold;
}

void runThreshold(int** O, int threshold)
{
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++)
			//O[i][j] = (int)(M[i][j] > threshold);
			O[i][j] = (int)(M[i][j]);
}

void runLengthCodes(int** I, int* hist0, int* hist1)
{
	for (int i = 0; i <= nRows; i++)
	{
		hist0[i] = 0;
		hist1[i] = 0;
	}
	for (int j = 0; j < nCols; j++)
	{
		int l = 1;
		for (int i = 1; i < nRows; i++)
		{
			if (I[i][j]==I[i-1][j])
				l++;
			else
			{
				if (I[i][j]==0)
					hist1[l]++;
				else
					hist0[l]++;
				l = 1;
			}
		}
		if (I[nRows-1][j]==1)
			hist1[l]++;
		else
			hist0[l]++;
	}
}

int getMaxVector(int* v, int sz)
{
	int max = v[0];
	int idx = 0;
        for (int i = 1; i < sz; i++)
	{
		if (v[i] > max)
		{
			max = v[i];
			idx = i;
		}
	}
	return idx;
}

void extractParameters(int* hist0, int* hist1)
{
	lThickness = getMaxVector(hist0, nRows+1);
	lSpacing = getMaxVector(hist1, nRows+1);
}

int* getPattern(int& len)
{
	len = 5*lThickness+5*lSpacing;
	int *pattern = (int*)malloc(sizeof(int)*len);
	
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < lThickness; j++)
		{
			pattern[i*(lThickness+lSpacing)+j] = 1;
		}
		for (int j = 0; j < lSpacing; j++)
		{
			pattern[i*(lThickness+lSpacing)+j+lThickness] = 0;
		}
	}

	return pattern;
}

void printMatrix(int** M, int dim1, int dim2)
{
	return;
	for (int i = 0; i < dim1; i++)
	{
		printf("\n");
		for (int j = 0; j < dim2; j++)
		{
			printf("%d ",M[i][j]);
		}
	}
}

void printMatrix2File(int** M, int dim1, int dim2, const char* fName)
{
	FILE *fout = fopen(fName,"wt");
	//fprintf(fout,"%d %d\n",dim1,dim2);
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2-1; j++)
		{
			fprintf(fout,"%d ",M[i][j]);
		}
		fprintf(fout,"%d\n",M[i][dim2-1]);
	}
	fclose(fout);
}

int** getMatches(int** I)
{
	int** Matches = newMatrix(nRows, nCols);
	int delta = (lThickness < lSpacing ? lThickness : lSpacing-2);
	int delta2= (delta+1)/2;
	//int delta2=0;
	int matches;

	for (int i = delta2; i < nRows-5*(lThickness+lSpacing); i ++)
	{
		for (int j = 0; j < nCols; j++)
		{
			matches = 0;
			for (int k = 0; k < 5; k++)
			{
				for (int l = -delta2; l < lThickness+delta2;l++)
					if (I[i+k*(lThickness+lSpacing)+l][j]==1)
					{
						matches++;
						break;
					}
				if (k < 4)
				{
					for (int l = 0; l < lSpacing;l++)
						if (I[i+k*(lThickness+lSpacing)+delta2+l][j]==0)
						{
							matches++;
						    break;
						}
				}
			}
			Matches[i][j] = matches + 2*I[i][j];
		}
	}

	for (int i = 0; i < delta2; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			Matches[i][j] = 0;
                }
	}

	for (int i = nRows-5*(lThickness+lSpacing); i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			Matches[i][j] = 0;
		}
	}

	return Matches;
}

void finalStablePaths(int** I, int** Matches, int lThickness, int lSpacing)
{
        int penalty = lThickness/2;//*3/2;
        int maxim, iMatch;
        int **W = newMatrix(nRows,nCols);
	int **P = newMatrix(nRows,nCols);
        //int **Out = newMatrix(nRows,nCols);

        int delta = (lThickness)/2;

        if (lThickness > lSpacing/2)
        {            
            delta = lSpacing/4;
            delta = (delta > 0 ? delta : 1 );

            //penalty = delta;
            penalty = 4*lThickness;
            //penalty = 5;
            delta = 1;
        }
        //delta = 3;
        //penalty = 3*lThickness;
        penalty = 3*10*lThickness;

        //delta = 4;
        if (lThickness > 4 && lThickness < 6 && lThickness < lSpacing/2)
        {
            penalty = lThickness*3/2;
            penalty = 3*10*lThickness;
            delta = 1;
        }

        printf("\ndelta = %d\n",delta);
        //adapt costs

        /*int line = 0;

        for (int i = 0; i < nRows; i++)
        {
            if (Matches[i][nCols-1]==-30000)
            {
                line++;
                int j = i;
                while (j>=0 && Matches[j][nCols-1]!=-20000)
                {
                    P[j][nCols-1]=line;
                    j--;
                }
                while (j<nRows && Matches[j][nCols-1]!=-20000)
                {
                    P[j][nCols-1]=line;
                    j++;
                }
            }
        }*/

        //for (int i = lThickness; i < nRows-lThickness-lThickness; i++)
        for (int i = 3*lSpacing; i < nRows-lThickness-3*lSpacing; i++)
        {
            if (Matches[i][nCols-1] == -30000)
                    P[i][nCols-1] = 1;

                for (int j = 0; j < nCols; j++)
                {
                        if (Matches[i][j]==-20000)
                                continue;

                        Matches[i][j] = 0;

                        //for (int k = -delta; k < lThickness + delta; k ++)
                        //for (int k = -lSpacing/8; k < lThickness + lSpacing/8; k ++)
                        //for (int k = -lThickness; k < lThickness+lThickness; k ++)
                        for (int k = 0; k < lThickness; k ++)
                        {
                                if (k < 0 || k >= lThickness)
                                        //Matches[i][j] += (I[i+k][j]==0)+(I[i+k-lThickness-lSpacing][j]==0)+(I[i+k+lThickness+lSpacing][j]==0);
                                        Matches[i][j] += (I[i+k][j]==0);//+(I[i+k-lThickness-lSpacing][j]==0)+(I[i+k+lThickness+lSpacing][j]==0);
                                else
                                        Matches[i][j] += 40*(I[i+k][j]==1)+10*(I[i+k-lThickness-lSpacing][j]==1)+10*(I[i+k+lThickness+lSpacing][j]==1);
                        }
                }
        }

        //run forward matching
	for (int j=1; j < nCols; j++)
	{
                for (int i=2; i < nRows-lThickness-lSpacing; i++)
		{
			maxim = -1000000;
                        for (int k = -2; k <= 2; k++)
                    //for (int k = -1; k <= 1; k++)
			{
				if (maxim < W[i+k][j-1] - abs(k)*penalty && Matches[i+k][j-1] >= 0)
				{
					maxim = W[i+k][j-1] - abs(k)*penalty;
				}
			}
			W[i][j] = Matches[i][j] + maxim;
		}
	}


       // for (int i = 0; i < nRows; i++)
       //     for (int j = 0; j < nCols; j ++)
       //         I[i][j] = 0;

	// trace back
        int line = 0;
	for (int i=2; i < nRows-lThickness-lSpacing-2; i++)
	{
		if (Matches[i][nCols-1] <=-10000 || P[i][nCols-1]==0)
			continue;
                /*line = P[i][nCols-1];
                int idmax = i, l=0;

                while (P[i+l][nCols-1]==line)
                {
                    P[i+l][nCols-1]=0;
                    if (W[i+l][nCols-1] >= W[idmax][nCols-1])
                        idmax = i+l;
                }

                //i = idmax;
*/
		iMatch = i;
		P[iMatch][nCols-1] = 1;
                line++;
		printf("Line %d starts at row %d and has evidence = %d\n",line,iMatch,W[iMatch][nCols-1]);

		for (int j = nCols-1; j >= 1; j--)
		{
                        int t1, t2, iM = iMatch;
                        int maxLen = lThickness+2*delta+lThickness%2;

			t1 = 0;
                        while (t1 < 2*maxLen && iM+t1 < nRows && I[iM+t1][j]==I[iM][j])
			{
				t1++;
			}

                        t2 = 0;
                        while (t2 < 2*maxLen && iM-t2 >= 0 && I[iM-t2][j]==I[iM][j])
			{
				t2++;
			}			

			if (I[iM][j]==1)
			{                            
                                if (t1+t2-1 <= maxLen)
				{
                                        for (int t = -t2; t < t1; t++)					
						I[iM+t][j]=0;
				}
			}
			else
                        {
                                if (t1 <= lThickness+delta+lThickness%2)
				{
                                        while (t1 < 2*maxLen && iM+t1 < nRows && I[iM+t1][j]==1)
						t1++;
                                        if (t1-1 <= lThickness+delta+lThickness%2)
					{
                                                while (t1 >= 0)
						{
                                                        I[iM+t1][j] = 0;
                                                        t1--;
						}
					}

                                }

                                if (t2-1 <= delta+lThickness%2)
                                {
                                        while (t2 < 2*maxLen && iM-t2 >= 0 && I[iM-t2][j]==1)
                                                t2++;
                                        if (t2-1 <= delta+lThickness%2)
                                        {
                                                while (t2 >= 0)
                                                {
                                                        I[iM-t2][j] = 0;
                                                        t2--;
                                                }
                                        }
                                }

                                if (I[iM-delta-lThickness%2][j] == 0 && I[iM+lThickness+delta+lThickness%2][j]==0)
                                    for (int t = -delta-lThickness%2; t < lThickness+delta+lThickness%2; t++)
                                        I[iM+t][j] = 0;

                                /*if (t1 <= lThickness+1)
                                {
                                        while (t1 < 2*lThickness && iM+t1 < nRows && I[iM+t1][j]==1)
                                                t1++;
                                        if (t1 <= lThickness+2)
                                        {
                                                while (t2 >= 0)
                                                {
                                                        I[iM+t2][j] = 0;
                                                        t2--;
                                                }
                                        }

                                }

                                if (t2 <= lThickness+1)
				{
                                        while (t2 < 2*lThickness && iM-t2 >= 0 && I[iM-t2][j]==1)
						t2++;
                                        if (t2 <= lThickness+2)
					{
						while (t2 >= 0)
						{
							I[iM-t2][j] = 0;
							t2--;
						}
					}
                                }*/
                        }
/*
                        if (I[iMatch-delta-1][j]==0 && I[iMatch+lThickness+delta][j]==0)
                            for (int k=iMatch-delta-1; k <= iMatch+lThickness+delta; k++)
                                I[k][j]=0;
*/
                   //     for (int p = 0; p < lThickness; p++)
                   //         I[iM+p][j]=1;

                        //for (int k = -1; k <= 1; k++)
                        for (int k = -2; k <= 2; k++)
			{
				if (W[iMatch][j] == Matches[iMatch][j] + W[iMatch+k][j-1] - abs(k)*penalty)
				{
					iMatch = iMatch+k;
					break;
				}
			}
                        if (Matches[iMatch][j-1] < 0 || W[iMatch][j-1]<0)
			{
				printf("\nProblem! Matches[%d][%d]=%d\n",iMatch,j-1,Matches[iMatch][j-1]);
				break;
			}
		    P[iMatch][j-1] = 1;
			
			if (j==1)
			{
				printf("\nTraced back (after forward) to %d row\n",iMatch);
			}

		}
	}


	delMatrix(W,nRows,nCols);
	delMatrix(P,nRows,nCols);
}

void runStaveMatching()
{        
        I = newMatrix(nRows, nCols);

        int threshold = 0;//OtsuThreshold();
        //printf("\nOtsu Threshold = %d", threshold);
	//runThreshold(I,threshold+51);
	runThreshold(I,threshold+0);

	printf("\n Matrix binarized");
	printMatrix(I,nRows, nCols);


        int* hist0 = new int[nRows+1];        
        int* hist1 = new int[nRows+1];

	runLengthCodes(I, hist0, hist1);	
	extractParameters(hist1, hist0);

        int l1 = hist1[lThickness]+(hist1[lThickness-1] > hist1[lThickness+1] ? hist1[lThickness-1] : hist1[lThickness+1]);
	int l0 = hist0[lSpacing]+(hist0[lSpacing-1] > hist0[lSpacing+1] ? hist0[lSpacing-1] : hist0[lSpacing+1]);
	float r1 = hist1[lThickness]/(float)(nCols*(nRows/(lSpacing+lThickness)));
	float r0 = hist0[lSpacing]/(float)(nCols*(nRows/(lSpacing+lThickness)));
	float lr1 = l1/(float)(nCols*(nRows/(lSpacing+lThickness)));
	float lr0 = l0/(float)(nCols*(nRows/(lSpacing+lThickness)));

        delete hist0;
        delete hist1;

	printf("\nThickness = %d, ratio = %.4f, ratio = %.4f",lThickness,r1, lr1);
	printf("\nSpacing = %d, ratio = %.4f, ratio = %.4f",lSpacing,r0, lr0);
	printf("\nRatio = %.4f",r1+r0);

        int lPattern = 5*(lThickness+lSpacing);

        printf("\nlPattern : %d, nRows=%d nCols=%d ", lPattern, nRows, nCols);

	int **W = newMatrix(nRows, nCols);
	int **P = newMatrix(nRows, nCols);

	int **W2 = newMatrix(nRows, nCols);
	int **Matches = getMatches(I);

	int matches, maxim;
	int iMatch;

        int penalty = 13;

        int ok = 1;

        int min_evidence = 0;

       // for (int i = 0; i < nRows; i++)
       //     for (int j = 0; j < nCols; j++)
       //         I[i][j]=0;


	while (ok)
	{
		printf("\nRound\n\n");

// forward matching
	for (int j=1; j < nCols; j++)
	{
		for (int i=2; i < nRows-lPattern+1; i++)
		{
			W[i][j] = 0;
			if (Matches[i][j] <= -10000)
				continue;
			maxim = 0;
			for (int k = -1; k <= 1; k++)
				//for (int k = -2; k <= 2; k++)
			{
				if (maxim < W[i+k][j-1] - abs(k)*penalty)
				{
					maxim = W[i+k][j-1] - abs(k)*penalty;
				}
			}
			W[i][j] = Matches[i][j] + maxim;
		}
	}
	printf("\n Matrix forward");
	printMatrix(W,nRows, nCols);

	printf("\n Matrix matches");
	printMatrix(Matches,nRows, nCols);

// tracing back
	for (int i=2; i < nRows-lPattern+1; i++)
	{
		iMatch = i;
                if (Matches[i][nCols-1] <=-10000 || W[iMatch][nCols-1] <= nCols*4+nCols*3/2 || W[iMatch][nCols-1] <= min_evidence-min_evidence/5)// || W[iMatch][nCols-1]<6900)
			continue;

                if (W[iMatch][nCols-1] > min_evidence)
                {
                    min_evidence = W[iMatch][nCols-1];
                    //printf("\n min evidence= %d\n", min_evidence);
                }


		P[iMatch][nCols-1] = 1;
		for (int j = nCols-1; j >= 1; j--)
		{
			for (int k = -1; k <= 1; k++)
			//for (int k = -2; k <= 2; k++)
			{
				if (W[iMatch][j] == Matches[iMatch][j] + W[iMatch+k][j-1] - abs(k)*penalty)
				{
					iMatch = iMatch+k;
					break;
				}
			}
			if (P[iMatch][j-1] == 1 || Matches[iMatch][j-1] <= -10000)
			{
				break;
			}
		    P[iMatch][j-1] = 1;
			if (j==1)
			{
				printf("\nTraced back (after forward) to %d row\n",iMatch);
                                i += lPattern/2;
			}
		}
	}
	printf("\n Matrix trace");
	printMatrix(P,nRows, nCols);
// backward matching
        for (int j=nCols-2; j >= 0; j--)
        {
                for (int i=2; i < nRows-lPattern+1; i++)
                {
                        W2[i][j] = 0;
                        if (Matches[i][j] <= -10000)
                                continue;
                        maxim = 0;
                        for (int k = -1; k <= 1; k++)
                                //for (int k = -2; k <= 2; k++)
                        {
                                if (maxim < W2[i+k][j+1] - abs(k)*penalty)
                                {
                                        maxim = W2[i+k][j+1] - abs(k)*penalty;
                                }
                        }
                        W2[i][j] = Matches[i][j] + maxim;
                }
        }
        printf("\n Matrix back");
        printMatrix(W2,nRows, nCols);
// tracing
        for (int j = 1; j < nCols; j++)
        {
                for (int i = 0; i < nRows; i++)
                {
                        P[i][j] = 0;
                }
        }
        int line = 0;
        for (int i = 10; i < nRows; i++)
        {
                if (P[i][0]==0 || Matches[i][0]<=-10000 || W2[i][0] <= nCols*4+nCols*3/2)
                        continue;

                // trace back

                int iM = i, stopped = 0;
                for (int j = 1; j < nCols; j++)
                {
                        //for (int k = -2; k <= 2; k++)
                        for (int k = -1; k <= 1; k++)
                        {
                                if (W2[iM][j-1] == Matches[iM+k][j] + W2[iM+k][j] - abs(k)*penalty)
                                {
                                        iM = iM+k;
                                        break;
                                }
                        }
                        if (P[iM][j] == 1)
                        {
                            stopped = 1;
                                break;
                        }
                }
                if (stopped || abs(W2[i][0]-W[iM][nCols-1]) > 100)
                {
                    printf("\nLine candidate starting at row %d ending at row %d has %d!=%d, it is not stable\n", i, iM, W2[i][0], W[iM][nCols-1]);
                    continue;
                }

                line ++;
                printf("\nLine %d starts at row %d and has W = %d",line,i,W2[i][0]);
                iMatch = i;
                P[iMatch][0] = 1;
                for (int j = 1; j < nCols; j++)
                {
                        int dim1 = (iMatch-lPattern/2 > 0 ? iMatch-lPattern/2:0);
                        int dim2 = (iMatch+lPattern < nRows ? iMatch+lPattern:nRows);
                        for (int k = dim1; k < dim2; k++)
                        {
                                if (Matches[k][j-1]>-10000)
                                        Matches[k][j-1] = -10000;
                        }

                        //remove stave ?
                        for (int k = 0; k < 5; k++)
                        {
                                //Matches[iMatch+lThickness/2+k*(lThickness+lSpacing)-(lThickness+lSpacing)/2][j] = - 20000;
                                //Matches[iMatch+lThickness/2+k*(lThickness+lSpacing)-(lThickness+lSpacing)/2-1][j] = - 20000;
                                //Matches[iMatch+lThickness/2+k*(lThickness+lSpacing)-(lThickness+lSpacing)/2-2][j] = - 20000;

                                Matches[iMatch+lThickness/40+k*(lThickness+lSpacing)-(lThickness+lSpacing)/2][j] = - 20000;
                                Matches[iMatch+lThickness/40+k*(lThickness+lSpacing)-(lThickness+lSpacing)/2-1][j] = - 20000;
                                Matches[iMatch+lThickness/40+k*(lThickness+lSpacing)-(lThickness+lSpacing)/2-2][j] = - 20000;

                                if (k==4)
                                {
                                        //Matches[iMatch+lThickness/2+k*(lThickness+lSpacing)+(lThickness+lSpacing)/2][j] = - 20000;
                                        //Matches[iMatch+lThickness/2+k*(lThickness+lSpacing)+(lThickness+lSpacing)/2+1][j] = - 20000;
                                        //Matches[iMatch+lThickness/2+k*(lThickness+lSpacing)+(lThickness+lSpacing)/2+2][j] = - 20000;
                                        Matches[iMatch+lThickness/40+k*(lThickness+lSpacing)+(lThickness+lSpacing)/2][j] = - 20000;
                                        Matches[iMatch+lThickness/40+k*(lThickness+lSpacing)+(lThickness+lSpacing)/2+1][j] = - 20000;
                                        Matches[iMatch+lThickness/40+k*(lThickness+lSpacing)+(lThickness+lSpacing)/2+2][j] = - 20000;

                                  //      I[iMatch+lThickness/2+k*(lThickness+lSpacing)+(lThickness+lSpacing)/2][j] = 1;
                                  //      I[iMatch+lThickness/2+k*(lThickness+lSpacing)+(lThickness+lSpacing)/2+1][j] = 1;

                                }                                

                                //I[iMatch+lThickness/2+k*(lThickness+lSpacing)-(lThickness+lSpacing)/2][j] = 1;
                                //I[iMatch+lThickness/2+k*(lThickness+lSpacing)-(lThickness+lSpacing)/2-1][j] = 1;
                                //I[iMatch+k*(lThickness+lSpacing)][j] = 1;
                                //I[iMatch+k*(lThickness+lSpacing)+1][j] = 1;

                        }


                        //for (int k = -2; k <= 2; k++)
                        for (int k = -1; k <= 1; k++)
                        {
                                if (W2[iMatch][j-1] == Matches[iMatch+k][j] + W2[iMatch+k][j] - abs(k)*penalty)
                                {
                                        iMatch = iMatch+k;
                                        break;
                                }
                        }
                        if (P[iMatch][j] == 1)
                        {
                                break;
                        }
                        P[iMatch][j] = 1;
                }
                int dim1 = (iMatch-lPattern/2 > 0 ? iMatch-lPattern/2:0);
                int dim2 = (iMatch+lPattern < nRows ? iMatch+lPattern:nRows);
                for (int k = dim1; k < dim2; k++)
                {
                        if (Matches[k][nCols-1]>-10000)
                                Matches[k][nCols-1] = -10000;
                }

                                        //remove stave ?
                        for (int k = 0; k < 5; k++)
                        {
                                //Matches[iMatch+lThickness/2+k*(lThickness+lSpacing)-(lThickness+lSpacing)/2][nCols-1] = - 20000;
                                //Matches[iMatch+lThickness/2+k*(lThickness+lSpacing)-(lThickness+lSpacing)/2-1][nCols-1] = - 20000;
                                //Matches[iMatch+lThickness/2+k*(lThickness+lSpacing)-(lThickness+lSpacing)/2-2][nCols-1] = - 20000;

                                Matches[iMatch+lThickness/40+k*(lThickness+lSpacing)-(lThickness+lSpacing)/2][nCols-1] = - 20000;
                                Matches[iMatch+lThickness/40+k*(lThickness+lSpacing)-(lThickness+lSpacing)/2-1][nCols-1] = - 20000;
                                //Matches[iMatch+lThickness/40+k*(lThickness+lSpacing)-(lThickness+lSpacing)/2-2][nCols-1] = - 20000;

                                Matches[iMatch+lThickness/2+k*(lThickness+lSpacing)][nCols-1] = -30000;
                                Matches[i+lThickness/2+k*(lThickness+lSpacing)][0] = -30000;
                                if (k == 4)
                                {
                                        //Matches[iMatch+lThickness/2+k*(lThickness+lSpacing)+(lThickness+lSpacing)/2][nCols-1] = - 20000;
                                        //Matches[iMatch+lThickness/2+k*(lThickness+lSpacing)+(lThickness+lSpacing)/2+1][nCols-1] = - 20000;
                                        //Matches[iMatch+lThickness/2+k*(lThickness+lSpacing)+(lThickness+lSpacing)/2+2][nCols-1] = - 20000;

                                        Matches[iMatch+lThickness/40+k*(lThickness+lSpacing)+(lThickness+lSpacing)/2][nCols-1] = - 20000;
                                        Matches[iMatch+lThickness/40+k*(lThickness+lSpacing)+(lThickness+lSpacing)/2+1][nCols-1] = - 20000;
                                        //Matches[iMatch+lThickness/40+k*(lThickness+lSpacing)+(lThickness+lSpacing)/2+2][nCols-1] = - 20000;

                                }
                        }

                printf("\nLine %d stops at row %d and has W = %d",line,iMatch,W[iMatch][nCols-1]);
        }
        printf("\nP final");
	printMatrix(P,nRows, nCols);
	ok = line;
	} 

        finalStablePaths(I, Matches, lThickness, lSpacing);


	delMatrix(W,nRows, nCols);
	delMatrix(W2,nRows, nCols);
	delMatrix(Matches, nRows, nCols);
	delMatrix(P, nRows, nCols);

}

int main(int argc, char **argv)
{
    image<uchar> *in;

    if (argc != 3)
    {
        printf("usage: %s in(pgm) out(pgm)\n",argv[0]);
        exit(1);
    }
    time_t seconds = time(NULL);

    // load input

    //in = loadPGM(argv[1]);
    in = loadPBM(argv[1]);

    nCols   = in->width();
    nRows = in->height();

    M = (int**)malloc((size_t)sizeof(int*)*nRows);
    for (int i = 0; i < nRows; i++)
    {
            M[i] = (int*) malloc((size_t)(sizeof(int)*nCols) );
            for (int j = 0; j < nCols; j++)
                M[i][j] = (uchar)(imRef(in, j, i)==0);
    }

    delete in;

    printf("\n Matrix original");
    printMatrix(M,nRows, nCols);
    runStaveMatching();
    delMatrix(M,nRows, nCols);

    image<uchar> *out = new image<uchar>(nCols, nRows);

    for (int i = 0; i < nRows; i++)
    {
        for (int j = 0; j < nCols; j++)
                imRef(out, j, i)= (uchar)(I[i][j]==0);
    }

    // save output
    //savePGM(out, argv[2])
    savePBM(out, argv[2]);

    delete out;
    delMatrix(I,nRows, nCols);

    time_t seconds2 = time(NULL);
    printf ("\ntime - %.4f seconds\n", float(seconds2-seconds));
    printf ("Thickness:%d, space:%d\n",lThickness,lSpacing);
    return 0;
}
