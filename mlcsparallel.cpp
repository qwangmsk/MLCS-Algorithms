/* mlcsparallel.cpp -- Implementation of Quick-DP
 * 
 * Compile
 * 
 * g++ -pthread mlcsparallel.cpp -o dc.exe
 * 
 * 
 * Copyright (c) 2009.
 * 
 * All rights reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software, to deal in the Software without restriction, including without
 * limitation the rights to use, copy, modify, merge, and distribute of the
 * Software, provided that the above copyright notice(s) and this permission notice
 * appear in all copies of the Software and that both the above copyright notice(s)
 * and this permission notice appear in supporting documentation.
 *
 * USE ENTIRELY AT YOUR OWN RISK.
 *
 * Citation
 * [1] Q. Wang, D. Korkin, and Y. Shang, "A Fast Multiple Longest Common Subsequence
 *      (MLCS) Algorithm," IEEE Transactions on Knowledge and Data Engineering, 23(3):321-34.
 * [2] Q. Wang, D. Korkin, and Y. Shang, "Efficient Dominant Point Algorithms for the
 *	Multiple Longest Common Subsequence (MLCS) Problem," Intl. Joint Conf. on
 *	    Artificial Intelligence (IJCAI), 2009. 
 *
 * Contact: 
 *
 *	Qingguo Wang, qwp4b@mail.missouri.edu
 *
 */
#ifdef MY_WIN_THREAD
#include "test.h"
#include <conio.h>
#else

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#endif

#include <vector>


/*-------- Constants. Begin ------------*/
#define		DIVIDEANDCONQUER
#define  	MAXPOSSIBLESTR		10000
#define  	MAXFILENAMELEN		128
#define  	THREADTHRESHOLD		40

char	 	ABC_str_20[] 	= 	{'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
char	 	ABC_str_4[]  	= 	{'A','C','G','T'};
char	 	ABC_idx_20[] 	= 	{ 0,0,1,2,3,4,5,6,7,7,8,9,10,11,11,12,13,14,15,16,16,17,18,18,19,19};
char	 	ABC_idx_4[] 	= 	{ 0,0,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3};

int	  		MAXDSIZE	    =   4096;
/*-------- Constants. End --------*/

/*-------- Variables. Begin ------------*/
int			MAXSTRSIZE;
int			MAXSTR;
int			MAXCPU;
int			ABC_len;
char*		ABC_str;
char*		ABC_idx;

char**		data;
int*		str_lngth;
int			cur_n_str = 0;

char		data_file_fullname[MAXFILENAMELEN];

int			cur_k;
bool**		f = 0;	

int*					M_size;
std::vector <int**>		M;
std::vector <int**>		cur_T;
std::vector <int**>		cur_D;
std::vector <int**>		buffer;
std::vector <int>		size_cur_D;

pthread_mutex_t			mutex = PTHREAD_MUTEX_INITIALIZER;
int						threadNum;
/*-------- Variables. End ------------*/


/*-------- Template. Begin ------------*/
template <typename T> 
T** Allocate2DArray(int nRows, int nCols)
{
    //(step 1) allocate memory for array of elements of column
    T **ppi = new T*[nRows];

    //(step 2) allocate memory for array of elements of each row
    T *curPtr = new T [nRows * nCols];

    // Now point the pointers in the right place
    for( int i = 0; i < nRows; ++i)
    {
        *(ppi + i) = curPtr;
         curPtr += nCols;
    }
    return ppi;
}

template <typename T> 
void Free2DArray(T** Array)
{
    delete [] *Array;
    delete [] Array;
}

/*-------- Template. End ------------*/

/*------- Quicksort. Start ----------*/

int inline ComparePoint(int left[], int right[], int dim)
{
	for(int i=0; i<=dim; i++){
		int result = left[i]-right[i];

		if (result)
			return result;
	}
	return 0;
}
int PartitionQuickSort(int**  list, int* indices, int dim, int left, int right, int pivotIndex)
{
	//swap list[pivotIndex] and list[right]  // Move pivot to end
	int temp = indices[pivotIndex];
	indices[pivotIndex] = indices[right];
	indices[right] = temp;

	//storeIndex := left
	int storeIndex = left;
	
	//for i from left to right-1
	for(int i= left; i<=right-1; i++)
	{
		//if list[i] < pivotValue
		if (ComparePoint(list[indices[i]], list[indices[right]], dim) < 0) //list[indices[i]][dim] < list[indices[right]][dim]
		{
			//swap list[storeIndex] and list[i]
			temp = indices[storeIndex];
			indices[storeIndex] = indices[i];
			indices[i] = temp;

			//storeIndex := storeIndex + 1
			storeIndex += 1;
		}
	}
	//swap list[right] and list[storeIndex]  // Move pivot to its final place
	temp = indices[storeIndex];
	indices[storeIndex] = indices[right];
	indices[right] = temp;

	return storeIndex;
}

void QuickSort(int**  list, int* indices, int dim, int left, int right)
{
	 if (right > left)
	 {
		 int pivotIndex = left;
		 int pivotNewIndex = PartitionQuickSort(list, indices, dim, left, right, pivotIndex);
		 QuickSort(list, indices, dim, left, pivotNewIndex - 1);
		 QuickSort(list, indices, dim, pivotNewIndex + 1, right);
	 }
}


/*-------- Quicksort. end -----------*/

/*------- PartitionConquer. start ---------*/
int PartitionConquer(int**  list, int* indices, int dim, int left, int right, int pivotValue)
{
	//storeIndex := left
	int storeIndex = left;
	
	//for i from left to right-1
	for(int i= left; i<=right; i++)
	{
		//if list[i] < pivotValue
		if (list[indices[i]][dim] < pivotValue)
		{
			//swap list[storeIndex] and list[i]
			int temp = indices[storeIndex];
			indices[storeIndex] = indices[i];
			indices[i] = temp;

			//storeIndex := storeIndex + 1
			storeIndex += 1;
		}
	}

	return storeIndex;
}



class CFindParentThread {
public:
	pthread_t m_thread;

	int  i_from, i_to;	

	int* par_size;	
	int***   parent;

	int** tmp_par; 
	bool* tmp_f; 

	CFindParentThread() {
		par_size = new int[ABC_len];

		tmp_par = Allocate2DArray<int>(ABC_len, MAXSTR);
		tmp_f = new bool[ABC_len];
	}
	~CFindParentThread(){
		delete[] par_size;

		Free2DArray<int>(tmp_par);
		delete[] tmp_f;
	}
		
	void Initialize(int from, int to){
		i_from = from;
		i_to = to;
		for (int i=0; i<ABC_len; i++)
			par_size[i] = 0;			
	}

	void FindParents(int p[], int cur_MAXSTR){
		int i,j,k;
		int result = 0;
		for(i=0; i<ABC_len; i++)
		{ 
			 tmp_f[i] = false;

			 j=0;
			 while (j<cur_MAXSTR)
			 {
				 //k = cur_T[j]->FindParents(p[j], i);
				// k = FindParents(int p, int Str_idx){
				//	return D[p+1][Str_idx];
				//}

				 k = cur_T[j][p[j]+1][i];

				 if  (k!=str_lngth[j]+1)
				 {
					 tmp_par[i][j]= k;// ??  Not +1 ???
				 } /* if */
				 else
				 {
					 break;
				 }//else 
				 j++; 
			 }// while
			 //printf("\n");

			 if (j==cur_MAXSTR)
			 {
				 tmp_f[i] = true;
				 result++;
			 }
		}//for i

		if (result>1)
		{
			 for(i=0; i<ABC_len-1; i++)
			 {
				 if (!tmp_f[i]) continue;
				 
				 int *Si = tmp_par[i];
				 for(j=i+1; j<ABC_len; j++)
				 {
					 if (!tmp_f[j]) continue;
					  
					 int *Sj = tmp_par[j];
					 int fl0=0;
					 int fl1=0;
					 int fl2=0; 
					 for(k=0; k<cur_MAXSTR; k++)
					 {
						 if (Sj[k] == Si[k])
						 {
							 fl0++;
							 fl1++;
							 fl2++;
						 }//if
						 else if (Sj[k] < Si[k])
							 fl1++;
						 else //if (Sj[k] >= Si[k])
							 fl2++;
					 }//for k

					 if (fl0 == cur_MAXSTR) 
						tmp_f[j]=false;
					 else
					 {  
						 if (fl1==cur_MAXSTR)
							 tmp_f[i] = 0;
						 else if (fl2==cur_MAXSTR)
							 tmp_f[j] = 0;
					 } /* else */
				 } /* for j  */
			 }/* for i */
		}

		for(i=0; i<ABC_len; i++)
		{ 
			int** par = parent[i];
			if(tmp_f[i])
			{
				memcpy(par[par_size[i]++], tmp_par[i], cur_n_str*sizeof(int));
			}
		}
	}

	static void* StartThread(void *ptr)
	{
		CFindParentThread* pSlaveThread =  (CFindParentThread*)ptr;

		int i, *cur_p;
		for(i=pSlaveThread->i_from; i<pSlaveThread->i_to; i++) 
		{
			cur_p = cur_D[cur_k][i];
			pSlaveThread->FindParents(cur_p, cur_n_str);
		} /* for i */ 
		
		//Union
		pthread_mutex_lock(&mutex);

		for(i=0; i<ABC_len; i++)
		{ 
			int **source = pSlaveThread->parent[i];
			int **dest   = M[i];

			int size_cur_Dk=M_size[i];
			M_size[i] += pSlaveThread->par_size[i];

			memcpy(dest[size_cur_Dk], source[0], pSlaveThread->par_size[i]*cur_n_str*sizeof(int));
		}

		pthread_mutex_unlock(&mutex);		

		return 0;
	}
};	



class CDivConqMinThread {
public:
	pthread_t m_thread;

	int**  m_M;
	int*   m_M_size;
	bool*  f;
	int*   stack;
	int    stack_size;
	int*   isA;

	int    *q_indices, *q_indicesA, *q_indicesB;
	int	   q_dim, q_left, q_right;

	CDivConqMinThread* left_son;
	CDivConqMinThread* mid_son;
	CDivConqMinThread* right_son;

	CDivConqMinThread() {
		stack = 0;
		stack_size = 0;
		left_son = 0;
		right_son = 0;
	}
		
	void InitThread(bool *ff, int**  MM, int* MM_size, int*  sstack, int* iisA)
	{
		f = ff;
		m_M = MM;
		m_M_size = MM_size;
		stack = sstack;
		isA = iisA;
	}

	void InitThread(CDivConqMinThread* parent, int*  sstack, int* iisA)
	{
		f        = parent->f;
		m_M      = parent->m_M;
		m_M_size = parent->m_M_size;
		stack    = sstack;
		isA      = iisA;
	}

	int* NEW_MEM(int n, int* init_val) 
	{
		if (stack_size+n>MAXDSIZE*MAXSTR)
		{
			printf("memory problem!\n");
			exit(-1);
		}

		if (init_val){
			memcpy(stack+stack_size, init_val, n*sizeof(int));
		}

		int temp = stack_size;
		stack_size += n;
		return stack+temp;
	}

	void RELEASE_MEM(int n) 
	{
		stack_size-=n;
	}


	/*------- PartitionDivide. start ---------*/
	int PartitionDivide(int* indices, int dim, int left, int right, int& pivotIndex)
	{
		//pivotValue := m_M[pivotIndex]
		int pivotValue = m_M[indices[pivotIndex]][dim];

		int* indices1  = NEW_MEM(right-left, false);
		int* indices2  = NEW_MEM(right-left, false);
		int storeIdx  = left;
		int storeIdx1  = 0;
		int storeIdx2  = 0;

		for(int i = left; i <= right; i++)
		{
			if (m_M[indices[i]][dim] < pivotValue)
				indices[storeIdx++] = indices[i];
			else if (m_M[indices[i]][dim] == pivotValue)
				indices1[storeIdx1++] = indices[i];
			else
				indices2[storeIdx2++] = indices[i];
		}

		pivotIndex = storeIdx;
		if (storeIdx1)
			memcpy(indices+storeIdx, indices1, storeIdx1*sizeof(int));
		if (storeIdx2)
			memcpy(indices+storeIdx+storeIdx1, indices2, storeIdx2*sizeof(int));

		RELEASE_MEM((right-left)+(right-left));
		return storeIdx1;
	}


	int FindFirstK(int* indices, int dim, int left, int right, int& k)
	{
		int t_left  = left;
		int t_right = right;

		while (t_right > t_left)
		{
			int pivotIndex = (t_left+t_right)/2;
			int eqlItemCnt = PartitionDivide(indices, dim, t_left, t_right, pivotIndex);

			if (pivotIndex == t_left)
			{
				if (pivotIndex+eqlItemCnt < k){
					t_left += eqlItemCnt;
					continue;
				}else{
					if (t_left == left){
						if (eqlItemCnt == right-t_left+1)
							return 7;
						else{
							k = t_left+eqlItemCnt;
							return 2;
						}
					}else{
						k = (abs(t_left-k) < abs(t_left+eqlItemCnt-k)) ? t_left : t_left+eqlItemCnt;
						return 0;
					}
				}
			}
			else if (eqlItemCnt == t_right-pivotIndex+1)
			{
				if (pivotIndex > k){
					t_right = pivotIndex-1;
					continue;
				}else{
					k = pivotIndex;
					if (t_right == right)
						return 1;
					else
						return 0;
				}
			}

			if (pivotIndex == k)
				return 0;
			else if (pivotIndex > k)
				t_right = pivotIndex-1;
			else
				t_left  = pivotIndex;
		}
		return 0;
	}
	/*-------- PartitionDivide. end ----------*/

	void Union(int* indicesA, int sizeA, int* indicesB, int sizeB, int dim)
	{
		if (dim==1)
		{
			int i_A=0, i_B=0;			
			while (i_B < sizeB && m_M[indicesA[0]][0]>m_M[indicesB[i_B]][0])
			{
				i_B++;
			}

			int minA_d1 = m_M[indicesA[0]][1];
			while(i_A < sizeA && i_B < sizeB)
			{
				while (i_A+1 < sizeA && m_M[indicesA[i_A+1]][0]<=m_M[indicesB[i_B]][0])
				{
					i_A++;
					minA_d1 = minA_d1<m_M[indicesA[i_A]][1] ? minA_d1 : m_M[indicesA[i_A]][1];
				}

				if (minA_d1 <= m_M[indicesB[i_B]][1])
					f[indicesB[i_B]] = 0;
				i_B++;
			}
			return;
		}

		int i, idx;
		if (sizeA == 1)
		{
			if (m_M[indicesB[sizeB-1]][0]<m_M[indicesA[0]][0])
				return;

			for(; dim>=0 && sizeB>0; dim--){
				idx = 0;
				for(i=0; i<sizeB; i++){
					if (m_M[indicesB[i]][dim]>=m_M[indicesA[0]][dim])
						indicesB[idx++] = indicesB[i];
				}
				sizeB = idx;
			}
			
			for(i=0; i<sizeB; i++)	f[indicesB[i]]=0;
			return;
		}
		else if (sizeB == 1)
		{
			if (m_M[indicesA[0]][0]>m_M[indicesB[0]][0])
				return;

			for(i=0; i<sizeA; i++){
				int j;
				for(j=0; j<=dim; j++){
					if (m_M[indicesA[i]][j]>m_M[indicesB[0]][j])
						break;
				}
				if (j>dim){
					f[indicesB[0]]=0;
					return;
				}
			}
			return;
		}

		// divide B
		int* indicesB_bk  = NEW_MEM(sizeB, indicesB);

		int sizeB1 = sizeB>>1;
		int retB = FindFirstK(indicesB_bk, dim, 0, sizeB-1, sizeB1);

		for(i=0; i<sizeB1; i++){
			isA[indicesB_bk[i]] = 1;
		}
		for(i=sizeB1; i<sizeB; i++){
			isA[indicesB_bk[i]] = 0;
		}
		idx = 0;
		for(i=0; i<sizeB; i++){
			if (isA[indicesB[i]])
				indicesB_bk[idx++] = indicesB[i];
		}
		for(i=0; i<sizeB; i++){
			if (!isA[indicesB[i]])
				indicesB_bk[idx++] = indicesB[i];
		}


		// divide A
		int* indicesA_bk  = NEW_MEM(sizeA, indicesA);

		int sizeA1 = PartitionConquer(m_M, indicesA_bk, dim, 0, sizeA-1, m_M[indicesB_bk[sizeB1]][dim]);


		for(i=0; i<sizeA1; i++){
			isA[indicesA_bk[i]] = 1;
		}
		for(i=sizeA1; i<sizeA; i++){
			isA[indicesA_bk[i]] = 0;
		}
		idx = 0;
		for(i=0; i<sizeA; i++){
			if (isA[indicesA[i]])
				indicesA_bk[idx++] = indicesA[i];
		}
		for(i=0; i<sizeA; i++){
			if (!isA[indicesA[i]])
				indicesA_bk[idx++] = indicesA[i];
		}

		if (left_son && sizeA>THREADTHRESHOLD && sizeB>THREADTHRESHOLD)
		{
			if (sizeA1 && sizeB1)
			{
				//Union(indicesA_bk, sizeA1, indicesB_bk, sizeB1, dim);	
				left_son->stack_size = 0;
				left_son->q_indicesA = indicesA_bk;
				left_son->q_indicesB = indicesB_bk;
				left_son->q_left = sizeA1;
				left_son->q_right = sizeB1;
				left_son->q_dim = dim;
				pthread_create(&left_son->m_thread,  NULL, left_son->StartUnion,  (void*)left_son );
			}

			if (sizeA>sizeA1 && sizeB>sizeB1)
			{
				//Union(indicesA_bk+sizeA1, sizeA-sizeA1, indicesB_bk+sizeB1, sizeB-sizeB1, dim);
				right_son->stack_size = 0;
				right_son->q_indicesA = indicesA_bk+sizeA1;
				right_son->q_indicesB = indicesB_bk+sizeB1;
				right_son->q_left = sizeA-sizeA1;
				right_son->q_right = sizeB-sizeB1;
				right_son->q_dim = dim;
				pthread_create(&right_son->m_thread,  NULL, right_son->StartUnion,  (void*)right_son );
			}

			if (sizeA1 && sizeB>sizeB1)
			{
				//Union(indicesA_bk, sizeA1, indicesB_bk+sizeB1, sizeB-sizeB1, dim-1);
				mid_son->stack_size = 0;
				mid_son->q_indicesA = indicesA_bk;
				mid_son->q_indicesB = indicesB_bk+sizeB1;
				mid_son->q_left = sizeA1;
				mid_son->q_right = sizeB-sizeB1;
				mid_son->q_dim = dim-1;
				pthread_create(&mid_son->m_thread,  NULL, mid_son->StartUnion,  (void*)mid_son );
			}

			if (sizeA1 && sizeB1)
				pthread_join(left_son->m_thread,  NULL);
			if (sizeA>sizeA1 && sizeB>sizeB1)
				pthread_join(right_son->m_thread, NULL);
			if (sizeA1 && sizeB>sizeB1)
				pthread_join(mid_son->m_thread, NULL);
		}
		else
		{
			if (sizeA1 && sizeB1)
				Union(indicesA_bk, sizeA1, indicesB_bk, sizeB1, dim);	
			if (sizeA>sizeA1 && sizeB>sizeB1)
				Union(indicesA_bk+sizeA1, sizeA-sizeA1, indicesB_bk+sizeB1, sizeB-sizeB1, dim);
			if (sizeA1 && sizeB>sizeB1)
				Union(indicesA_bk, sizeA1, indicesB_bk+sizeB1, sizeB-sizeB1, dim-1);
		}

		RELEASE_MEM(sizeB+sizeA);
	}

	void Divide(int* indices, int dim, int size)
	{
		//Termination condition
		int i, idx;
		if (size == 1){
			f[indices[0]] = 1;
			return;
		}else if (dim==1){
			int minY=m_M[indices[0]][dim];
			f[indices[0]] = 1;
			for(i=1; i<size; i++){
				if (m_M[indices[i]][dim] < minY){
					f[indices[i]] = 1;
					minY = m_M[indices[i]][dim];
				}else{
					f[indices[i]] = 0;
				}
			}
			return;
		}

		int* indices_bk  = NEW_MEM(size, indices);

		int sizeA = size>>1;
		int ret   = FindFirstK(indices_bk, dim, 0, size-1, sizeA);

		if (ret & 4) 
		{
			Divide(indices, dim-1, size);	
			RELEASE_MEM(size);
			return;
		}

		for(i=0; i<sizeA; i++){
			isA[indices_bk[i]] = 1;
		}
		for(i=sizeA; i<size; i++){
			isA[indices_bk[i]] = 0;
		}
		idx = 0;
		for(i=0; i<size; i++){
			if (isA[indices[i]])
				indices_bk[idx++] = indices[i];
		}
		for(i=0; i<size; i++){
			if (!isA[indices[i]])
				indices_bk[idx++] = indices[i];
		}

		if (!sizeA){
			Divide(indices_bk, dim-(ret == 1), size);
			return;
		}else if (size==sizeA){
			Divide(indices_bk, dim-(ret == 2), size);	
			return;
		}else if (left_son && size>THREADTHRESHOLD)
		{
			left_son->stack_size = 0;
			right_son->stack_size = 0;

			left_son->NEW_MEM(sizeA, indices_bk);
			right_son->NEW_MEM(size-sizeA, indices_bk+sizeA);

			left_son->stack[sizeA] = dim-(ret == 2);
			right_son->stack[size-sizeA] = dim-(ret == 1);

			pthread_create(&left_son->m_thread,  NULL, left_son->StartDivide,  (void*)left_son );
			pthread_create(&right_son->m_thread, NULL, right_son->StartDivide, (void*)right_son);
			pthread_join(left_son->m_thread,  NULL);
			pthread_join(right_son->m_thread, NULL);

			memcpy(indices_bk, left_son->stack, sizeA*sizeof(int));
			memcpy(indices_bk+sizeA, right_son->stack, (size-sizeA)*sizeof(int));
		}
		else
		{
			Divide(indices_bk, dim-(ret == 2), sizeA);	
			Divide(indices_bk+sizeA, dim-(ret == 1), size-sizeA);
		}

		int storeA = 0;
		for(i=0; i<sizeA; i++){
			if (f[indices_bk[i]]){
				indices_bk[storeA++] = indices_bk[i];
			}
		}

		int storeB = sizeA;
		for(i=sizeA; i<size; i++){
			if (f[indices_bk[i]]){
				indices_bk[storeB++] = indices_bk[i];
			}
		}

		if (storeA>0 && storeB>sizeA)
		{
			Union(indices_bk, storeA, indices_bk+sizeA, storeB-sizeA, dim-1);
		}

		RELEASE_MEM(size);
	}	

	static void* StartDivide(void *ptr)
	{
		CDivConqMinThread* pSlaveThread =  (CDivConqMinThread*)ptr;
		pSlaveThread->Divide(pSlaveThread->stack,pSlaveThread->stack[pSlaveThread->stack_size],pSlaveThread->stack_size);
		return 0;
	}

	static void* StartUnion(void *ptr)
	{
		CDivConqMinThread* pSlaveThread =  (CDivConqMinThread*)ptr;
		pSlaveThread->Union(pSlaveThread->q_indicesA,pSlaveThread->q_left,pSlaveThread->q_indicesB,pSlaveThread->q_right,pSlaveThread->q_dim);
		return 0;
	}

	static void* QuickSort(void *ptr)
	{
		CDivConqMinThread* pThread =  (CDivConqMinThread*)ptr;

		if (!pThread->left_son || (pThread->q_right-pThread->q_left)<THREADTHRESHOLD)
			::QuickSort(pThread->m_M, pThread->q_indices, pThread->q_dim, pThread->q_left, pThread->q_right);
		else
		{
			int pivotIndex = pThread->q_left;
			int pivotNewIndex = PartitionQuickSort(pThread->m_M, pThread->q_indices, pThread->q_dim, pThread->q_left, pThread->q_right, pivotIndex);
			 
			pThread->left_son->q_indices = pThread->q_indices;
			pThread->left_son->q_dim = pThread->q_dim;
			pThread->left_son->q_left = pThread->q_left;
			pThread->left_son->q_right = pivotNewIndex - 1;
			pThread->right_son->q_indices = pThread->q_indices;
			pThread->right_son->q_dim = pThread->q_dim;
			pThread->right_son->q_left = pivotNewIndex + 1;
			pThread->right_son->q_right = pThread->q_right;
			 
			pthread_create(&pThread->left_son->m_thread,  NULL, pThread->left_son->QuickSort,  (void*)pThread->left_son );
			pthread_create(&pThread->right_son->m_thread, NULL, pThread->right_son->QuickSort, (void*)pThread->right_son);
			pthread_join(pThread->left_son->m_thread,  NULL);
			pthread_join(pThread->right_son->m_thread, NULL);
		}

		return 0;

	}

	static void* StartMinimize(void *ptr)
	{
		CDivConqMinThread* pThread =  (CDivConqMinThread*)ptr;
		
		int size = *(pThread->m_M_size);
		if (size<2) return 0;

		memset(pThread->f,1,size);
		pThread->f[size] = 0;

		int j;
		int *indices = pThread->NEW_MEM(size, 0);
		for (j=0; j<size; j++)
		{
			indices[j]=j;
		}

		//::QuickSort(pThread->m_M, indices, cur_n_str-1, 0, size-1);
		pThread->q_indices = indices;
		pThread->q_dim = cur_n_str-1;
		pThread->q_left = 0;
		pThread->q_right = size-1;
		CDivConqMinThread::QuickSort(pThread);

		pThread->Divide(indices, cur_n_str-1, size);

		pThread->RELEASE_MEM(size);

		int storeIdx=0;
		j=size;
		int **source = pThread->m_M;
		while (pThread->f[storeIdx] && storeIdx < j) storeIdx++;
		while(storeIdx < j) 
		{
			while (pThread->f[storeIdx] && storeIdx < j) storeIdx++;
			while (!pThread->f[j] && storeIdx < j) j--;
			if (storeIdx < j){
				memmove(source[storeIdx++], source[j--], cur_n_str*sizeof(int));
				while (pThread->f[storeIdx] && storeIdx <= j) storeIdx++;
			}
		}
		*(pThread->m_M_size) = storeIdx;
		return 0;
	}
};	

CFindParentThread* parent_finder;
CDivConqMinThread* divcon_thread;

	/*---------Get Minima. END---------*/


/*******************
* Level 1   E N D  *
*******************/
/*---------MAIN PROGRAM---------*/

bool ReSizeBuffer();

void FindMLCS(){
	
	int i, i1;  
	int size_cur_Dk = 0;
	int mlcs_size = 0;
	int flag_child = 0;
	int cur_D_member = 0;

	char* mlcs = new char[MAXSTRSIZE+1];
	
	long DominPntNum = 1;
	
	//long start_time = clock();
	time_t time_start,time_end;
	time(&time_start);

	cur_D.push_back(Allocate2DArray<int>(1,cur_n_str));
	for (i=0; i<cur_n_str; i++) 
	{
		cur_D[0][0][i] = -1;
	}/* for i */

	size_cur_D.push_back(1);
	size_cur_Dk=1;
	cur_k=0;

	//clock_t partime = 0, mintime = 0;
	// main loop
	while ((size_cur_Dk != 0))
	{
		for(i=0; i<ABC_len; i++)
			M_size[i] = 0;
	
		//clock_t start_t = clock();

		//Find parents
		if (size_cur_Dk <= MAXCPU){
			threadNum = size_cur_Dk;
			for(i=0; i<threadNum; i++) {			
				parent_finder[i].Initialize(i, i+1);
			}			
		}
		else
		{
			threadNum = MAXCPU;
			int job_per_thread = size_cur_Dk/MAXCPU;		
			int rest_job = size_cur_Dk%MAXCPU;	

			for(i=0; i<rest_job; i++) {			
				parent_finder[i].Initialize(i*(job_per_thread+1), (i+1)*(job_per_thread+1));
			}			
			int allocated = rest_job*(job_per_thread+1);
			for(i=rest_job; i<threadNum; i++) {			
				parent_finder[i].Initialize(allocated+(i-rest_job)*job_per_thread, 
										  allocated+(i+1-rest_job)*job_per_thread);
			}			
		}

		if (threadNum==1)
		{
			CFindParentThread::StartThread((void*)&parent_finder[0]);

			for(i=0; i<ABC_len; i++){	
				CDivConqMinThread::StartMinimize((void*)&divcon_thread[i]);
			}
		}
		else
		{
			for(i=0; i<threadNum; i++) {	
				pthread_create( &parent_finder[i].m_thread, NULL, CFindParentThread::StartThread, (void*)&parent_finder[i]);
			}
			for(i=0; i<threadNum; i++) {
				pthread_join(parent_finder[i].m_thread, NULL);
			}	

			for(i=0; i<ABC_len; i++) {	
				pthread_create( &divcon_thread[i].m_thread, NULL, CDivConqMinThread::StartMinimize, (void*)&divcon_thread[i]);
			}		
			for(i=0; i<ABC_len; i++) {
				pthread_join(divcon_thread[i].m_thread, NULL);
			}
		}

		//end_t = clock();
		//mintime = end_t-start_t;
		//printf("%d ,       %d\n", partime, mintime);		
		
		//printf("Finished computating minima!\n");

		size_cur_Dk = 0;
		
		for(i=0; i<ABC_len; i++)
		{
			if (M_size[i] > 0) size_cur_Dk += M_size[i];
		}

		if (size_cur_Dk>0)
		{
			cur_k+=1;
			cur_D.push_back(Allocate2DArray<int>(size_cur_Dk,cur_n_str));
			//printf("Parent set reduced to the set of minima of size %d\n",size_cur_Dk);

			//printf("%d  \n", cur_k);

			size_cur_Dk = 0;
			for(i=0; i<ABC_len; i++)
			{
				if (M_size[i] <= 0) continue;
				int **source = M[i];
				memcpy(cur_D[cur_k][size_cur_Dk], source[0], M_size[i]*cur_n_str*sizeof(int));
				size_cur_Dk += M_size[i];

			}
			size_cur_D.push_back(size_cur_Dk);				
			
			DominPntNum += size_cur_Dk;

		}
		
		if (MAXDSIZE < size_cur_Dk*ABC_len)
		{
			if (ABC_len < 20)
				MAXDSIZE *= ABC_len;
			else
				MAXDSIZE *= ABC_len/2;

			if (!ReSizeBuffer())
			{
				printf("\n\nError: out of memory!\n");
				exit(1);
			}
			//printf("k=%d\t\t D=%d\t\t M=%d\n",cur_k,size_cur_Dk,MAXDSIZE);
		}
		//printf("Dominant set of level %d: %d\n",cur_k,size_cur_Dk);
		//***System.out.println("----------------------------------------\n");
	} /* while */


	//getting the MLCS itself
	mlcs_size = cur_k;
	cur_D_member = 0;
	mlcs[mlcs_size-1] = data[0][cur_D[mlcs_size][cur_D_member][0]];
	//printf("The mlcs[%d]='%c'\n", mlcs_size, mlcs[mlcs_size-1]);
	while (mlcs_size>1)
	{
		flag_child = 0;
		for (i = 0; i<size_cur_D[mlcs_size-1]; i++)
		{
			for (i1=0; i1< cur_n_str; i1++) {
				if (cur_D[mlcs_size-1][i][i1]>= cur_D[mlcs_size][cur_D_member][i1]) {
					break;
				}
			}
			if (i1 >= cur_n_str)
			{
				mlcs_size--;
				cur_D_member = i;
				flag_child = 1;
				break;
			}
		}

		mlcs[mlcs_size-1] = data[0][cur_D[mlcs_size][cur_D_member][0]];
		//printf("The mlcs[%d]='%c'\n", mlcs_size, mlcs[mlcs_size-1]);
	}//while

	mlcs[cur_k] = 0;
	//printf("The MLCS is '%s', of length %d",+mlcs[cur_k-1],cur_k-1);	
			
	//float elapsedTime = (clock()-start_time)/CLOCK_PER_SEC;
	time(&time_end);
	double elapsedTime = difftime(time_end,time_start);

	//printf("\nmlcs is: %s\n",mlcs);
	printf("\n[MLCS](length=%d)\n%s\n",cur_k, mlcs);	

	//Generate Superposed Allignment
	char** buf =	Allocate2DArray<char>(cur_n_str+1, 2*MAXPOSSIBLESTR); 
	
	int* data_index = new int[cur_n_str];
    for (i=0; i< cur_n_str; i++) {
    	data_index[i] = 0;
    }
	
	mlcs_size = cur_k;
	int buf_index  = 0;	
	
	for (int mlcs_index = 0; mlcs_index <= mlcs_size; mlcs_index++){
   	    while(true){
    	    for (i=0; i< cur_n_str; i++) {
    	    	if (data_index[i] < str_lngth[i]){
    	    		break;
    	    	}
    	    }

    	    // have finished processing all strings
    	    if (i == cur_n_str)
    	    	break;
    	    
    	    for (i=0; i< cur_n_str; i++) {
    	    	if (data_index[i] >= str_lngth[i] || 
    	    		data[i][data_index[i]] != mlcs[mlcs_index]) break;
    	    }
    	    
    	    if (i< cur_n_str){	//No alignment.
        	    for (i=0; i < cur_n_str; i++) {
        	    	if (data_index[i] < str_lngth[i] && 
        	    			data[i][data_index[i]] != mlcs[mlcs_index]) {
        	    		buf[i][buf_index] = data[i][data_index[i]];
        	    		data_index[i]++;
        	    	}else{
        	    		buf[i][buf_index] = '-';       	    		
        	    	}
        	    }
				buf[i][buf_index] = ' '; 
	    	    buf_index++;    	    	
     	    }else{				//Find an alignment.
           	    for (i=0; i < cur_n_str; i++) {
         	    	buf[i][buf_index] = data[i][data_index[i]];
        	    	data_index[i]++;
        	    }
				buf[i][buf_index] = '*'; 
           	    buf_index++;
           	    break;
    	    }
	    }
	}
	
    for (i=0; i < cur_n_str+1; i++) {
    	buf[i][buf_index] = 0;
    }
	
	printf("\n[Alignment]\n");
	printf("%s\n",buf[cur_n_str]);	
    for (i=0; i< cur_n_str; i++) {
		printf("%s\n",buf[i]);	
    }
	printf("\n[Number of dominant points]\n");


	printf("Total  : %d\n",DominPntNum);		
	printf("Average: %d/%d=%d\n",DominPntNum,cur_k,(int)(DominPntNum/cur_k));		

	printf("\n[Running time]\n%6.2f seconds\n",elapsedTime);			
	Free2DArray<char>(buf);
	delete[] data_index;

	delete[] mlcs;
}

bool ReSizeBuffer()
{
	int i, **PPTR;

	if (f)
	{
		Free2DArray<bool>(f);

		PPTR = M[0];
		Free2DArray<int>(PPTR);
		
		PPTR = buffer[0];
		Free2DArray<int>(PPTR);
	}

	f    = Allocate2DArray<bool>(ABC_len, MAXDSIZE);
	if (!f) goto MEM_FAIL;

	PPTR = Allocate2DArray<int>(MAXDSIZE*ABC_len, MAXSTR);
	if (!PPTR) goto MEM_FAIL;
	
	for (i=0; i<ABC_len; i++)
	{
		M[i] = PPTR;
		PPTR += MAXDSIZE;
	}

	if (MAXCPU<=ABC_len)
	{
		int buf_num = ABC_len*MAXCPU;
		if (MAXCPU == 1) buf_num *= 2;

		//buffer = Allocate3DArray<int>(buf_num, MAXDSIZE/ABC_len, MAXSTR);
		buffer.resize(buf_num);

		int   nCol = MAXDSIZE/ABC_len;
		int** PPTR = Allocate2DArray<int>(nCol*buf_num, MAXSTR);
		if (!PPTR) goto MEM_FAIL;

		for (i=0; i<buf_num; i++)
		{
			buffer[i] = PPTR;
			PPTR += nCol;
		}

		int*** temp = &buffer[0];
		for (i=0; i<MAXCPU; i++)
		{
			parent_finder[i].parent = temp;
			temp += ABC_len;
		}

		for (i=0; i<ABC_len; i++)
		{
			divcon_thread[i].InitThread(f[i], M[i], &M_size[i], *buffer[2*i], *buffer[2*i+1]);
		}
	}
	else
	{
		int buf_num = ABC_len*MAXCPU;
		if (buf_num < 120) buf_num = 120;

		//buffer = Allocate3DArray<int>(buf_num, MAXDSIZE/ABC_len, MAXSTR);
		buffer.resize(buf_num);

		int   nCol = MAXDSIZE/ABC_len;
		int** PPTR = Allocate2DArray<int>(nCol*buf_num, MAXSTR);
		if (!PPTR) goto MEM_FAIL;

		for (i=0; i<buf_num; i++)
		{
			buffer[i] = PPTR;
			PPTR += nCol;
		}

		int*** temp = &buffer[0];
		for (i=0; i<MAXCPU; i++)
		{
			parent_finder[i].parent = temp;
			temp += ABC_len;
		}

		for (i=0; i<ABC_len; i++)
		{
			divcon_thread[i].InitThread(f[i], M[i], &M_size[i], *buffer[2*i], *buffer[2*i+1]);
		}

		for (i=0; i<ABC_len; i++)
		{
			divcon_thread[i].left_son  = &divcon_thread[ABC_len+3*i];
			divcon_thread[i].right_son = &divcon_thread[ABC_len+3*i+1];	
			divcon_thread[i].mid_son   = &divcon_thread[ABC_len+3*i+2];	

			divcon_thread[i].left_son->InitThread(  &divcon_thread[i], 
													*buffer[2*ABC_len+6*i], 
													*buffer[2*ABC_len+6*i+1]);

			divcon_thread[i].right_son->InitThread( &divcon_thread[i], 
													*buffer[2*ABC_len+6*i+2], 
													*buffer[2*ABC_len+6*i+3]);

			divcon_thread[i].mid_son->InitThread(   &divcon_thread[i], 
													*buffer[2*ABC_len+6*i+4], 
													*buffer[2*ABC_len+6*i+5]);
		}

		if (MAXCPU>2*ABC_len && MAXCPU < 32)
		{
			for (i=0; i<3*ABC_len; i++)
			{
				divcon_thread[ABC_len+i].left_son  = &divcon_thread[4*ABC_len+3*i];
				divcon_thread[ABC_len+i].right_son = &divcon_thread[4*ABC_len+3*i+1];
				divcon_thread[ABC_len+i].right_son = &divcon_thread[4*ABC_len+3*i+2];

				divcon_thread[ABC_len+i].left_son->InitThread(  &divcon_thread[ABC_len+i], 
																*buffer[8*ABC_len+6*i], 
																*buffer[8*ABC_len+6*i+1]);

				divcon_thread[ABC_len+i].right_son->InitThread( &divcon_thread[ABC_len+i], 
																*buffer[8*ABC_len+6*i+2], 
																*buffer[8*ABC_len+6*i+3]);

				divcon_thread[ABC_len+i].mid_son->InitThread(   &divcon_thread[ABC_len+i], 
																*buffer[8*ABC_len+6*i+4], 
																*buffer[8*ABC_len+6*i+5]);
			}
		}
	}

	return true;

MEM_FAIL:
	return false;
}

int main(int argc, char *argv[]){	
	//***System.out.println("MLCS Version 1.0\n");
	//get input parameters and initialize global variables
	int i;
	if(argc<5){
		printf("Too few parameters! Please provide parameters in order as follows,\n");
		printf("(1) The number of threads.\n");
		printf("(2) The number of Strings.\n");
		printf("(3) Alphabet size (4 or 20).\n");
		printf("(4) Input file name.\n");
		exit(1);
	}

	MAXCPU = atoi(argv[1]);
	MAXSTR = atoi(argv[2]);

	ABC_len = atoi(argv[3]);
	if (ABC_len==4){
		ABC_str = ABC_str_4;
		ABC_idx = ABC_idx_4;
	}else if (ABC_len==20){
		ABC_str = ABC_str_20;
		ABC_idx = ABC_idx_20;
	}else{
		printf("Unsupported alphabet size!\n");
		exit(1);
	}

	strcpy(data_file_fullname, argv[4]);

	//check the content of the input sequence file
	FILE *my_file = fopen(data_file_fullname, "r");
	if (my_file == NULL)
	{
		printf("Empty string file: %s!\n", data_file_fullname);
		exit(1);
	}
	else
	{
		MAXSTRSIZE = 0;
		char cur_line[MAXPOSSIBLESTR];

		while (!feof(my_file) && (MAXSTR > cur_n_str))  
		{
			fscanf(my_file, "%s\n", cur_line);
			//sscanf(cur_line, "%s", cur_line);
			int len = strlen(cur_line);
			if (!len || cur_line[0]=='>') continue;
			if (len > MAXSTRSIZE)
				 MAXSTRSIZE = len;
			cur_n_str ++;
		}
		fclose(my_file);
	}
	
	if (cur_n_str == 0){
		printf("\n%d input sequence:\n", cur_n_str);
		exit(1);
	}else if (MAXSTR > cur_n_str){
		MAXSTR = cur_n_str;
	}

	//if (MAXSTR>1)
	//	printf("\n%d input sequences\n", MAXSTR);
	//else
	//	printf("\n%d input sequence\n", MAXSTR);

	//allocate memory for global variables
	data = Allocate2DArray<char>(MAXSTR, MAXSTRSIZE+1);
	str_lngth = new int[MAXSTR];

	cur_T.resize(MAXSTR);
	int **PPTR = Allocate2DArray<int>((MAXSTRSIZE+1)*MAXSTR, ABC_len);

	for (i = 0; i < MAXSTR; i++)
	{
		cur_T[i] = PPTR;
		PPTR += MAXSTRSIZE+1;
	}

	//read input strings and do preprocessing
	cur_n_str = 0;
	my_file = fopen(data_file_fullname, "r");
	while (!feof(my_file) && (cur_n_str < MAXSTR))  
	{
		fscanf(my_file, "%s\n", data[cur_n_str]);
		int len = strlen(data[cur_n_str]);
		if (!len || data[cur_n_str][0]=='>')continue;
		
		str_lngth[cur_n_str]=len;
		//printf("sequence %d: '%s'\n", cur_n_str+1, data[cur_n_str]);    
				
		for (i=0; i<ABC_len; i++)
			cur_T[cur_n_str][len][i] = len+1;

		for (int j=len-1; j>=0; j--)
		{
			for (i=0; i<ABC_len; i++)
				cur_T[cur_n_str][j][i] = cur_T[cur_n_str][j+1][i];

			if (data[cur_n_str][j]>='a')
				data[cur_n_str][j]=data[cur_n_str][j]-'a'+'A';

			cur_T[cur_n_str][j][ABC_idx[data[cur_n_str][j]-'A']] = j;
		}

		cur_n_str ++;
	}
	fclose(my_file);

	M_size = new int[ABC_len];
	M.resize(ABC_len);

	parent_finder = new CFindParentThread[MAXCPU];
	divcon_thread = new CDivConqMinThread[(MAXCPU<=ABC_len? ABC_len : 80)];

	ReSizeBuffer();

	//find MLCS of strings
	FindMLCS();

	delete[] M_size;
	delete[] str_lngth;
	delete[] parent_finder;
	delete[] divcon_thread;

	//release memory
	Free2DArray<char>(data);
	Free2DArray<bool>(f);

	PPTR = M[0];
	Free2DArray<int>(PPTR);
	
	PPTR = buffer[0];
	Free2DArray<int>(PPTR);

	for (i=0; i<cur_D.size(); i++){
		PPTR = cur_D[i];
		Free2DArray<int>(PPTR);
	}

	PPTR = cur_T[0];
	Free2DArray<int>(PPTR);

	return 0;
}
