// heuristic.cpp : Defines the entry point for the console application.
//
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "tree.h"


/*-------- Global Variables. Begin ------------*/

char	ABC_str_20[] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
char	ABC_str_4[]  = {'A','C','G','T'};
char	ABC_idx_20[] = { 0,0,1,2,3,4,5,6,7,7,8,9,10,11,11,12,13,14,15,16,16,17,18,18,19,19};
char	ABC_idx_4[]  = { 0,0,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3};

int	MAXSTRSIZE;
int	MAXSTR;
int	ABC_len;
char*	ABC_str;
char*	ABC_idx;

char**	data;
int*	str_lngth;
int	cur_n_str = 0;
char	data_file_fullname[128];

std::vector <short**>	cur_T;
std::vector <short**>	cur_M;

/*-------- AVL Tree.  Begin ------------*/

int NUMDIM  = 1;
int NUMBYTE = 1;

typedef struct _Node {
	union {
		unsigned char cXYZ[16]; 
		unsigned short sXYZ[8];

		#ifdef __int64
			__int64 lXYZ[2];
		#else
			unsigned long lXYZ[2];
		#endif
	} uXYZ;

	short			h;
	int			fg;
	bool			expended;
	_Node			*parent;
	TREE_ENTRY(_Node)	linkage;
} Node;

typedef TREE_HEAD(_Tree, _Node) Tree;

TREE_DEFINE(_Node, linkage);


_Node *Node_new(_Node *rhs) {
	_Node *self = (_Node *)malloc(sizeof(_Node));
	if (rhs) {
		memcpy(self, rhs, sizeof(_Node));
	}
	return self;
}

void Node_free(_Node *rhs) {
	if (rhs) free(rhs);
}

int Node_compare(_Node *lhs, _Node *rhs) {
	if (NUMDIM == 1) {
		return lhs->uXYZ.lXYZ[0] > rhs->uXYZ.lXYZ[0] ? 1 : (lhs->uXYZ.lXYZ[0] < rhs->uXYZ.lXYZ[0] ? -1 : 0);
	} else {
		if (lhs->uXYZ.lXYZ[0] == rhs->uXYZ.lXYZ[0]) {
			return lhs->uXYZ.lXYZ[1] > rhs->uXYZ.lXYZ[1] ? 1 : (lhs->uXYZ.lXYZ[1] < rhs->uXYZ.lXYZ[1] ? -1 : 0);
		} else {
			return lhs->uXYZ.lXYZ[0] > rhs->uXYZ.lXYZ[0] ? 1 : -1;
		}
	}

	/*for(int i=0; i<MAXSTR; i++){
		int result = lhs->xyz[i]-rhs->xyz[i];

		if (result)
			return result;
	}*/
	return 0;
}

/*-------- Heap . Begin ------------*/

std::vector <_Node*> myHeap(1000, 0);
int heap_sz = 0;

void SiftDown() {
	int parent = 0;
	while(parent<heap_sz) {
		int l = 2 * parent + 1;
		int r = l + 1;
		int flag = parent;
	 
		if (l<heap_sz && (myHeap[l]->fg > myHeap[flag]->fg)) {
			flag = l;
		}
		if (r<heap_sz && (myHeap[r]->fg > myHeap[flag]->fg)) {
			flag = r;
		}
		if (flag == parent) {
			break;
		} else {
			Node* temp = myHeap[flag];
			myHeap[flag] = myHeap[parent];
			myHeap[parent] = temp;
			parent = flag;
		}
	}
}

void SiftUp() {
	int son = heap_sz-1;
	while(son > 0) {
		int parent = son / 2;
	 
		if (myHeap[parent]->fg > myHeap[son]->fg) {
			break;
		} else {
			Node* temp = myHeap[son];
			myHeap[son] = myHeap[parent];
			myHeap[parent] = temp;
			son = parent;
		}
	}
}

void HeapInsert(_Node* p) {
	if (heap_sz == myHeap.size())
		myHeap.push_back(p);
	else
		myHeap[heap_sz] = p;
	heap_sz++;

	SiftUp();
}

_Node* HeapRemove() {
	_Node* ret = myHeap[0];
	myHeap[0] = myHeap[heap_sz-1];
	heap_sz--;

	SiftDown();
	return ret;
}

/*-------- functions. Begin ------------*/

void ComputeM() {
	int i;
	for (i=0; i<(MAXSTR+1)/2; i++) {
		int i1 = i << 1;
		int i2 = 1+ (i << 1);
		if (i2>=MAXSTR) i2=0;

		short** Mi = Allocate2DArray<short>(str_lngth[i1]+1, str_lngth[i2]+1);
		cur_M.push_back(Mi);

		int x, y;
		for (x=0; x<=str_lngth[i1]; x++)
			Mi[x][str_lngth[i2]] = 0;

		for (x=0; x<=str_lngth[i2]; x++)
			Mi[str_lngth[i1]][x] = 0;

		for (x=str_lngth[i1]-1; x>=0; x--) {
			for (y=str_lngth[i2]-1; y>=0; y--) {
				if (data[i1][x] == data[i2][y])
					Mi[x][y] = Mi[x+1][y+1]+1;
				else
					Mi[x][y] = max(Mi[x+1][y], Mi[x][y+1]);
			}
		}
	}
}

void ComputeF(_Node* p) {
	int h = cur_M[0][0][0];
	for (int i=0; i<(MAXSTR+1)/2; i++) {
		int i1 = i << 1;
		int i2 = 1+ (i << 1);
		if (i2>=MAXSTR) i2=0;
		if (NUMBYTE == 1)
			h = min(h, cur_M[i][p->uXYZ.cXYZ[i1]+1][p->uXYZ.cXYZ[i2]+1]);
		else
			h = min(h, cur_M[i][p->uXYZ.sXYZ[i1]+1][p->uXYZ.sXYZ[i2]+1]);
	}

	p->h = h;
	h = p->fg & 0xffff;
	p->fg = ((h + p->h) << 16) | h;
}

int CommonSeq(_Node* p) {
	if (!p) {
		return 0;
	}
	int mlcs_len = CommonSeq(p->parent);
		
	if (NUMBYTE == 1)
		printf("%c",data[0][p->uXYZ.cXYZ[0]]);
	else
		printf("%c",data[0][p->uXYZ.sXYZ[0]]);
	
	return mlcs_len+1;
}

/*-------- Main loop. Begin ------------*/

void Mainloop() {	
	long closedset_sz = 0;
	long openset_sz   = 0;
	
	//long start_time = clock();
	time_t time_start,time_end;
	time(&time_start);

	Tree tree= TREE_INITIALIZER(Node_compare);
	
	_Node p0 = {0};
	p0.uXYZ.lXYZ[0] = 0;
	p0.uXYZ.lXYZ[1] = 0;
	ComputeF(&p0);

	_Node* temp = Node_new(&p0);
	temp->uXYZ.lXYZ[0] = 0;
	temp->uXYZ.lXYZ[1] = 0;

	_Node* p = Node_new(temp);
	TREE_INSERT(&tree, _Node, linkage, p);
	HeapInsert(p);
	openset_sz++;

	while(true) {
		_Node* p = HeapRemove();
		p->expended = true;
		closedset_sz++;
		openset_sz--;

		if (p->h == 0) 	{
			printf("\nThe MLCS is:\n");
			int l = CommonSeq(p);
			printf("\n\n|MLCS|=%d\n", l);
			break;
		}

		for(int i=0; i<ABC_len; i++) { 
			 int j=0;

			 if (NUMBYTE == 1) {
				 while (j<MAXSTR) {
					 temp->uXYZ.cXYZ[j] = cur_T[j][p->uXYZ.cXYZ[j]+1][i];

					 if  (temp->uXYZ.cXYZ[j]==str_lngth[j]+1)
						 break;
					 j++; 
				 }
			 } else {
				 while (j<MAXSTR) {
					 temp->uXYZ.sXYZ[j] = cur_T[j][p->uXYZ.sXYZ[j]+1][i];

					 if  (temp->uXYZ.sXYZ[j]==str_lngth[j]+1)
						 break;
					 j++; 
				 }
			 }

			 if (j < MAXSTR)
				 continue;

			_Node *q = TREE_FIND(&tree, _Node, linkage, temp);

			if (q) { //node is already in tree
				if (!q->expended && (q->fg & 0xffff) < (p->fg & 0xffff)+1) {
					q->parent = p;
					q->fg = ((q->h+(p->fg & 0xffff)+1) << 16) | ((p->fg & 0xffff)+1);
				}
			} else {
				q = Node_new(temp);
				q->parent = p;
				q->fg = p->fg+1;
				ComputeF(q);

				TREE_INSERT(&tree, _Node, linkage, q);
				HeapInsert(q);
				openset_sz++;
			}
		}

	}

	printf("|closedset|=%d\n|openset|=%d\n\n", closedset_sz, openset_sz+1);

	time(&time_end);
	double elapsedTime = difftime(time_end,time_start);

	printf("The time is %6.2f seconds\n",elapsedTime);			
}

/*-------- Memory Management. Begin ------------*/

int main(int argc, char* argv[]) {
	//***System.out.println("MLCS Version 1.0\n");
	//get input parameters and initialize global variables
	if(argc<3){
		printf("Too few parameters! Please provide parameters in order as follows,\n");
		printf("(1) The number of Strings.\n");
		printf("(2) Alphabet size (4 or 20).\n");
		printf("(3) Input file name.\n");
		exit(1);
	}

	MAXSTR = atoi(argv[1]);

	ABC_len = atoi(argv[2]);
	if (ABC_len==4){
		ABC_str = ABC_str_4;
		ABC_idx = ABC_idx_4;
	}else if (ABC_len==20){
		ABC_str = ABC_str_20;
		ABC_idx = ABC_idx_20;
	}else{
		printf("Unsupported alphabet size!\n");
		return 1;
	}

	strcpy(data_file_fullname, argv[3]);

	//check the content of the input sequence file
	FILE *my_file = fopen(data_file_fullname, "r");
	if (my_file == NULL) {
		printf("Empty string file: %s!\n", data_file_fullname);
		return 1;
	} else {
		MAXSTRSIZE = 0;
		char cur_line[10000];

		while (!feof(my_file)) {
			fscanf(my_file, "%s\n", cur_line);
			//sscanf(cur_line, "%s", cur_line);
			int len = strlen(cur_line);
			if (!len) continue;
			if (len > MAXSTRSIZE)
				 MAXSTRSIZE = len;
			cur_n_str ++;
		}
		fclose(my_file);
	}
	
	if (cur_n_str == 0){
		printf("\nNumber of input strings is %d.\n", cur_n_str);
		return 1;
	}else if (MAXSTR > cur_n_str){
		MAXSTR = cur_n_str;
	}
	printf("\nNumber of input strings is %d.\n", MAXSTR);

	//allocate memory for global variables
	data = Allocate2DArray<char>(MAXSTR, MAXSTRSIZE+1);
	str_lngth = new int[MAXSTR];

	//read input strings and do preprocessing
	cur_n_str = 0;
	int i, j;
	my_file = fopen(data_file_fullname, "r");
	while (!feof(my_file) && (cur_n_str < MAXSTR)) {
		fscanf(my_file, "%s\n", data[cur_n_str]);
		int len = strlen(data[cur_n_str]);
		if (!len)continue;
		
		str_lngth[cur_n_str]=len;
		//printf("string N%d: '%s'\n", cur_n_str, data[cur_n_str]);    
				
		cur_T.push_back(Allocate2DArray<short>(len+1, ABC_len));
		
		for (i=0; i<ABC_len; i++) {
			cur_T[cur_n_str][len][i] = len+1;
		}

		for (j=len-1; j>=0; j--) {
			for (i=0; i<ABC_len; i++) {
				cur_T[cur_n_str][j][i] = cur_T[cur_n_str][j+1][i];
			}
			cur_T[cur_n_str][j][ABC_idx[data[cur_n_str][j]-'A']] = j;
		}

		cur_n_str ++;
	}
	fclose(my_file);

	ComputeM();

	//find MLCS of strings
	Mainloop();

	//release memory
	Free2DArray<char>(data);

	delete[] str_lngth;
	
	for (i=0; i<cur_T.size(); i++){
		short** temp = cur_T[i];
		Free2DArray<short>(temp);
	}

	for (i=0; i<(MAXSTR+1)/2; i++) {
		short** temp = cur_M[i];
		Free2DArray<short>(temp);
	}

	printf("Done!\n");
	return 0;
}
