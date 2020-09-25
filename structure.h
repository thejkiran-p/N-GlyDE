#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cstring>
#include <iterator>
#include <algorithm>
#include <pthread.h>

#define Level2Base 0x40
#define Level1Base 0x4000000 // 4^13

#define QueryFilePath  "/var/www/html/GlycoPred/queries/"
#define ResultFilePath "/var/www/html/GlycoPred/results/"

using namespace std;

typedef struct
{
	string seq;
	string ac_number;
} ProteinItem_t;

typedef struct
{
	int tid;
	string InputFilename;
} ThreadParameter_t;

// Global Variables
extern int PositiveNum;
extern pthread_mutex_t Lock;
extern vector<string> ProteinList;
extern vector<bool> ThreadStatusVec;
extern map<string, vector<int> > HomologyMap; //first:database protein, second: training protein list (homology)
extern vector<ThreadParameter_t> ThreadParameterVec;

// GetData.cpp
extern void BuildHomologyMap();
extern void GetDBProteinList(string filename);
extern vector<string> GetMyHomologs(string folder, string protein_name);

// Tools.cpp
extern bool CheckFileExistence(const char* filename);
extern float PairwiseAlignment(string seq1, string seq2);
extern string ChangeFileSubName(string QueryFileName, string SubName);

// Prediction.cpp
extern void *HHBlitSearch(void *Arg);
extern void *SecondStagePred(void *Arg);
extern void FirstStagePrediction(string QueryFileName);
