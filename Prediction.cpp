#include <cmath>
#include "structure.h"

void *HHBlitSearch(void *Arg)
{
	int tid;
	ThreadParameter_t *MyPara;

	fstream file;
	string cmd, str, ProteinSeq, InputFileName, QueryBaseName, NSP_FileName, ResultFileName;

	MyPara = (ThreadParameter_t*)Arg; tid = (*MyPara).tid;
	InputFileName = (string)(*MyPara).InputFilename;

	QueryBaseName = InputFileName.substr(0, (int)InputFileName.find_last_of('.'));
	ResultFileName = QueryBaseName + ".hhb"; NSP_FileName = QueryBaseName + ".nsp";

	file.open(InputFileName.c_str(), ios_base::in);
	getline(file, str); getline(file, ProteinSeq);
	file.close();

	cmd = "/usr/local/hhsuite/bin/hhblits -cpu 4 -i " + InputFileName + " -d /home/arith/tools/hhsuite-3.0.1-Source/database/uniprot20_2016_02 -n 1 -hide_cons -p 50 -aliw 5000 -v 0 -o " + ResultFileName;
	system(cmd.c_str());

	cmd = "/home/arith/Protemoics/Glycosylation/netsurfp-1.0/netsurfp -i " + InputFileName + " -t FASTA -a -o " + NSP_FileName;
	system(cmd.c_str());

	ChangeFileSubName(InputFileName, ".s1");
	pthread_mutex_lock(&Lock);
	ThreadStatusVec[tid] = false;
	pthread_mutex_unlock(&Lock);
	return (void*)(1);
}

void FirstStagePrediction(string QueryFileName)
{
	fstream file;
	vector<int> vec;
	string filename;
	int i, tid, score;
	vector<string> HomologVec;
	vector<string>::iterator iter;
	double PositiveScore, NegativeScore, TotalScore;

	filename = QueryFileName.substr(0, QueryFileName.find_last_of('.')) + ".hhb";
	HomologVec = GetMyHomologs("", filename);

	for (iter = HomologVec.begin(); iter != HomologVec.end(); iter++) if (HomologyMap[*iter].size() > 1) copy(HomologyMap[*iter].begin(), HomologyMap[*iter].end(), back_inserter(vec));
	sort(vec.begin(), vec.end()); vec.push_back(0); tid = 0; score = 0;

	PositiveScore = NegativeScore = 0;
	for (i = 0; i < (int)vec.size(); i++)
	{
		if (vec[i] == tid) score++;
		else
		{
			if (score > 5)
			{
				if (tid < PositiveNum) PositiveScore += pow(score, 2.0);
				else NegativeScore += pow(score, 2.0);
			}
			tid = vec[i]; score = 1;
		}
	}
	if ((TotalScore = PositiveScore + NegativeScore) > 0)
	{
		PositiveScore = (PositiveScore / TotalScore);
		NegativeScore = (NegativeScore / TotalScore);
		//printf("Positive=%.2f, Negative=%.2f, Total=%.2f\n", PositiveScore, NegativeScore, TotalScore);
	}
	filename = QueryFileName.substr(0, QueryFileName.find_last_of('.')) + ".1st";
	file.clear(); file.open(filename.c_str(), ios_base::out);
	file << PositiveScore << " " << NegativeScore << endl;
	file.close();
}

void *SecondStagePred(void *Arg)
{
	int tid;
	ThreadParameter_t *MyPara;

	fstream file;
	string cmd, str, ProteinSeq, InputFileName, QueryBaseName, ProfileName, SecondStageFileName, ResultFileName;

	MyPara = (ThreadParameter_t*)Arg; tid = (*MyPara).tid;
	InputFileName = (string)(*MyPara).InputFilename;
	QueryBaseName = InputFileName.substr(0, (int)InputFileName.find_last_of('.'));
	
	ProfileName = QueryBaseName + ".prf"; SecondStageFileName = QueryBaseName + ".2nd"; ResultFileName = QueryBaseName + ".res";

	//generate profile
	cmd = "/usr/bin/python GenQueryProfile.py " + QueryBaseName;
	system(cmd.c_str());

	//run svm-predict
	cmd = "./svm-predict -q -b 1 " + ProfileName + " prediction.model " + SecondStageFileName;
	system(cmd.c_str());

	//generate final prediction
	cmd = "/usr/bin/python GenPredResult.py " + QueryBaseName;
	system(cmd.c_str());

	cmd = "chmod 644 " + ResultFileName;
	system(cmd.c_str());

	cmd = "mv " + ResultFileName + " " + ResultFilePath;
	system(cmd.c_str());

	cmd = "rm " + QueryBaseName + "*";
	system(cmd.c_str());

	pthread_mutex_lock(&Lock);
	ThreadStatusVec[tid] = false;
	pthread_mutex_unlock(&Lock);
	return (void*)(1);
}
