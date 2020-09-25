#include "structure.h"
#include <dirent.h>
#include <sys/types.h>
#include <unistd.h>

#define NUM_THREADS 10

int PositiveNum;
pthread_mutex_t Lock;
vector<string> ProteinList;
vector<bool> ThreadStatusVec;
map<string, vector<int> > HomologyMap;
vector<ThreadParameter_t> ThreadParameterVec;

int main(int argc, char* argv[])
{
	int tid;
	pthread_t ThreadArr[NUM_THREADS];
	
	//DIR *dirp;
	stringstream ss;
	//struct dirent *dp;
	fstream file, f, msgFile;
	map<int, int> MyTemplateMap;
	vector<double> PredictionResult;
	string str, cmd, QueryFileName, QueryType, JobID, EmailAddr;

	GetDBProteinList("/home/arith/Protemoics/Glycosylation/Positives.fa");
	PositiveNum = (int)ProteinList.size();
	GetDBProteinList("/home/arith/Protemoics/Glycosylation/Negatives.fa");
	BuildHomologyMap();

	ThreadStatusVec.clear(); ThreadStatusVec.resize(NUM_THREADS);
	ThreadParameterVec.clear(); ThreadParameterVec.resize(NUM_THREADS);

	while(true)
	{
		// job
		system("/bin/ls -tr /var/www/html/GlycoPred/queries/*.job > job_list.txt");
		file.clear(); file.open("job_list.txt", ios_base::in);
		while (!file.eof())
		{
			getline(file, QueryFileName); if (QueryFileName == "") break;
			if (QueryFileName == "/var/www/html/GlycoPred/queries/null.job") continue;
			//cout << QueryFileName << endl;
			bool bDone = false;
			int iFinish = 0, iTotal = 0;
			f.clear(); f.open(QueryFileName.c_str()); getline(f, JobID); getline(f, EmailAddr);

			if (EmailAddr != "none" && EmailAddr.find('@')>0)
			{
				string ResFileName;
				while (!f.eof())
				{
					getline(f, str); if (str == "") break;
					ResFileName = "/var/www/html/GlycoPred/results/" + str.substr(0, 17) + ".res";
					iTotal++;
					if (CheckFileExistence(ResFileName.c_str())) iFinish++;
				}
				if (iTotal>0 && iFinish == iTotal)
				{
					bDone = true;
					msgFile.clear(); msgFile.open("msg.txt", ios_base::out);
					msgFile << "All your sequences have been done prediction." << endl << "Please go to the URL: http://bioapp.iis.sinica.edu.tw/GlycoPred/MakeSummary.php?job=" << JobID << " to check out the prediction results.\n" << endl;
					msgFile.close();
					cmd = "/usr/bin/mutt -s Glycosylation_Prediction_Result " + EmailAddr + " < msg.txt";
					system(cmd.c_str());
				}
			}
			else bDone = true;

			f.close();
			if (bDone)
			{
				cmd = "rm " + QueryFileName;
				system(cmd.c_str());
			}
		}
		file.close();
		system("/bin/sleep 1s");

		// query
		system("/bin/ls -tr /var/www/html/GlycoPred/queries/*.s2 > job_list.txt");
		file.clear(); file.open("job_list.txt", ios_base::in);
		while (!file.eof())
		{
			getline(file, QueryFileName); if (QueryFileName == "") break;
			if (QueryFileName == "/var/www/html/GlycoPred/queries/null.s2") continue;
			//cerr << QueryFileName << endl;
			for (tid = 0; tid < NUM_THREADS; tid++)
			{
				if (ThreadStatusVec[tid] == false)
				{
					pthread_mutex_lock(&Lock);
					ThreadStatusVec[tid] = true;
					pthread_mutex_unlock(&Lock);
					ThreadParameterVec[tid].tid = tid;
					QueryFileName = ChangeFileSubName(QueryFileName, ".query");
					ThreadParameterVec[tid].InputFilename = QueryFileName;
					//cerr << "Send " << ThreadParameterVec[tid].InputFilename << " to thread " << tid << endl;
					pthread_create(&ThreadArr[tid], NULL, SecondStagePred, &ThreadParameterVec[tid]);
					break;
				}
			}
			if (tid == NUM_THREADS) break;
		}
		file.close();

		// homo
		system("/bin/ls -tr /var/www/html/GlycoPred/queries/*.s1 > job_list.txt");
		file.clear(); file.open("job_list.txt", ios_base::in);
		while (!file.eof())
		{
			getline(file, QueryFileName); if (QueryFileName == "") break;
			if (QueryFileName != "/var/www/html/GlycoPred/queries/null.s1")
			{
				FirstStagePrediction(QueryFileName);
				ChangeFileSubName(QueryFileName, ".s2");
			}
		}
		file.close();
		system("/bin/sleep 1s");

		// submit
		system("/bin/ls -tr /var/www/html/GlycoPred/queries/*.submit > job_list.txt");
		file.clear(); file.open("job_list.txt", ios_base::in);
		while (!file.eof())
		{
			getline(file, QueryFileName); if (QueryFileName == "") break;
			if (QueryFileName == "/var/www/html/GlycoPred/queries/null.submit") continue;
			//cerr << QueryFileName << endl;
			for (tid = 0; tid < NUM_THREADS; tid++)
			{
				if (ThreadStatusVec[tid] == false)
				{
					//system("sleep 1s");
					pthread_mutex_lock(&Lock);
					ThreadStatusVec[tid] = true;
					pthread_mutex_unlock(&Lock);
					ThreadParameterVec[tid].tid = tid;
					QueryFileName = ChangeFileSubName(QueryFileName, ".tmp_1");
					ThreadParameterVec[tid].InputFilename = QueryFileName;
					//cerr << "Send " << ThreadParameterVec[tid].InputFilename << " to thread " << tid << endl;
					pthread_create(&ThreadArr[tid], NULL, HHBlitSearch, &ThreadParameterVec[tid]);
					break;
				}
			}
			if (tid == NUM_THREADS) break;
		}
		file.close();
		system("/bin/sleep 5s");
	}

	return 0;
}
