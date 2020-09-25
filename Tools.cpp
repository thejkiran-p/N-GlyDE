#include "structure.h"

bool CheckFileExistence(const char* filename)
{
	bool bExist = true;

	FILE *file = fopen(filename, "r");
	if (file) fclose(file);
	else bExist = false;

	return bExist;
}

string ChangeFileSubName(string QueryFileName, string SubName)
{
	string cmd, NewFileName;

	NewFileName = QueryFileName.substr(0, (int)QueryFileName.find_last_of('.'));
	NewFileName = NewFileName + SubName;

	cmd = "mv " + QueryFileName + " " + NewFileName;
	system(cmd.c_str());
	
	return NewFileName;
}

float PairwiseAlignment(string str1, string str2)
{
	string str;
	fstream file;
	int i, identity;

	file.open("pairwise_alignment.seq", ios_base::out);
	file << ">1" << endl << str1 << endl << ">2" << endl << str2 << endl;
	file.close();

	system("/home/arith/clustalw kbloc_alignment.seq > /dev/null");

	file.clear(); file.open("kbloc_alignment.aln", ios_base::in);
	getline(file, str); getline(file, str); getline(file, str);

	for(identity=0, i=0;i<(int)str.length();i++) if(str[i] == '*') identity++;

	return 100.0*identity/(int)str.length();
}
