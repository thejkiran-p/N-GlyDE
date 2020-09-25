#include "structure.h"

void GetDBProteinList(string filename)
{
	string str;
	fstream file;

	file.open(filename.c_str(), ios_base::in);
	while (!file.eof())
	{
		getline(file, str); if (str == "") break;
		if(str[0] == '>') ProteinList.push_back(str.substr(1));
	}
	file.close();
}

vector<string> GetMyHomologs(string folder, string protein_name)
{
	int p;
	float prob;
	fstream file;
	vector<string> HomologVec;
	string str, tmp, fname, homolog_name;

	fname = folder + protein_name;

	file.open(fname.c_str(), ios_base::in);
	do
	{
		getline(file, str);
	} while (str != "");

	getline(file, str); getline(file, str);

	while (str != "")
	{
		tmp = str.substr(35, 5); prob = atof(tmp.c_str());
		if (prob > 50.0)
		{
			p = str.find_first_of('|', 7); homolog_name = str.substr(7, p - 7);
			HomologVec.push_back(homolog_name);
		}
		getline(file, str);
	}
	file.close();
	sort(HomologVec.begin(), HomologVec.end());
	HomologVec.erase(unique(HomologVec.begin(), HomologVec.end()), HomologVec.end());

	return HomologVec;
}

void BuildHomologyMap()
{
	string name;
	int protein_id, num;
	vector<string> HomologVec;
	vector<string>::iterator iter;

	num = (int)ProteinList.size();
	for (protein_id = 0; protein_id < num; protein_id++)
	{
		HomologVec = GetMyHomologs("/home/arith/Protemoics/Glycosylation/HHBlits/", ProteinList[protein_id]);
		//fprintf(stderr, "\r%d: %s --> %d\t\t\t", protein_id, ProteinList[protein_id].c_str(), (int)HomologVec.size());
		for (iter = HomologVec.begin(); iter != HomologVec.end(); iter++) HomologyMap[*iter].push_back(protein_id);
	}
}
