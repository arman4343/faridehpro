#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <unordered_map>
#include "My_Functions.h"

using namespace std;

//****************************************************************************************
int ReadNumberOfFragments(int fragmentLength, int overlappingResidues)
{

	int frLn, ovRd, tabPos1, tabPos2;
	string line;
	ifstream FragmentsFile;

	//Reads the interfaces list
	try
	{
		FragmentsFile.open(rootDir + "NumberOfFragments.txt");

		if (FragmentsFile.fail())
		{
			cout << "\nError reading the fragments information file.";
			return -1;
		}

	}
	catch (std::ifstream::failure &FileExcep)
	{
		cout << endl << FileExcep.what();
		return -1;
	}
	catch (...)
	{
		cout << "\nOther exception thrown.";
		return -1;
	}

	getline(FragmentsFile, line);	//To skip the column titles

	while (getline(FragmentsFile, line))
	{
		tabPos1 = line.find("\t", 0);
		frLn = atoi(line.substr(0, tabPos1).c_str());
		tabPos2 = line.find("\t", tabPos1 + 1);
		ovRd = atoi(line.substr(tabPos1 + 1, tabPos2).c_str());

		if (fragmentLength == frLn && overlappingResidues == ovRd)
		{
			FragmentsFile.close();
			return atoi(line.substr(tabPos2 + 1, line.length()).c_str());
		}
	}

	FragmentsFile.close();
	return 0;
}

//****************************************************************************************
int ReadNumberOfClusters(int fragmentLength, int overlappingResidues, int expectedMatchedPoints, int binSize)
{

	int frLn, ovRd, bnSz, mtSz, tabPos1, tabPos2;
	string line;
	ifstream clustersInfoFile;

	//Reads the interfaces list
	try
	{
		clustersInfoFile.open(rootDir + "NumberOfClusters.txt");

		if (clustersInfoFile.fail())
		{
			cout << "\nError reading the clusters information file.";
			return -1;
		}

	}
	catch (std::ifstream::failure &FileExcep)
	{
		cout << endl << FileExcep.what();
		return -1;
	}
	catch (...)
	{
		cout << "\nOther exception thrown.";
		return -1;
	}

	getline(clustersInfoFile, line);	//To skip the column titles

	while (getline(clustersInfoFile, line))
	{
		tabPos1 = line.find("\t", 0);
		frLn = atoi(line.substr(0, tabPos1).c_str());
		tabPos2 = line.find("\t", tabPos1 + 1);
		ovRd = atoi(line.substr(tabPos1 + 1, tabPos2).c_str());
		tabPos1 = line.find("\t", tabPos2 + 1);
		mtSz = atoi(line.substr(tabPos2 + 1, tabPos1).c_str());
		tabPos2 = line.find("\t", tabPos1 + 1);
		bnSz = atoi(line.substr(tabPos1 + 1, tabPos2).c_str());

		if (fragmentLength == frLn && overlappingResidues == ovRd && expectedMatchedPoints == mtSz && binSize == bnSz)
		{
			clustersInfoFile.close();
			return atoi(line.substr(tabPos2 + 1, line.length()).c_str());
		}
	}

	clustersInfoFile.close();
	return 0;
}


//****************************************************************************************
/*
Just to test the correctness of Geometric Hashing code
*/
void TestGeometricCalculations(void)
{

	/*		ATOM    305  CA  SER A  38 - 6.021  35.657   8.563  1.00 41.39   C	ATOM    311  CA  LEU A  39 - 9.078  37.570   9.943  1.00 40.77           C	ATOM    319  CA  ASN A  40 - 10.190  41.078  10.154  1.00 55.79           C	ATOM    327  CA  ALA A  41 - 12.706  42.490  12.323  1.00 64.31           C
	*/

	/*
	ATOM   4659  CA AALA B 133      45.293 -53.074 194.667  0.60 65.88           C  	ATOM   4669  CA AARG B 134      42.233 -55.298 194.504  0.60 73.34           C  	ATOM   4691  CA BPHE B 135      42.024 -56.189 193.815  0.40 72.72           C  
	*/

	
	/*int x1 = -12.706, y1 = 42.490, z1 = 12.323;
	int x2 = -9.078, y2 = 37.570, z2 = 9.943;
	int x3 = -10.190, y3 = 41.078, z3 = 10.154;*/

	/*float x1 = 45.293, y1 = -53.074, z1 = 194.667;
	float x2 = 42.233, y2 = -55.298, z2 = 194.504;
	float x3 = 42.024, y3 = -56.189, z3 = 193.815;*/

	float x1 = 0, y1 = 0, z1 = 0;
	float x2 = 0, y2 = 1, z2 = 1;
	float x3 = 0, y3 = 2, z3 = 2;

	TranslationParameter selectedRS;
	selectedRS = CalculateGeoTranslation(x3, y3, z3, x1, y1, z1, x2, y2, z2, 2);
	cout << "Selected three points to create the model are: (3,1,2)\n";
	
	if (selectedRS.Rx != 1000)
	{
		cout << "Point 2 is in: " << selectedRS.p2.x << "\t" << selectedRS.p2.y << "\t" << selectedRS.p2.z << "\n";
		cout << "Point 3 is in: " << selectedRS.p3.x << "\t" << selectedRS.p3.y << "\t" << selectedRS.p3.z << "\n";
	}

	//Point p;
	//p.x = -6.021, p.y = 35.657, p.z = 8.563;
	//p = CalculateNewPoint(selectedRS, p, 3);

	//cout << "Point 1 is in: " << p.x << "\t" << p.y << "\t" << p.z << "\n";
	

}

//*************************************************************************************
/*
This functions extracts some fragments from the fragments file which belong to the same cluster
and save them as text files. Then we can visualize them in Chimera to see how much similar they are.
*/

int ExtractFragmentsInClusters(list<int> fragNos, int fragmentLength, int overlappingResidues)
{

	list<int> fragmentNos = fragNos;

	//open file
	FILE * fragmentsFile;
	long fragmentsFileSize;
	long readsize = 0;
	char *fragmentsBuffer;
	FILE * destinationFile;
	ofstream clusteringFile;

	fragmentsBuffer = (char *)malloc(fragmentLength * 85);
	if (fragmentsBuffer == 0)
	{
		cout << "Error allocate the buffer.\r\n";
		return -1;
	}

	try
	{
		fragmentsFile = fopen((rootDir + "Protein_Fragments/FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "/ContinuousFragments.txt").c_str(), "r");
		if (fragmentsFile == NULL)
		{
			cout << "\nError reading the fragments file.";
			return -1;
		}
	}
	catch (ifstream::failure &FileExcep)
	{
		cout << FileExcep.what() << endl;
		return -1;
	}

	int fragmentNo = 0;
	while (!fragmentNos.empty()) 
	{

		fragmentNo = fragmentNos.front();

		fseek(fragmentsFile, (fragmentLength * 85) * fragmentNo, SEEK_SET);
		readsize = fread(fragmentsBuffer, 1, fragmentLength * 85, fragmentsFile);

		if (fragmentLength * 85 != readsize)
		{
			cout << "\nError reading the all of fragments file.";
			return -1;
		}


		try
		{
			clusteringFile.open((rootDir + "Random_Fragments/Fragment_" + to_string(fragmentNo) + ".txt").c_str(), ios::out | ios::app);
			if (clusteringFile.fail())
			{
				cout << "\nError creating the fragment information file.";
				return -1;
			}
			clusteringFile.write(fragmentsBuffer, (fragmentLength * 85));
		}
		catch (ifstream::failure &FileExcep)
		{
			cout << FileExcep.what() << endl;
			return -1;
		}
		catch (...)
		{
			cout << "\nOther exception thrown." << endl;
			return -1;
		}

		clusteringFile.close();

		fragmentNos.pop_front();
	}

	fclose(fragmentsFile);

	return 0;
}

//*************************************************************************************
void TestComparingHashTables()
{
	unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> hashTable1;
	unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> hashTable2;

	Point p;

	p.x = -7, p.y = -8, p.z = -9;
	AddToHashTable(hashTable1, p, 2, 5, 3);
	p.x = 2, p.y = 7, p.z = 3;
	AddToHashTable(hashTable1, p, 1, 2, 3);
	AddToHashTable(hashTable1, p, 2, 4, 3);
	p.x = 10, p.y = 8, p.z = -9;
	AddToHashTable(hashTable1, p, 2, 4, 3);
	p.x = 12, p.y = 17, p.z = 13;
	AddToHashTable(hashTable1, p, 5, 3, 4);
	p.x = 13, p.y = 17, p.z = 13;
	AddToHashTable(hashTable1, p, 2, 3, 3);

	//print the hash table
	cout << "-------------- Printing Hash Table ---------------" << endl;
	for (auto it1 = hashTable1.begin(); it1 != hashTable1.end(); ++it1)
	{
	cout << it1->first << ":";
	for (auto it2 = hashTable1[it1->first].begin(); it2 != hashTable1[it1->first].end(); ++it2)
	{
	cout << "   " << it2->first << ":";
	for (auto it3 = hashTable1[it1->first][it2->first].begin(); it3 != hashTable1[it1->first][it2->first].end(); ++it3)
	{
	cout << "   " << it3->first << ":";
	for (int it4 = 0; it4 < it3->second.size(); it4++)
	{
	cout << it3->second[it4] << " , ";
	}
	cout << endl;
	}
	cout << endl;
	}
	cout << endl;
	}
	
	p.x = 13, p.y = 17, p.z = 13;
	AddToHashTable(hashTable2, p, 5, 3, 4);
	p.x = -7, p.y = -8, p.z = -9;
	AddToHashTable(hashTable2, p, 5, 3, 4);
	p.x = 2, p.y = 7, p.z = 3;
	AddToHashTable(hashTable2, p, 5, 4, 3);
	
	//print the hash table
	cout << "-------------- Printing Hash Table ---------------" << endl;
	for (auto it1 = hashTable2.begin(); it1 != hashTable2.end(); ++it1)
	{
		cout << it1->first << ":";
		for (auto it2 = hashTable2[it1->first].begin(); it2 != hashTable2[it1->first].end(); ++it2)
		{
			cout << "   " << it2->first << ":";
			for (auto it3 = hashTable2[it1->first][it2->first].begin(); it3 != hashTable2[it1->first][it2->first].end(); ++it3)
			{
				cout << "   " << it3->first << ":";
				for (int it4 = 0; it4 < it3->second.size(); it4++)
				{
					cout << it3->second[it4] << " , ";
				}
				cout << endl;
			}
			cout << endl;
		}
		cout << endl;
	}
	
	if (CompareTwoHashTables(hashTable1, hashTable2, 2))
		cout << "\nHash tables are similar.";
	else
		cout << "\nHash tables are NOT similar.";


}


//*************************************************************************************
void PrintHashTable(unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>>& hashTable)
{
	
	cout << "\n-------------- Printing Hash Table ---------------" << endl;

	for (auto it1 = hashTable.begin(); it1 != hashTable.end(); ++it1)
	{
		//cout << it1->first << ":";
		for (auto it2 = hashTable[it1->first].begin(); it2 != hashTable[it1->first].end(); ++it2)
		{
			//cout << "   " << it2->first << ":";
			for (auto it3 = hashTable[it1->first][it2->first].begin(); it3 != hashTable[it1->first][it2->first].end(); ++it3)
			{
				cout << it1->first << ", " << it2->first << ", " << it3->first << "     ";
				for (int it4 = 0; it4 < it3->second.size(); it4++)
				{
					cout << it3->second[it4] << "  ";
				}
				cout << endl;
			}
			//cout << endl;
		}
		//cout << endl;
	}

}


