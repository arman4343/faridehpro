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
	string dataPath = "D:\\PhD\\My_Thesis\\Second_Step\\Data\\";

	//Reads the interfaces list
	try
	{
		FragmentsFile.open(dataPath + "NumberOfFragments.txt");

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

	string dataPath = "D:\\PhD\\My_Thesis\\Second_Step\\Data\\";
	//Reads the interfaces list
	try
	{
		clustersInfoFile.open(dataPath + "NumberOfClusters.txt");

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

	/*		ATOM    305  CA  SER A  38 - 6.021  35.657   8.563  1.00 41.39           C	ATOM    311  CA  LEU A  39 - 9.078  37.570   9.943  1.00 40.77           C	ATOM    319  CA  ASN A  40 - 10.190  41.078  10.154  1.00 55.79           C	ATOM    327  CA  ALA A  41 - 12.706  42.490  12.323  1.00 64.31           C
	*/

	cout << "Selected three points to create the model are: (4,2,3)\n";
	int x1 = -12.706, y1 = 42.490, z1 = 12.323;
	int x2 = -9.078, y2 = 37.570, z2 = 9.943;
	int x3 = -10.190, y3 = 41.078, z3 = 10.154;

	TranslationParameter selectedRS;
	selectedRS = CalculateGeoTranslation(x1, y1, z1, x2, y2, z2, x3, y3, z3, 3);
	cout << "Point 2 is in: " << selectedRS.p2.x << "\t" << selectedRS.p2.y << "\t" << selectedRS.p2.z << "\n";
	cout << "Point 3 is in: " << selectedRS.p3.x << "\t" << selectedRS.p3.y << "\t" << selectedRS.p3.z << "\n";

	Point p;
	p.x = -6.021, p.y = 35.657, p.z = 8.563;
	p = CalculateNewPoint(selectedRS, p, 3);

	cout << "Point 1 is in: " << p.x << "\t" << p.y << "\t" << p.z << "\n";

}

//*************************************************************************************
/*
This functions extracts some fragments from the fragments file which belong to the same cluster
and save them as text files. Then we can visualize them in Chimera to see how much similar they are.
*/

int ExtractFragments(std::list<int> fragNos, int fragmentLength, int overlappingResidues)
{

	std::list<int> fragmentNos = fragNos;

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
		fragmentsFile = fopen((rootDir + "Protein_Fragments" + dirSeperator + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + dirSeperator + "ContinuousFragments.txt").c_str(), "r");
		if (fragmentsFile == NULL)
		{
			cout << "\nError reading the fragments file.";
			return -1;
		}
	}
	catch (std::ifstream::failure &FileExcep)
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
			clusteringFile.open((rootDir + "Random_Fragments" + dirSeperator + "Fragment_" + to_string(fragmentNo) + ".txt").c_str(), ios::out | ios::app);
			if (clusteringFile.fail())
			{
				cout << "\nError creating the fragment information file.";
				return -1;
			}
			clusteringFile.write(fragmentsBuffer, (fragmentLength * 85));
		}
		catch (std::ifstream::failure &FileExcep)
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

