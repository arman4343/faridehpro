#include <omp.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <unordered_map>
#include <ctime>
#include "My_Functions.h"


#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4503) 

using namespace std;
using std::string;
string rootDir;
string dirSeperator;

int main(void)
{

	#ifdef linux 
		//linux code goes here
		rootDir = "/scratch/users/fhalakou13/";
		dirSeperator = "/";
	#elif _WIN32
		// windows code goes here
		rootDir = "D:\\PhD\\My_Thesis\\Second_Step\\Data\\Continuous_Fragments_Results\\";
		dirSeperator = "\\";
	#else
		//default code goes here
		rootDir = "D:\\PhD\\My_Thesis\\Second_Step\\Data\\Continuous_Fragments_Results\\";
		dirSeperator = "\\";
	#endif

	int fragmentLength = 6;
	int overlappingResidues = 5;
	int	expectedMatchedPoints = 5; //one point in each fragment is always at 0, 0, 0 coordinate
	int binSize = 2;	//used for clustering
	int includeRatio = 70;
	int result = 0;
	int minNumberOfClusterMembers = 5;
	//int numberOfClusters, numberOfFragments;
	//int numberOfInterfaces = 45222;

	//----------------------------Extract overlapping fragments from PRISM interfaces--------------------------------
	
	/*
	time_t startTime = time(0);
	cout << "\nExtracting fragments started at " << time(&startTime);

	numberOfFragments = Extract_NonContinuous_Fragments(fragmentLength, overlappingResidues);

	time_t endTime = time(0);
	cout << "\nExtracting fragments finished at " << time(&endTime);
	cout << "\nElapsed time: " << to_string(difftime(endTime, startTime));
	*/
		
	//----------------------------Cluster the fragments--------------------------------

	/*
	time_t startTime = time(0);
	cout << "\nClustering fragments started at " << time(&startTime);

	numberOfFragments = ReadNumberOfFragments(fragmentLength, overlappingResidues);
	numberOfClusters = ClusterFragments_v1_Parallel(fragmentLength, numberOfFragments, overlappingResidues, expectedMatchedPoints, binSize, minNumberOfClusterMembers);

	time_t endTime = time(0);
	cout << "\nClustering fragments finished at " << time(&endTime);
	cout << "\nElapsed time: " << to_string(difftime(endTime, startTime));
	*/

	//----------------------------Create the descriptive vectors for the PRISM interfaces--------------------------------

	/*
	time_t startTime = time(0);
	cout << "\nCreating interface descriptor vectors started at " << time(&startTime);
	

	//numberOfClusters = ReadNumberOfClusters(fragmentLength, overlappingResidues, expectedMatchedPoints, binSize);
	//cout << numberOfClusters;
	result = CreateInterfaceDescriptors_v2_LessRotations(fragmentLength, overlappingResidues, binSize, expectedMatchedPoints, numberOfClusters);

	time_t endTime = time(0);
	cout << "\nCreating interface descriptor vectors finished at " << time(&endTime);
	cout << "\nElapsed time: " << to_string(difftime(endTime, startTime));
	*/
	
	//----------------------------Create the descriptive vectors for the protein surfaces--------------------------------
	
	/*
	// Running PRISM to extract the surfaces of the PDBs listed in the pair_list
	string jobName = "prismJob_6Angstrome";
	//string filename = "prism.py pair_list template_default " + jobName;
	//string command = "python ";
	//command += filename;
	//system(command.c_str());

	//Creating the descriptor vectors of the PDB surfaces
	time_t startTime = time(0);
	cout << "\nCreating descriptor vectors for proteins started at " << time(&startTime);

	numberOfClusters = ReadNumberOfClusters(fragmentLength, overlappingResidues, expectedMatchedPoints, binSize);
	result = CreateProteinDescriptors_v2(jobName, fragmentLength, overlappingResidues, binSize, expectedMatchedPoints, numberOfClusters);

	time_t endTime = time(0);
	cout << "\nCreating descriptor vectors for proteins finished at " << time(&endTime);
	cout << "\nElapsed time: " << to_string(difftime(endTime, startTime));
	*/

	//----------------------------Compare the protein descriptive vectors with interface descriptor vectors--------------------------------
	
	
	time_t startTime = time(0);
	cout << "\nComparing protein surfaces with interfaces started at " << time(&startTime);
	
	CompareProteinWithInterfaces_v2(fragmentLength, overlappingResidues, binSize, expectedMatchedPoints, includeRatio);
	
	time_t endTime = time(0);
	cout << "\nComparing protein surfaces with interfaces finished at " << time(&endTime);
	cout << "\nElapsed time: " << to_string(difftime(endTime, startTime)); 
	


	//--------------------------- Debugging -----------------------------------------
	
	//TestGeometricCalculations();
	
	//numberOfClusters = ReadNumberOfClusters(fragmentLength, overlappingResidues, expectedMatchedPoints, binSize);
	//result = CalculateCosineDistance(fragmentLength, overlappingResidues, binSize, expectedMatchedPoints, numberOfInterfaces, numberOfClusters);
	
	//list<int> fragmentNos = {174, 1945, 1038772};
	//result = ExtractFragments(fragmentNos,fragmentLength, overlappingResidues);

	//--------------------------------------------------------------------------------


	getchar();
	return 0;
}
