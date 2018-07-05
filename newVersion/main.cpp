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

#define numberOfInterfaces 45222


using namespace std;
using std::string;
string rootDir;
string slash;

int main(void)
{

	#ifdef linux 
		//linux code goes here
		rootDir = "/home/fhalakou/My_Thesis/Data/New_Continuous_Fragments_Results/";
		slash = "/";
	#elif _WIN32
		// windows code goes here
		rootDir = "D:/PhD/My_Thesis/Second_Step/Data/New_Continuous_Fragments_Results/";
		//rootDir = "N:\\Farideh\\";
		slash = "/";
	#else
		//default code goes here
		rootDir = "D:/PhD/My_Thesis/Second_Step/Data/New_Continuous_Fragments_Results/";
		slash = "/";
	#endif


	int fragmentLength = 4;
	int overlappingResidues = 3;
	int	expectedMatchedPoints = 3; //one point in each fragment is always at 0, 0, 0 coordinate
	int binSize = 2;	//used for clustering
	//int minNumberOfClusterMembers = 5;

	float topXresults = 0.1;

	int result = 0;
	int numberOfClusters, numberOfFragments;
	//int includeRatio = 100;

	//----------------------------Extract overlapping fragments from PRISM interfaces--------------------------------
	
	/*
	time_t startTime = time(0);
	cout << "\nExtracting fragments started... ";

	numberOfFragments = Extract_NonContinuous_Fragments(fragmentLength, overlappingResidues);

	time_t endTime = time(0);
	cout << "\nExtracting fragments finished.";
	cout << "\nElapsed time: " << to_string(difftime(endTime, startTime));
	*/
		
	//----------------------------Cluster the fragments----------------------------------------------------------------

	/*
	time_t startTime = time(0);
	cout << "\nClustering fragments started...";

	numberOfFragments = ReadNumberOfFragments(fragmentLength, overlappingResidues);
	cout << "\nNumber of fragments: " << numberOfFragments;
	numberOfClusters = ClusterFragments_v1_Parallel(fragmentLength, numberOfFragments, overlappingResidues, expectedMatchedPoints, binSize, minNumberOfClusterMembers);

	time_t endTime = time(0);
	cout << "\nClustering fragments finished.";
	cout << "\nElapsed time: " << to_string(difftime(endTime, startTime));
	*/

	//----------------------------Create the descriptive vectors for the PRISM interfaces--------------------------------

	/*
	time_t startTime = time(0);
	cout << "\nCreating interface descriptor vectors started...";
	

	//numberOfClusters = ReadNumberOfClusters(fragmentLength, overlappingResidues, expectedMatchedPoints, binSize);
	//cout << numberOfClusters;
	result = CreateInterfaceDescriptors_v2_Parallel(fragmentLength, overlappingResidues, binSize, expectedMatchedPoints);

	time_t endTime = time(0);
	cout << "\nCreating interface descriptor vectors finished.";
	cout << "\nElapsed time: " << to_string(difftime(endTime, startTime));
	*/
	
	//----------------------------Create protein surface descriptors and compare with interfaces--------------------------------
	
	
	// Running PRISM to extract the surfaces of the PDBs listed in the pair_list
	string jobName = "prismJob_6Angstrome";
	//string filename = "prism.py pair_list template_default " + jobName;
	//string command = "python ";
	//command += filename;
	//system(command.c_str());

	/*
	//numberOfClusters = ReadNumberOfClusters(fragmentLength, overlappingResidues, expectedMatchedPoints, binSize);
	//result = CreateProteinDescriptors_v2(jobName, fragmentLength, overlappingResidues, binSize, expectedMatchedPoints, numberOfClusters);
	//result = CompareProteinWithInterfaces_v2(fragmentLength, overlappingResidues, binSize, expectedMatchedPoints, includeRatio);
	*/


	time_t startTime = time(0);
	cout << "\nComparing protein surfaces with interfaces started...";
	
	//result = CompareProteinsWithInterfaces_v3(jobName, fragmentLength, overlappingResidues, expectedMatchedPoints, binSize);
	int numberOfResults = numberOfInterfaces * topXresults;
	//cout << endl << numberOfResults;
	result = CreateProteinDescriptorAndCompare(jobName, numberOfResults, fragmentLength, overlappingResidues, expectedMatchedPoints, binSize);
	
	time_t endTime = time(0);
	cout << "\nComparing protein surfaces with interfaces finished at " << time(&endTime);
	cout << "\nElapsed time: " << to_string(difftime(endTime, startTime)); 


	//--------------------------- Debugging -----------------------------------------
	
	//TestGeometricCalculations();
	//TestComparingHashTables();
	
	//numberOfClusters = ReadNumberOfClusters(fragmentLength, overlappingResidues, expectedMatchedPoints, binSize);
	//result = CalculateCosineDistance(fragmentLength, overlappingResidues, binSize, expectedMatchedPoints, numberOfInterfaces, numberOfClusters);
	
	//list<int> fragmentNos = {1193, 48937, 66192};
	//result = ExtractFragmentsInClusters(fragmentNos,fragmentLength, overlappingResidues);

	//--------------------------------------------------------------------------------


	getchar();
	return 0;
}
