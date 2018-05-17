/*
This code creates the nonredundant library of small protein fragments from PRISM interface library.
Written by Farideh Halakou
*/

#include <omp.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <math.h>
#include <unordered_map>
#include "My_Functions.h"



#define PI 3.14159265

using namespace std;

#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4503) 
//****************************************************************************
/*
This process extracts the overlapping fragments from all the interfaces
and writes them to a file with fixed length of 85 characters for each line.
*/

int Extract_NonContinuous_Fragments(int fragmentLength, int overlappingResidues)
{
	int numberOfFragments = 0;
	string line;
	string *fragment;
	fragment = new string[fragmentLength];
	int lineInFragment = -1;
	ifstream interfacesFile;
	FILE * fragmentsFile;
	//bool continuousResidues;
	char wrt_buffer[85];


	//read the PRISM interfaces file
	try
	{
		interfacesFile.open("D:\\PhD\\prism\\prism_standalone\\template\\All_Prism_Interfaces.txt"); //All PRISM interfaces

		if (interfacesFile.fail())
		{
			cout << "\nError reading the interfaces file.";
			return -1;
		}

	}
	catch (std::ifstream::failure &FileExcep)
	{
		cout << FileExcep.what() << endl;
		return -1;
	}
	catch (...)
	{
		cout << "\nOther exception thrown.";
		return -1;
	}

	//Create the fragments file
	try
	{
		fragmentsFile = fopen((rootDir + dirSeperator + "NonContinuous_Fragments_Results" + dirSeperator + "Protein_Fragments" + dirSeperator + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + dirSeperator + "NonContinuousFragments.txt").c_str(), "wb");
		if (fragmentsFile == NULL)
		{
			cout << "\nError creating the fragments file.";
			return -1;
		}
	}
	catch (std::ifstream::failure &FileExcep)
	{
		cout << FileExcep.what() << endl;
		return -1;
	}
	catch (...)
	{
		cout << "\nOther exception thrown.";
		return -1;
	}


	while (getline(interfacesFile, line))	//each line has an interface name or residue information
	{

		//check if it is a line showin an interface name or not
		if (line.find("Interface Name") == 0)
		{
			lineInFragment = 0;
			//continuousResidues = true;
		}
		else
		{
			//cout << "\n" << lineInFragment;
			fragment[lineInFragment] = line;
			lineInFragment++;
			
			if (lineInFragment == fragmentLength)
			{
				for (int i = 0; i < fragmentLength; i++)
				{
					memset(wrt_buffer, 0, 85);
					memcpy(wrt_buffer, fragment[i].c_str(), fragment[i].size());
					if (fwrite(wrt_buffer, sizeof(char), 85, fragmentsFile)<85)
					{
						cout << "Error writing in file.\n";
						return -1;
					}
				}

				numberOfFragments++;

				lineInFragment = overlappingResidues;
				for (int i = 0; i < overlappingResidues; i++)
					fragment[i] = fragment[fragmentLength - overlappingResidues + i];

			}
		}
	}

	interfacesFile.close();
	fclose(fragmentsFile);
	cout << "\nNumber of fragments: " << numberOfFragments;
	return numberOfFragments;

}


//****************************************************************************
/*
This process extracts the overlapping fragments from all the interfaces 
and writes them to a file with fixed length of 85 characters for each line. 
*/

int ExtractFragments_v1(int fragmentLength, int overlappingResidues)
{

	int numberOfFragments = 0;
	string line;
	string *fragment;
	fragment = new string[fragmentLength];
	int lineInFragment = -1;
	ifstream interfacesFile;
	FILE * fragmentsFile;
	bool continuousResidues;
	char wrt_buffer[85];


	//read the PRISM interfaces file
	try
	{
		interfacesFile.open("D:\\PhD\\prism\\prism_standalone\\template\\tempAll_Prism_Interfaces.txt"); //All PRISM interfaces

		if (interfacesFile.fail())
		{
			cout << "\nError reading the interfaces file.";
			return -1;
		}

	}
	catch (std::ifstream::failure &FileExcep)
	{
		cout << FileExcep.what() << endl;
		return -1;
	}
	catch (...)
	{
		cout << "\nOther exception thrown.";
		return -1;
	}

	//Create the fragments file
	try
	{
		fragmentsFile = fopen(("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Protein_Fragments\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "\\fragments.txt").c_str(), "wb");
		if (fragmentsFile == NULL)
		{
			cout << "\nError creating the fragments file.";
			return -1;
		}
	}
	catch (std::ifstream::failure &FileExcep)
	{
		cout << FileExcep.what() << endl;
		return -1;
	}
	catch (...)
	{
		cout << "\nOther exception thrown.";
		return -1;
	}


	while (getline(interfacesFile, line))	//each line has an interface name or residue information
	{

		//check if it is a line showin an interface name or not
		if (line.find("Interface Name") == 0)
		{
			lineInFragment = 0;
			continuousResidues = true;
		}
		else
		{
			//cout << "\n" << lineInFragment;
			fragment[lineInFragment] = line;


			if (lineInFragment == 0)
				lineInFragment++;
			else if (lineInFragment > 0)
			{
				int previousResidueNumber = atoi(fragment[lineInFragment - 1].substr(22, 5).c_str());
				//cout << endl << previousResidueNumber;
				int nextResidueNumber = atoi(fragment[lineInFragment].substr(22, 5).c_str());
				//cout << endl << nextResidueNumber;

				if (previousResidueNumber < nextResidueNumber - 1)
				{
					continuousResidues = false;
					fragment[0] = fragment[lineInFragment];
					lineInFragment = 1;
				}
				else
					lineInFragment++;
			}

			if (lineInFragment == fragmentLength)	//If it is a continuous fragment of size fragmentLength
			{
				for (int i = 0; i < fragmentLength; i++)
				{
					memset(wrt_buffer, 0, 85);
					memcpy(wrt_buffer, fragment[i].c_str(), fragment[i].size());
					if (fwrite(wrt_buffer, sizeof(char), 85, fragmentsFile)<85)
					{
						cout << "Error writing on file \r\n";
						return -1;
					}
				}


				numberOfFragments++;

				lineInFragment = overlappingResidues;
				for (int i = 0; i < overlappingResidues; i++)
					fragment[i] = fragment[fragmentLength - overlappingResidues + i];

			}
		}
	}
	interfacesFile.close();
	fclose(fragmentsFile);
	cout << "\nNumber of fragments: " << numberOfFragments;
	return numberOfFragments;
}

//****************************************************************************
/*
This process randomly selects a seed fragment as a new cluster and finds,
in a parallel manner, the similar fragments to the seed. It continues till
all the fragments are assigned to a cluster.
This version reads the whole fragments file to the RAM.
*/

int ClusterFragments_v1_Parallel(int fragmentLength, int numberOfFragments, int overlappingResidues, int expectedMatchedPoints, int binSize, int minNumberOfClusterMembers)
{

	int numberOfClusters = 0;
	int numberOfSimilarFragments = 0;
	string *lines;
	lines = new string[fragmentLength];
	string similarFragments, line;
	int seedFragmentNumber;
	//char rd_buffer[85];
	string dataPath = "D:\\PhD\\My_Thesis\\Second_Step\\Data\\";

	unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> seedHashTable;
	//unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> nonSeedHashTable;
	vector <int> unclusteredFragments;	//A vector of fragments which are all unclustered yet
	vector <int> remainedUnclustered;

	FILE * fragmentsFile;
	long fragmentsFileSize;
	long readsize = 0;
	char *fragmentsBuffer;
	ofstream clusteringFile;
	fstream unclusteredFragmentsFile;


	//Reads the unclustered fragments file or creates it if doesn't exist
	try
	{
		unclusteredFragmentsFile.open((dataPath + "Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\UnclusteredFragments.txt").c_str(), fstream::in | fstream::out | fstream::app);
		if (unclusteredFragmentsFile.fail())	//If the file doesn't exist it means that it is the first run do we create the list including all fragments
		{
			cout << "\nError creting or reading the unclustered fragments file.";
			return -1;
		}
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


	//Reads the fragments file
	try
	{
		fragmentsFile = fopen((dataPath + "Protein_Fragments\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "\\fragments.txt").c_str(), "r");
		if (fragmentsFile == NULL)
		{
			cout << "\nError reading the fragments file.";
			return -1;
		}
		//Read the file to RAM
		else
		{

			fseek(fragmentsFile, 0, SEEK_END);
			fragmentsFileSize = ftell(fragmentsFile);
			fragmentsBuffer = (char *)malloc(fragmentsFileSize + 1);
			fseek(fragmentsFile, 0, SEEK_SET);
			readsize = fread(fragmentsBuffer, 1, fragmentsFileSize, fragmentsFile);
			fclose(fragmentsFile);

			if (readsize != fragmentsFileSize)
			{
				cout << "\nError reading the all of fragments file.";
				return -1;
			}
		}
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



	//Reads the clustering information file or creates it if doesn't exist
	try
	{
		clusteringFile.open((dataPath + "Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\ClusteringLog.txt").c_str(), ios::out | ios::app);
		if (clusteringFile.fail())
		{
			cout << "\nError creating the clustering information file.";
			return -1;
		}
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



	if (unclusteredFragmentsFile.peek() == EOF)	//if the file is empty, so it is the first run and all fragments are not clustered yet
	{

		for (int i = 0; i < numberOfFragments; i++)
		{
			unclusteredFragments.push_back(i);
			//unclusteredFragmentsFile << i << "\t";
		}
	}
	else
	{

		while (getline(unclusteredFragmentsFile, line));	//Read till the last line which has the updated unclustered fragments list	
		stringstream stream(line);
		string temp;

		stream >> temp;
		stream >> temp;
		int fragmentNum;
		while (stream >> fragmentNum)
		{
			unclusteredFragments.push_back(fragmentNum);
			//cout << fragmentNum;
			//getchar();
		}
		unclusteredFragmentsFile << "\n";
	}

	unclusteredFragmentsFile.close();



	srand((unsigned int)time(NULL));
	char *lOprationPosition;

	//Clustering the fragments
	while (unclusteredFragments.size() > 0)
	{
		int seed = rand() % unclusteredFragments.size();	//Randomly select one fragment to be the seed
		//seed = 743643;
		seedFragmentNumber = unclusteredFragments[seed];
		cout << "\nSeed fragment: " << seedFragmentNumber << "\t";
		//cout << "\nunclusteredFragments.size: " << unclusteredFragments.size() << "\n";
		similarFragments = " ";
		numberOfSimilarFragments = 0;

		unclusteredFragments.erase(unclusteredFragments.begin() + seed);

		int positionInFile = seedFragmentNumber * (fragmentLength * 85);
		//cout << "\nPosition in File: " << positionInFile;
		lOprationPosition = fragmentsBuffer + positionInFile;
		for (int index = 0; index < fragmentLength; index++)
		{
			lines[index].assign(lOprationPosition, lOprationPosition + 85);
			lOprationPosition += 85;
		}

		seedHashTable.clear();

		//Select three consecutive points as Reference Set (RS) and create the hashtable
		for (int c = 0; c < fragmentLength - 2; c++)
			for (int i = c; i < c + 3; i++)
				for (int j = c; j < c + 3; j++)
					for (int k = c; k < c + 3; k++)
						if (i != j && i != k && j != k)	//three different points to create RS
						{

							float x1 = (float)atof(lines[i].substr(30, 8).c_str());
							float y1 = (float)atof(lines[i].substr(38, 8).c_str());
							float z1 = (float)atof(lines[i].substr(46, 8).c_str());

							float x2 = (float)atof(lines[j].substr(30, 8).c_str());
							float y2 = (float)atof(lines[j].substr(38, 8).c_str());
							float z2 = (float)atof(lines[j].substr(46, 8).c_str());

							float x3 = (float)atof(lines[k].substr(30, 8).c_str());
							float y3 = (float)atof(lines[k].substr(38, 8).c_str());
							float z3 = (float)atof(lines[k].substr(46, 8).c_str());

							TranslationParameter selectedRS;
							selectedRS = CalculateGeoTranslation(x1, y1, z1, x2, y2, z2, x3, y3, z3, binSize);


							if (selectedRS.Rx != 0 || selectedRS.Ry != 0 || selectedRS.Rz != 0)
							{
								//Add the second and third points of RS to the hash table
								AddToHashTable(seedHashTable, selectedRS.p2, i + 1, j + 1, k + 1);
								AddToHashTable(seedHashTable, selectedRS.p3, i + 1, j + 1, k + 1);

								//Add non-RS points of the fragment to the hash table
								for (int f = 0; f < fragmentLength; f++)
								{
									if (f != i && f != j && f != k)
									{
										Point p;
										p.x = (float)atof(lines[f].substr(30, 8).c_str());
										p.y = (float)atof(lines[f].substr(38, 8).c_str());
										p.z = (float)atof(lines[f].substr(46, 8).c_str());
										p = CalculateNewPoint(selectedRS, p, binSize);
										AddToHashTable(seedHashTable, p, i + 1, j + 1, k + 1);
									}
								}
							}
						}

		

		//Compare the seed hashtable with all other unclustered fragments hashtables

		int fragmentNoGlobal = 0;
		time_t t = time(0);

#pragma omp parallel //num_threads(16)
		{
			//private variabales
			int fragmentNo = 0;
			unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> nonSeedHashTable2;
			int index2;
			bool matched2;
			int c2, i2, j2, k2, f2;
			char *lOprationPosition2;
			string *lines2;
			lines2 = new string[fragmentLength];
			int quit = 0;

			while (true)
			{
				//Set the fragmentNo in critical section
#pragma omp critical
				{
					//printf(" world(%d) \n", omp_get_num_threads());

					fragmentNo = fragmentNoGlobal;
					fragmentNoGlobal++;
					if (time(0) - t > 1)
					{
						cout << "unclusteredFragments.size: " << unclusteredFragments.size() << "\n" << "fragmentNo" << fragmentNo << "\n";
						t = time(0);
					}
				}


				if (fragmentNo >= (int)unclusteredFragments.size())
				{
					break;
				}

				int positionInFile2 = unclusteredFragments[fragmentNo] * (fragmentLength * 85);
				//cout << endl << "Fragment no: " << unclusteredFragments[fragmentNo] << "\tPosition in file: " << positionInFile;
				lOprationPosition2 = fragmentsBuffer + positionInFile2;
				for (index2 = 0; index2 < fragmentLength; index2++)
				{
					lines2[index2].assign(lOprationPosition2, lOprationPosition2 + 85);
					lOprationPosition2 += 85;
				}


				matched2 = false;	//if the fragment hashtable and seed hash table matches
				for (c2 = 0; c2 < fragmentLength - 2; c2++)	// we select three consecutive points as Reference Set (RS)
				{
					for (i2 = c2; i2 < c2 + 3; i2++)
					{
						for (j2 = c2; j2 < c2 + 3; j2++)
						{
							for (k2 = c2; k2 < c2 + 3; k2++)
								if (i2 != j2 && i2 != k2 && j2 != k2)	//three different points to create RS
								{

									nonSeedHashTable2.clear();

									float x1Private = (float)atof(lines2[i2].substr(30, 8).c_str());
									float y1Private = (float)atof(lines2[i2].substr(38, 8).c_str());
									float z1Private = (float)atof(lines2[i2].substr(46, 8).c_str());

									float x2Private = (float)atof(lines2[j2].substr(30, 8).c_str());
									float y2Private = (float)atof(lines2[j2].substr(38, 8).c_str());
									float z2Private = (float)atof(lines2[j2].substr(46, 8).c_str());

									float x3Private = (float)atof(lines2[k2].substr(30, 8).c_str());
									float y3Private = (float)atof(lines2[k2].substr(38, 8).c_str());
									float z3Private = (float)atof(lines2[k2].substr(46, 8).c_str());

									TranslationParameter selectedRS;
									selectedRS = CalculateGeoTranslation(x1Private, y1Private, z1Private, x2Private, y2Private, z2Private, x3Private, y3Private, z3Private, binSize);

									if (selectedRS.Rx != 0 || selectedRS.Ry != 0 || selectedRS.Rz != 0)
									{
										AddToHashTable(nonSeedHashTable2, selectedRS.p2, i2 + 1, j2 + 1, k2 + 1);
										AddToHashTable(nonSeedHashTable2, selectedRS.p3, i2 + 1, j2 + 1, k2 + 1);

										for (f2 = 0; f2 < fragmentLength; f2++)
										{
											if (f2 != i2 && f2 != j2 && f2 != k2)
											{
												Point p;
												p.x = (float)atof(lines2[f2].substr(30, 8).c_str());
												p.y = (float)atof(lines2[f2].substr(38, 8).c_str());
												p.z = (float)atof(lines2[f2].substr(46, 8).c_str());
												p = CalculateNewPoint(selectedRS, p, binSize);
												AddToHashTable(nonSeedHashTable2, p, i2 + 1, j2 + 1, k2 + 1);
											}
										}

										//compares two hash tables
										matched2 = CompareTwoHashTables(seedHashTable, nonSeedHashTable2, expectedMatchedPoints);
										if (matched2) // if the seed fragment and the new one are similar
											break;
									}
									else {
										//cout << x1Private << y1Private << z1Private << x2Private << y2Private << z2Private << x3Private << y3Private << z3Private << "\r\n";
									}
								}
							if (matched2) break;
						}// for j
						if (matched2) break;
					}// for i
					if (matched2) break;
				}//for c


#pragma omp critical
				{
					if (matched2)
					{
						numberOfSimilarFragments++;
						similarFragments += to_string(unclusteredFragments[fragmentNo]) + "\t";
						//cout << "\n" << similarFragments;
					}
					else
					{
						remainedUnclustered.push_back(unclusteredFragments[fragmentNo]);
					}
				}
			}
		}

		cout << "Number of cluster members: " << numberOfSimilarFragments + 1;

		unclusteredFragments.clear();
		unclusteredFragments = remainedUnclustered;
		remainedUnclustered.clear();

		unclusteredFragmentsFile.open((dataPath + "Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\UnclusteredFragments.txt").c_str(), fstream::out | fstream::app);
		unclusteredFragmentsFile << "SeedFragment:" << seedFragmentNumber << "\t NumberOfUnclusteredFragments:" << unclusteredFragments.size();
		for (int i = 0; i < (int)unclusteredFragments.size(); i++)	//Writes the unclustered fragments to file
			unclusteredFragmentsFile << "\t" << unclusteredFragments[i];
		unclusteredFragmentsFile << "\n";
		unclusteredFragmentsFile.close();


		if (numberOfSimilarFragments + 1 >= minNumberOfClusterMembers)
		{
			clusteringFile << "\nCluster_" << ++numberOfClusters << "\t" << seedFragmentNumber << "\t" << similarFragments;
			int res = WriteHashTableToFile((dataPath + "Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\Clusters.txt").c_str(), numberOfClusters, seedHashTable);
			if (res == -1)
				break;
		}
	}

	clusteringFile.close();
	cout << "\nNumber of clusters: " << numberOfClusters;
	return numberOfClusters;

}


//********************************************************************************
/*
This process randomly selects a seed fragment as a new cluster and finds,
in a parallel manner, the similar fragments to the seed. It continues till
all the fragments are assigned to a cluster. It reads the whole fragments 
file to the RAM.
This version has less rotations (1,2,3 and 3,2,1).
*/

int ClusterFragments_v2_Parallel(int fragmentLength, int numberOfFragments, int overlappingResidues, int expectedMatchedPoints, int binSize, int minNumberOfClusterMembers)
{

	int numberOfClusters = 0;
	int numberOfSimilarFragments = 0;
	string *lines;
	lines = new string[fragmentLength];
	string similarFragments, line;
	int seedFragmentNumber;
	//char rd_buffer[85];

	unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> seedHashTable;
	//unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> nonSeedHashTable;
	vector <int> unclusteredFragments;	//A vector of fragments which are all unclustered yet
	vector <int> remainedUnclustered;

	FILE * fragmentsFile;
	long fragmentsFileSize;
	long readsize = 0;
	char *fragmentsBuffer;
	ofstream clusteringFile;
	fstream unclusteredFragmentsFile;
	string dataPath = "NonContinuous_Fragments_Results\\";


	//Reads the unclustered fragments file or creates it if doesn't exist
	try
	{
		unclusteredFragmentsFile.open((rootDir + dataPath + "Less_Rotations_Clustering_Results" + dirSeperator + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + dirSeperator + "UnclusteredFragments.txt").c_str(), fstream::in | fstream::out | fstream::app);
		if (unclusteredFragmentsFile.fail())	//If the file doesn't exist it means that it is the first run do we create the list including all fragments
		{
			cout << "\nError creting or reading the unclustered fragments file.";
			return -1;
		}
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


	//Reads the fragments file
	try
	{
		fragmentsFile = fopen((rootDir + "Protein_Fragments" + dirSeperator + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + dirSeperator + "fragments.txt").c_str(), "r");
		if (fragmentsFile == NULL)
		{
			cout << "\nError reading the fragments file.";
			return -1;
		}
		//Read the file to RAM
		else
		{
			fseek(fragmentsFile, 0, SEEK_END);
			fragmentsFileSize = ftell(fragmentsFile);
			fragmentsBuffer = (char *)malloc(fragmentsFileSize + 1);
			fseek(fragmentsFile, 0, SEEK_SET);
			readsize = fread(fragmentsBuffer, 1, fragmentsFileSize, fragmentsFile);
			fclose(fragmentsFile);
			if (readsize != fragmentsFileSize)
			{
				cout << "\nError reading the all of fragments file.";
				return -1;
			}
		}
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



	//Reads the clustering information file or creates it if doesn't exist
	try
	{
		clusteringFile.open((rootDir + "Clustering_Results" + dirSeperator + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + dirSeperator + "ClusteringLog.txt") , ios::out | ios::app);
		if (clusteringFile.fail())
		{
			cout << "\nError creating the clustering information file.";
			return -1;
		}
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



	if (unclusteredFragmentsFile.peek() == EOF)	//if the file is empty, so it is the first run and all fragments are not clustered yet
	{

		for (int i = 0; i < numberOfFragments; i++)
		{
			unclusteredFragments.push_back(i);
			//unclusteredFragmentsFile << i << "\t";
		}
	}
	else
	{

		while (getline(unclusteredFragmentsFile, line));	//Read till the last line which has the updated unclustered fragments list	
		stringstream stream(line);
		string temp;

		stream >> temp;
		stream >> temp;
		int fragmentNum;
		while (stream >> fragmentNum)
		{
			unclusteredFragments.push_back(fragmentNum);
			//cout << fragmentNum;
			//getchar();
		}
		unclusteredFragmentsFile << "\n";
	}

	unclusteredFragmentsFile.close();



	srand((unsigned int)time(NULL));
	char *lOprationPosition;

	//Clustering the fragments
	while (unclusteredFragments.size() > 0)
	{
		int seed = rand() % unclusteredFragments.size();	//Randomly select one fragment to be the seed
		//seed = 743643;
		seedFragmentNumber = unclusteredFragments[seed];
		cout << "\nSeed fragment: " << seedFragmentNumber << "\t";
		//cout << "\nunclusteredFragments.size: " << unclusteredFragments.size() << "\n";
		similarFragments = " ";
		numberOfSimilarFragments = 0;

		unclusteredFragments.erase(unclusteredFragments.begin() + seed);

		int positionInFile = seedFragmentNumber * (fragmentLength * 85);
		//cout << "\nPosition in File: " << positionInFile;
		lOprationPosition = fragmentsBuffer + positionInFile;
		for (int index = 0; index < fragmentLength; index++)
		{
			lines[index].assign(lOprationPosition, lOprationPosition + 85);
			lOprationPosition += 85;
		}

		seedHashTable.clear();

		//Select three consecutive points as Reference Set (RS) and create the hashtable
		for (int c = 0; c < fragmentLength - 2; c++)
			/*for (int i = c; i < c + 3; i++)
			for (int j = c; j < c + 3; j++)
			for (int k = c; k < c + 3; k++)
			if (i != j && i != k && j != k)	//three different points to create RS*/
		{
			int i = c;
			int j = c + 1;
			int k = c + 2;
			float x1 = (float)atof(lines[i].substr(30, 8).c_str());
			float y1 = (float)atof(lines[i].substr(38, 8).c_str());
			float z1 = (float)atof(lines[i].substr(46, 8).c_str());


			float x2 = (float)atof(lines[j].substr(30, 8).c_str());
			float y2 = (float)atof(lines[j].substr(38, 8).c_str());
			float z2 = (float)atof(lines[j].substr(46, 8).c_str());

			float x3 = (float)atof(lines[k].substr(30, 8).c_str());
			float y3 = (float)atof(lines[k].substr(38, 8).c_str());
			float z3 = (float)atof(lines[k].substr(46, 8).c_str());

			TranslationParameter selectedRS;
			selectedRS = CalculateGeoTranslation(x1, y1, z1, x2, y2, z2, x3, y3, z3, binSize);


			if (selectedRS.Rx != 0 || selectedRS.Ry != 0 || selectedRS.Rz != 0)
			{
				//Add the second and third points of RS to the hash table
				AddToHashTable(seedHashTable, selectedRS.p2, i + 1, j + 1, k + 1);
				AddToHashTable(seedHashTable, selectedRS.p3, i + 1, j + 1, k + 1);

				//Add non-RS points of the fragment to the hash table
				for (int f = 0; f < fragmentLength; f++)
				{
					if (f != i && f != j && f != k)
					{
						Point p;
						p.x = (float)atof(lines[f].substr(30, 8).c_str());
						p.y = (float)atof(lines[f].substr(38, 8).c_str());
						p.z = (float)atof(lines[f].substr(46, 8).c_str());
						p = CalculateNewPoint(selectedRS, p, binSize);
						AddToHashTable(seedHashTable, p, i + 1, j + 1, k + 1);
					}
				}
			}

			x1 = (float)atof(lines[k].substr(30, 8).c_str());
			y1 = (float)atof(lines[k].substr(38, 8).c_str());
			z1 = (float)atof(lines[k].substr(46, 8).c_str());

			x2 = (float)atof(lines[j].substr(30, 8).c_str());
			y2 = (float)atof(lines[j].substr(38, 8).c_str());
			z2 = (float)atof(lines[j].substr(46, 8).c_str());

			x3 = (float)atof(lines[i].substr(30, 8).c_str());
			y3 = (float)atof(lines[i].substr(38, 8).c_str());
			z3 = (float)atof(lines[i].substr(46, 8).c_str());

			//TranslationParameter selectedRS;
			selectedRS = CalculateGeoTranslation(x1, y1, z1, x2, y2, z2, x3, y3, z3, binSize);


			if (selectedRS.Rx != 0 || selectedRS.Ry != 0 || selectedRS.Rz != 0)
			{
				//Add the second and third points of RS to the hash table
				AddToHashTable(seedHashTable, selectedRS.p2, i + 1, j + 1, k + 1);
				AddToHashTable(seedHashTable, selectedRS.p3, i + 1, j + 1, k + 1);

				//Add non-RS points of the fragment to the hash table
				for (int f = 0; f < fragmentLength; f++)
				{
					if (f != i && f != j && f != k)
					{
						Point p;
						p.x = (float)atof(lines[f].substr(30, 8).c_str());
						p.y = (float)atof(lines[f].substr(38, 8).c_str());
						p.z = (float)atof(lines[f].substr(46, 8).c_str());
						p = CalculateNewPoint(selectedRS, p, binSize);
						AddToHashTable(seedHashTable, p, i + 1, j + 1, k + 1);
					}
				}
			}
		}

		

		//Compare the seed hashtable with all other unclustered fragments hashtables

		int fragmentNoGlobal = 0;
		time_t t = time(0);

#pragma omp parallel //num_threads(16)
		{
			//private variabales
			int fragmentNo = 0;
			unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> nonSeedHashTable2;
			int index2;
			bool matched2;
			int c2, i2, j2, k2, f2;
			char *lOprationPosition2;
			string *lines2;
			lines2 = new string[fragmentLength];
			int quit = 0;

			while (true)
			{
				//Set the fragmentNo in critical section
#pragma omp critical
				{
					//printf(" world(%d) \n", omp_get_num_threads());

					fragmentNo = fragmentNoGlobal;
					fragmentNoGlobal++;
					if (time(0) - t > 1)
					{
						cout << "unclusteredFragments.size: " << unclusteredFragments.size() << "\n" << "fragmentNo" << fragmentNo << "\n";
						t = time(0);
					}
				}


				if (fragmentNo >= (int)unclusteredFragments.size())
				{
					break;
				}

				int positionInFile2 = unclusteredFragments[fragmentNo] * (fragmentLength * 85);
				//cout << endl << "Fragment no: " << unclusteredFragments[fragmentNo] << "\tPosition in file: " << positionInFile;
				lOprationPosition2 = fragmentsBuffer + positionInFile2;
				for (index2 = 0; index2 < fragmentLength; index2++)
				{
					lines2[index2].assign(lOprationPosition2, lOprationPosition2 + 85);
					lOprationPosition2 += 85;
				}


				matched2 = false;	//if the fragment hashtable and seed hash table matches
				for (c2 = 0; c2 < fragmentLength - 2; c2++)	// we select three consecutive points as Reference Set (RS)
				{
	
					i2 = c2;
					j2 = c2 + 1;
					k2 = c2 + 2;

					nonSeedHashTable2.clear();

					float x1Private = (float)atof(lines2[i2].substr(30, 8).c_str());
					float y1Private = (float)atof(lines2[i2].substr(38, 8).c_str());
					float z1Private = (float)atof(lines2[i2].substr(46, 8).c_str());

					float x2Private = (float)atof(lines2[j2].substr(30, 8).c_str());
					float y2Private = (float)atof(lines2[j2].substr(38, 8).c_str());
					float z2Private = (float)atof(lines2[j2].substr(46, 8).c_str());

					float x3Private = (float)atof(lines2[k2].substr(30, 8).c_str());
					float y3Private = (float)atof(lines2[k2].substr(38, 8).c_str());
					float z3Private = (float)atof(lines2[k2].substr(46, 8).c_str());

					TranslationParameter selectedRS;
					selectedRS = CalculateGeoTranslation(x1Private, y1Private, z1Private, x2Private, y2Private, z2Private, x3Private, y3Private, z3Private, binSize);

					if (selectedRS.Rx != 0 || selectedRS.Ry != 0 || selectedRS.Rz != 0)
					{
						AddToHashTable(nonSeedHashTable2, selectedRS.p2, i2 + 1, j2 + 1, k2 + 1);
						AddToHashTable(nonSeedHashTable2, selectedRS.p3, i2 + 1, j2 + 1, k2 + 1);

						for (f2 = 0; f2 < fragmentLength; f2++)
						{
							if (f2 != i2 && f2 != j2 && f2 != k2)
							{
								Point p;
								p.x = (float)atof(lines2[f2].substr(30, 8).c_str());
								p.y = (float)atof(lines2[f2].substr(38, 8).c_str());
								p.z = (float)atof(lines2[f2].substr(46, 8).c_str());
								p = CalculateNewPoint(selectedRS, p, binSize);
								AddToHashTable(nonSeedHashTable2, p, i2 + 1, j2 + 1, k2 + 1);
							}
						}

						//compares two hash tables
						matched2 = CompareTwoHashTables(seedHashTable, nonSeedHashTable2, expectedMatchedPoints);
						if (matched2) // if the seed fragment and the new one are similar
							break;
					}

					nonSeedHashTable2.clear();

					x1Private = (float)atof(lines2[k2].substr(30, 8).c_str());
					y1Private = (float)atof(lines2[k2].substr(38, 8).c_str());
					z1Private = (float)atof(lines2[k2].substr(46, 8).c_str());

					x2Private = (float)atof(lines2[j2].substr(30, 8).c_str());
					y2Private = (float)atof(lines2[j2].substr(38, 8).c_str());
					z2Private = (float)atof(lines2[j2].substr(46, 8).c_str());

					x3Private = (float)atof(lines2[i2].substr(30, 8).c_str());
					y3Private = (float)atof(lines2[i2].substr(38, 8).c_str());
					z3Private = (float)atof(lines2[i2].substr(46, 8).c_str());

					//selectedRS;
					selectedRS = CalculateGeoTranslation(x1Private, y1Private, z1Private, x2Private, y2Private, z2Private, x3Private, y3Private, z3Private, binSize);

					if (selectedRS.Rx != 0 || selectedRS.Ry != 0 || selectedRS.Rz != 0)
					{
						AddToHashTable(nonSeedHashTable2, selectedRS.p2, i2 + 1, j2 + 1, k2 + 1);
						AddToHashTable(nonSeedHashTable2, selectedRS.p3, i2 + 1, j2 + 1, k2 + 1);

						for (f2 = 0; f2 < fragmentLength; f2++)
						{
							if (f2 != i2 && f2 != j2 && f2 != k2)
							{
								Point p;
								p.x = (float)atof(lines2[f2].substr(30, 8).c_str());
								p.y = (float)atof(lines2[f2].substr(38, 8).c_str());
								p.z = (float)atof(lines2[f2].substr(46, 8).c_str());
								p = CalculateNewPoint(selectedRS, p, binSize);
								AddToHashTable(nonSeedHashTable2, p, i2 + 1, j2 + 1, k2 + 1);
							}
						}

						//compares two hash tables
						matched2 = CompareTwoHashTables(seedHashTable, nonSeedHashTable2, expectedMatchedPoints);
						if (matched2) // if the seed fragment and the new one are similar
							break;
					}
				}

#pragma omp critical
				{
					if (matched2)
					{
						numberOfSimilarFragments++;
						similarFragments += to_string(unclusteredFragments[fragmentNo]) + "\t";
						//cout << "\n" << similarFragments;
					}
					else
					{
						remainedUnclustered.push_back(unclusteredFragments[fragmentNo]);
					}
				}
			}
		}

		cout << "Number of cluster members: " << numberOfSimilarFragments + 1;

		unclusteredFragments.clear();
		unclusteredFragments = remainedUnclustered;
		remainedUnclustered.clear();

		unclusteredFragmentsFile.open((rootDir + "Clustering_Results" + dirSeperator + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + dirSeperator + "UnclusteredFragments.txt").c_str(), fstream::out | fstream::app);
		unclusteredFragmentsFile << "SeedFragment:" << seedFragmentNumber << "\t NumberOfUnclusteredFragments:" << unclusteredFragments.size();
		for (int i = 0; i < (int)unclusteredFragments.size(); i++)	//Writes the unclustered fragments to file
			unclusteredFragmentsFile << "\t" << unclusteredFragments[i];
		unclusteredFragmentsFile << "\n";
		unclusteredFragmentsFile.close();


		if (numberOfSimilarFragments + 1 >= minNumberOfClusterMembers)
		{
			clusteringFile << "\nCluster_" << ++numberOfClusters << "\t" << seedFragmentNumber << "\t" << similarFragments;
			int res = WriteHashTableToFile((rootDir + "Clustering_Results" + dirSeperator + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + dirSeperator + "Clusters.txt").c_str(), numberOfClusters, seedHashTable);
			if (res == -1)
				break;
		}
	}

	clusteringFile.close();
	cout << "\nNumber of clusters: " << numberOfClusters;
	return numberOfClusters;

}



//****************************************************************************
/* This process writes a hash table into a file */
int WriteHashTableToFile(string fname, int numberOfClusters, unordered_map<float, std::unordered_map<float, std::unordered_map<float, std::vector <std::string>>>> hashTable)
{

	ofstream hashTableFile;

	try
	{
		hashTableFile.open(fname, ios::out | ios::app);
		if (hashTableFile.fail())
		{
			cout << "\nError creating or reading the clusters file.";
			return -1;
		}
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

	hashTableFile << "Cluster_" << to_string(numberOfClusters) << endl;

	for (auto it1 = hashTable.begin(); it1 != hashTable.end(); ++it1)
		for (auto it2 = hashTable[it1->first].begin(); it2 != hashTable[it1->first].end(); ++it2)
			for (auto it3 = hashTable[it1->first][it2->first].begin(); it3 != hashTable[it1->first][it2->first].end(); ++it3)
			{
				hashTableFile << it1->first << "\t" << it2->first << "\t" << it3->first;
				for (int it4 = 0; it4 < (int)it3->second.size(); it4++)
					hashTableFile << "\t" << it3->second[it4];
				hashTableFile << "\n";
			}

	hashTableFile.close();

}

//****************************************************************************
/* This process calculates the translational degree of the points based on the basis point */

TranslationParameter CalculateGeoTranslation(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, int binSize)
{
	TranslationParameter res;
	//define Ty1,Tz1,Tx1 to store roation degrees
	float Ty1 = 0, Tz1 = 0, Tx1 = 0;
	float Ry1 = 0, Rz1 = 0, Rx1 = 0;
	float degree = 0;
	float sint = 0, cost = 0;
	float tmpx = 0, tmpy = 0, tmpz = 0;

	float newx2 = x2 - x1, newy2 = y2 - y1, newz2 = z2 - z1;
	float newx3 = x3 - x1, newy3 = y3 - y1, newz3 = z3 - z1;

	if (
		newx2 * newz3 == newx3 * newz2 &&
		newx2 * newy3 == newx3 * newy2 &&
		newy2 * newz3 == newy3 * newz2
		)
	{
		cout << "Error: Three points are collinear." << endl;
		return res;
	}

	// calculate the degree for Y axis
	//cout << "Y axis" << endl;

	if (newx2 == 0)
	{
		if (newz2 == 0)
			degree = 0;
		else if (newz2 > 0)
			degree = (float)PI / 2;
		else
			degree = (float)(-1 * PI / 2);
	}
	else
	{
		degree = atan2(newz2, newx2);
		/*if (newx2 < 0)
		{
		degree = degree + PI;
		}*/

	}

	Ry1 = degree;
	Ty1 = (float)(degree * 180 / PI);
	//cout << Ty1 << endl;

	if (Ry1 != 0)
	{
		sint = sin(Ry1);
		cost = cos(Ry1);
		//cout << "sint : " << sint << " , cost : " << cost << endl;
		tmpx = newx2*cost + newz2*sint;
		tmpz = -1 * newx2 *sint + newz2 *cost;
		newx2 = (float)((int)(tmpx * 10000 + .5) / 10000.0);
		newz2 = (float)((int)(tmpz * 10000 + .5) / 10000.0);
		//newx2 = tmpx;
		//newz2 = tmpz;
		//cout << "After rotation newx2:" << newx2 << " ,newy2:" << newy2 << " ,newz2:" << newz2 << endl;

		tmpx = newx3 * cost + newz3 * sint;
		tmpz = -1 * newx3 *sint + newz3 *cost;
		newx3 = (float)((int)(tmpx * 10000 + .5) / 10000.0);
		newz3 = (float)((int)(tmpz * 10000 + .5) / 10000.0);
		//newx2 = tmpx;
		//newz2 = tmpz;
		//cout << "After rotation newx3:" << newx3 << " ,newy3:" << newy3 << " ,newz3:" << newz3 << endl;

	}


	/// calculate the degree for Z axis
	//cout << "Z axis" << endl;

	if (newx2 == 0)
	{
		if (newy2 == 0)
			degree = 0;
		else if (newy2 > 0)
			degree = (float)(PI / 2);
		else
			degree = (float)(-1 * PI / 2);
	}
	else
	{
		degree = atan2(newy2, newx2);
		/*if (newx2 < 0)
		{
		degree = degree + PI;
		}*/

	}

	Rz1 = degree;
	Tz1 = (float)(degree * 180 / PI);
	//cout << "degree:" << Tz1 << endl;

	if (Rz1 != 0)
	{
		sint = sin(-1 * Rz1);
		cost = cos(-1 * Rz1);
		//cout << "sint : " << sint << " , cost : " << cost << endl;
		tmpx = newx2*cost + -1 * newy2*sint;
		tmpz = newx2 *sint + newy2 *cost;
		newx2 = (float)((int)(tmpx * 10000 + .5) / 10000.0);
		newy2 = (float)((int)(tmpz * 10000 + .5) / 10000.0);
		//newx2 = tmpx;
		//newz2 = tmpz;
		//cout << "After rotation newx2:" << newx2 << " ,newy2:" << newy2 << " ,newz2:" << newz2 << endl;

		tmpx = newx3 * cost + -1 * newy3 * sint;
		tmpz = newx3 *sint + newy3 * cost;
		newx3 = (float)((int)(tmpx * 10000 + .5) / 10000.0);
		newy3 = (float)((int)(tmpz * 10000 + .5) / 10000.0);
		//newx2 = tmpx;
		//newz2 = tmpz;
		//cout << "After rotation newx3:" << newx3 << " ,newy3:" << newy3 << " ,newz3:" << newz3 << endl;

	}

	// calculate the degree for X axis
	//cout << "X axis" << endl;

	if (newy3 == 0)
	{
		if (newz3 == 0)
			degree = 0;
		else if (newz3 > 0)
			degree = PI / 2;
		else
			degree = -1 * PI / 2;
	}
	else
	{
		degree = atan2(newz3, newy3);
		/*if (newy3 < 0)
		{
		degree = degree + PI;
		}*/

	}

	Rx1 = degree;
	Tx1 = degree * 180 / PI;
	//cout << "degree:" << Tx1 << endl;

	if (Rx1 != 0)
	{
		sint = sin(-1 * Rx1);
		cost = cos(-1 * Rx1);
		//cout << "sint : " << sint << " , cost : " << cost << endl;

		tmpx = newy3 * cost + -1 * newz3 * sint;
		tmpz = newy3 *sint + newz3 *cost;
		newy3 = ((int)(tmpx * 10000 + .5) / 10000.0);
		newz3 = ((int)(tmpz * 10000 + .5) / 10000.0);
		//newx2 = tmpx;
		//newz2 = tmpz;
		//cout << "After rotation newx3:" << newx3 << " ,newy3:" << newy3 << " ,newz3:" << newz3 << endl;

	}

	res.Transfer.x = x1;
	res.Transfer.y = y1;
	res.Transfer.z = z1;

	res.Ry = Ry1;
	res.Rz = Rz1;
	res.Rx = Rx1;

	/*res.p2.x = roundf(newx2 / binSize);
	res.p2.y = roundf(newy2 / binSize);
	res.p2.z = roundf(newz2 / binSize);

	res.p3.x = roundf(newx3 / binSize);
	res.p3.y = roundf(newy3 / binSize);
	res.p3.z = roundf(newz3 / binSize);*/

	res.p2.x = (newx2 <= 0 && newx2 > -1) ? 0 : roundf(newx2 / binSize);
	res.p2.y = (newy2 <= 0 && newy2 > -1) ? 0 : roundf(newy2 / binSize);
	res.p2.z = (newz2 <= 0 && newz2 > -1) ? 0 : roundf(newz2 / binSize);

	res.p3.x = (newx3 <= 0 && newx3 > -1) ? 0 : roundf(newx3 / binSize);
	res.p3.y = (newy3 <= 0 && newy3 > -1) ? 0 : roundf(newy3 / binSize);
	res.p3.z = (newz3 <= 0 && newz3 > -1) ? 0 : roundf(newz3 / binSize);

	return res;

}

//****************************************************************************
/* This process calculates the new position for a point according to the Reference Set (RS) */
Point CalculateNewPoint(TranslationParameter TP, Point P, int binSize)
{
	Point res = P;
	float sint = 0, cost = 0;
	float tmpx = 0, tmpy = 0, tmpz = 0;

	res.x = res.x - TP.Transfer.x;
	res.y = res.y - TP.Transfer.y;
	res.z = res.z - TP.Transfer.z;
	//cout << "After Transfer X:" << res.x << " ,Y:" << res.y << " ,Z:" << res.z << endl;
	sint = sin(TP.Ry);
	cost = cos(TP.Ry);
	tmpx = res.x*cost + res.z*sint;
	tmpz = -1 * res.x *sint + res.z *cost;
	res.x = ((int)(tmpx * 10000 + .5) / 10000.0);
	res.z = ((int)(tmpz * 10000 + .5) / 10000.0);
	//cout << "After Ry X:" << res.x << " ,Y:" << res.y << " ,Z:" << res.z << endl;
	sint = sin(-1 * TP.Rz);
	cost = cos(-1 * TP.Rz);
	tmpx = res.x*cost + -1 * res.y*sint;
	tmpz = res.x *sint + res.y *cost;
	res.x = ((int)(tmpx * 10000 + .5) / 10000.0);
	res.y = ((int)(tmpz * 10000 + .5) / 10000.0);
	//cout << "After Rz X:" << res.x << " ,Y:" << res.y << " ,Z:" << res.z << endl;
	sint = sin(-1 * TP.Rx);
	cost = cos(-1 * TP.Rx);

	tmpx = res.y * cost + -1 * res.z * sint;
	tmpz = res.y *sint + res.z *cost;
	res.y = ((int)(tmpx * 10000 + .5) / 10000.0);
	res.z = ((int)(tmpz * 10000 + .5) / 10000.0);


	/*res.x = roundf(res.x);
	res.y = roundf(res.y);
	res.z = roundf(res.z);*/

	res.x = (res.x <= 0 && res.x > -1) ? 0 : roundf(res.x / binSize);
	res.y = (res.y <= 0 && res.y > -1) ? 0 : roundf(res.y / binSize);
	res.z = (res.z <= 0 && res.z > -1) ? 0 : roundf(res.z / binSize);
	//cout << "X:" << res.x << " ,Y:" << res.y << " ,Z:" << res.z << endl;
	return res;
}


//****************************************************************************
/* This process adds a new point to a hashtable  */

void AddToHashTable(unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>>& pointHashTable, Point p, int rsX, int rsY, int rsZ)
{

	if (pointHashTable.find(p.x) != pointHashTable.end() && pointHashTable[p.x].find(p.y) != pointHashTable[p.x].end() && pointHashTable[p.x][p.y].find(p.z) != pointHashTable[p.x][p.y].end())	// found
	{
		pointHashTable[p.x][p.y][p.z].insert(pointHashTable[p.x][p.y][p.z].end(), "(" + to_string(rsX) + "," + to_string(rsY) + "," + to_string(rsZ) + ")");
	}
	else // not found
	{
		pointHashTable[p.x][p.y][p.z] = { "(" + to_string(rsX) + "," + to_string(rsY) + "," + to_string(rsZ) + ")" };
	}
}

//****************************************************************************
/* This process compares two hash tables belong to two fragments */
bool CompareTwoHashTables(unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> centroidHashTable, unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> nonCentroidHashTable, int expectedMatchedPoints)
{

	unordered_map<string, int> matchesHashTable;	//stores each reference set model has how many matched points with the new fragment model
													//tries to find the x, y, z coordinates saved in nonCentroidHashTable in centroidHashTable

	for (auto it1 = nonCentroidHashTable.begin(); it1 != nonCentroidHashTable.end(); ++it1)
	{
		if (centroidHashTable.find(it1->first) != centroidHashTable.end())  //we have found x
		{
			for (auto it2 = nonCentroidHashTable[it1->first].begin(); it2 != nonCentroidHashTable[it1->first].end(); ++it2)
			{
				if (centroidHashTable[it1->first].find(it2->first) != centroidHashTable[it1->first].end())  //we have found y
				{
					for (auto it3 = nonCentroidHashTable[it1->first][it2->first].begin(); it3 != nonCentroidHashTable[it1->first][it2->first].end(); ++it3)
					{
						if (centroidHashTable[it1->first][it2->first].find(it3->first) != centroidHashTable[it1->first][it2->first].end())  //we have found z
						{
							for (int it4 = 0; it4 < it3->second.size(); it4++)
							{
								if (matchesHashTable.find(it3->second[it4]) != matchesHashTable.end())
								{
									//cout << "Adding to hash table.";
									matchesHashTable[it3->second[it4]]++;
								}
								else
								{
									//cout << "Initializing the hash table.";
									matchesHashTable[it3->second[it4]] = 1;
								}
							}
						}
						else  //not found z
						{
							return false;
						}
					}//for
				}
				else  //not found y
				{
					return false;
				}
			}//for it2
		}
		else  //not found x
		{
			//cout << "\nNot found x.";
			return false;
		}
	}//for it1

	int maxMatches = 0;

	for (auto it = matchesHashTable.begin(); it != matchesHashTable.end(); ++it)
	{
		//cout << matchesHashTable[it->first] << endl;
		if (matchesHashTable[it->first] > maxMatches)
			maxMatches = matchesHashTable[it->first];
	}
	//cout << "\nMax matches is: " << maxMatches;
	if (maxMatches >= expectedMatchedPoints)
		return true;
	else
		return false;

}

