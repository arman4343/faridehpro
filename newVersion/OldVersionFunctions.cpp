#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <map>
#include <ctime>
#include "My_Functions.h"


#define PI 3.14159265

using namespace std;

#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4503) 
//****************************************************************************
/*
This process creates the descriptor vector for PIFACE nonrepresentative interfaces
*/

int CreateProteinDescriptors_v1(int fragmentLength, int overlappingResidues, int binSize, int expectedMatchedPoints, int numberOfClusters)
{

	int *proteinDescriptorVector;
	proteinDescriptorVector = new int[numberOfClusters];

	string line, proteinFileName;
	int incrementValue = fragmentLength - overlappingResidues;

	ifstream proteinList, clustersFile;

	//Read the list of proteins
	try
	{
		proteinList.open("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Piface\\PifaceMixedList.txt"); //list of target piface interfaces

		if (proteinList.fail())
		{
			cout << "\nError reading the proteins list file.";
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


	//Read the clusters file
	try
	{
		clustersFile.open("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\Clusters.txt");
		if (clustersFile.fail())
		{
			cout << "\nError reading the cluster hash tables file.";
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
	getline(clustersFile, line);	//To pass the line having cluster number




	while (getline(proteinList, proteinFileName))	//each line has an interface name
	{

		//proteinFileName = proteinFileName.substr(0, proteinFileName.size() - 1);
		cout << proteinFileName << endl;

		ifstream proteinDescriptorVectorFile("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\Piface_Interfaces_Descriptor_Vectors\\" + proteinFileName + "_DescriptorVector.txt");
		if (!proteinDescriptorVectorFile)	//if the descriptor vector file of this protein doesn't exist
		{

			ofstream proteinDescriptorVectorFile("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\Piface_Interfaces_Descriptor_Vectors\\" + proteinFileName + "_DescriptorVector.txt");

			ifstream proteinFile(("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Piface\\Interface Structures\\" + proteinFileName + ".pdb").c_str());	//read the protein surface file
			if (proteinFile && proteinDescriptorVectorFile)
			{

				int counter = 0;
				while (getline(proteinFile, line))	//count the number of lines in the interface file
				{
					counter++;
				}
				proteinFile.close();
				//cout << counter << "\n";

				string *lines;
				lines = new string[counter];
				ifstream proteinFile(("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Piface\\Interface Structures\\" + proteinFileName + ".pdb").c_str());
				int lineCounter = 0;
				while (getline(proteinFile, line))	//read the lines of the interface file
				{
					lines[lineCounter++] = line;
				}
				proteinFile.close();

				for (int i = 0; i < numberOfClusters; i++) //initialize the interface descriptor vector 
				{
					proteinDescriptorVector[i] = 0;
				}

				unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> fragmentHashTable;

				for (int lineIndex1 = 0; lineIndex1 <= counter - (fragmentLength); lineIndex1 = lineIndex1 + incrementValue)	//extracts continuous fragments from interfaces
				{
					//cout << "\nincrement value: " << incrementValue;
					//cout << "\nlineIndex: " << lineIndex1;
					int previousResidueNumber = atoi(lines[lineIndex1].substr(22, 5).c_str());
					char previousResidueChainId = lines[lineIndex1][21];
					bool continuousResidues = true;

					for (int lineIndex2 = (lineIndex1 + 1); lineIndex2 < (lineIndex1 + fragmentLength); lineIndex2++)
					{

						int nextResidueNumber = atoi(lines[lineIndex2].substr(22, 5).c_str());
						char nextResidueChainId = lines[lineIndex2][21];

						if (previousResidueNumber < nextResidueNumber - 1 || previousResidueChainId != nextResidueChainId)
						{
							//cout << "\nnon continuous residues\t" << lineIndex2;
							continuousResidues = false;
							lineIndex1 = lineIndex2 - incrementValue;/////////////
							break;
						}
						previousResidueNumber = nextResidueNumber;
						previousResidueChainId = nextResidueChainId;
					}


					if (continuousResidues == true)	//We made sure that it is a continuous fragment. Now we want to find which fragment in our fragment library is similar to this fragment.
					{

						/*for (int tmp = lineIndex1; tmp < (lineIndex1 + fragmentLength); tmp++)
						cout << "\n" << lines[tmp];*/

						unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> clusterHashTable;
						bool matched = false;
						int scannedClusters = 1;

						while (!matched && scannedClusters <= numberOfClusters) //till we find which cluster hash table is similar to the fragment hashtable or none of them are similar
						{
							//cout << "\n" << scannedClusters;

							LoadHashTableFromFile(clustersFile, scannedClusters, clusterHashTable);


							for (int c = lineIndex1; c < (lineIndex1 + fragmentLength - 2); c++)	// we select three consecutive points as Reference Set (RS)
							{
								for (int i = c; i < c + 3; i++)
								{
									for (int j = c; j < c + 3; j++)
									{
										for (int k = c; k < c + 3; k++)
											if (i != j && i != k && j != k)	//three different points to create RS
											{

												fragmentHashTable.clear();

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
												//cout << "GeoTranslation started.\n";
												selectedRS = CalculateGeoTranslation(x1, y1, z1, x2, y2, z2, x3, y3, z3, binSize);
												//cout << "GeoTranslation success.\n";

												AddToHashTable(fragmentHashTable, selectedRS.p2, i + 1, j + 1, k + 1);
												AddToHashTable(fragmentHashTable, selectedRS.p3, i + 1, j + 1, k + 1);

												for (int f = lineIndex1; f < fragmentLength + lineIndex1; f++)
												{
													if (f != i && f != j && f != k)
													{
														Point p;
														p.x = (float)atof(lines[f].substr(30, 8).c_str());
														p.y = (float)atof(lines[f].substr(38, 8).c_str());
														p.z = (float)atof(lines[f].substr(46, 8).c_str());
														p = CalculateNewPoint(selectedRS, p, binSize);
														AddToHashTable(fragmentHashTable, p, i + 1, j + 1, k + 1);
													}
												}


												//compares two hash tables
												matched = CompareTwoHashTables(clusterHashTable, fragmentHashTable, expectedMatchedPoints);
												if (matched) // if the cluster hash table and the new one are similar
													break;
											}
										if (matched) break;
									}// for j
									if (matched) break;
								}// for i
								if (matched) break;
							}//for c

							if (matched)
								proteinDescriptorVector[scannedClusters - 1]++;

							scannedClusters++;
						}
						clustersFile.clear();
						clustersFile.seekg(0, ios::beg);
						getline(clustersFile, line);	//To pass the line having cluster number
					}
				}

				//proteinDescriptorVectorFile << interfaceFileName << "\t";
				for (int i = 0; i < numberOfClusters; i++) //Write the protein descriptor vector into file
				{
					proteinDescriptorVectorFile << proteinDescriptorVector[i] << "\t";
				}
				proteinDescriptorVectorFile << "\n";
				proteinDescriptorVectorFile.close();

			}
			else
				cout << "\nError reading the interface file or creating the interface descriptor vector file " + proteinFileName;
		}
	}
	clustersFile.close();
	proteinList.close();
}


//****************************************************************************
/*
This process creates descriptor vectors for protein in PRISM pair_list from which the surfaces are extracted.
It creates the descriptor vector as an array which is of course sparse.
*/
int CreateProteinDescriptors(string jobName, int fragmentLength, int overlappingResidues, int binSize, int expectedMatchedPoints, int numberOfClusters)
{

	int *proteinDescriptorVector;
	proteinDescriptorVector = new int[numberOfClusters];

	string line, proteinName;
	int incrementValue = fragmentLength - overlappingResidues;
	vector <string> proteinList;

	ifstream pairList, clustersFile;

	//Read the pair_list file
	try
	{
		pairList.open("D:\\PhD\\prism\\prism_standalone\\pair_list"); //list of template proteins

		if (pairList.fail())
		{
			cout << "\nError reading the pair list file.";
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


	//Read the clusters file
	try
	{
		clustersFile.open("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\Clusters.txt");
		if (clustersFile.fail())
		{
			cout << "\nError reading the cluster hash tables file.";
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
	getline(clustersFile, line);	//To pass the line having cluster number



	while (getline(pairList, line))	//each line has two protein names
	{
		int spac = line.find(" ", 0);

		proteinName = line.substr(0, spac);
		if (!any_of(proteinList.begin(), proteinList.end(), bind2nd(std::equal_to<std::string>(), proteinName)))	// not found
			proteinList.push_back(proteinName);

		proteinName = line.substr(spac + 1, line.length());
		if (!any_of(proteinList.begin(), proteinList.end(), bind2nd(std::equal_to<std::string>(), proteinName)))	// not found
			proteinList.push_back(proteinName);
	}


	for (int i = 0; i < (int)proteinList.size(); ++i)
	{

		proteinName = proteinList[i];

		ifstream proteinFile(("D:\\PhD\\prism\\prism_standalone\\jobs\\" + jobName + "\\surfaceExtract\\" + proteinName + ".asa.pdb").c_str());	//read the protein surface file
		if (proteinFile)
		{
			ifstream proteinDescriptorVectorFile("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\Protein_Descriptor_Vectors\\" + proteinName + "DescriptorVectorNonCompact.txt");
			if (!proteinDescriptorVectorFile)	//if the descriptor vector file of this protein doesn't exist
			{
				ofstream proteinDescriptorVectorFile("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\Protein_Descriptor_Vectors\\" + proteinName + "DescriptorVectorNonCompact.txt");

				if (proteinDescriptorVectorFile)
				{

					vector <string> surfaceResidueCoordinates;

					int lineCounter = 0;
					while (getline(proteinFile, line))	//read the lines of the protein file
					{
						if (line.find("ATOM") == 0)
							surfaceResidueCoordinates.push_back(line);
					}
					proteinFile.close();

					for (int i = 0; i < numberOfClusters; i++) //initialize the protein descriptor vector 
						proteinDescriptorVector[i] = 0;

					unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> fragmentHashTable;

					for (int lineIndex1 = 0; lineIndex1 <= (int)surfaceResidueCoordinates.size() - (fragmentLength); lineIndex1 = lineIndex1 + incrementValue)	//extracts continuous fragments from interfaces
					{

						int previousResidueNumber = atoi(surfaceResidueCoordinates[lineIndex1].substr(22, 5).c_str());
						bool continuousResidues = true;

						for (int lineIndex2 = (lineIndex1 + 1); lineIndex2 < (lineIndex1 + fragmentLength); lineIndex2++)
						{

							int nextResidueNumber = atoi(surfaceResidueCoordinates[lineIndex2].substr(22, 5).c_str());

							if (previousResidueNumber < nextResidueNumber - 1)
							{

								continuousResidues = false;
								lineIndex1 = lineIndex2 - incrementValue;
								break;
							}
							previousResidueNumber = atoi(surfaceResidueCoordinates[lineIndex2].substr(22, 5).c_str());
						}


						if (continuousResidues == true)	//We made sure that it is a continuous fragment. Now we want to find which fragment in our fragment library is similar to this fragment.
						{

							unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> clusterHashTable;
							bool matched = false;
							int scannedClusters = 1;

							while (!matched && scannedClusters <= numberOfClusters) //till we find which cluster hash table is similar to the fragment hashtable or none of them are similar
							{

								LoadHashTableFromFile(clustersFile, scannedClusters, clusterHashTable);

								for (int c = lineIndex1; c < (lineIndex1 + fragmentLength - 2); c++)	// we select three consecutive points as Reference Set (RS)
								{
									for (int i = c; i < c + 3; i++)
									{
										for (int j = c; j < c + 3; j++)
										{
											for (int k = c; k < c + 3; k++)
												if (i != j && i != k && j != k)	//three different points to create RS
												{

													fragmentHashTable.clear();

													float x1 = (float)atof(surfaceResidueCoordinates[i].substr(30, 8).c_str());
													float y1 = (float)atof(surfaceResidueCoordinates[i].substr(38, 8).c_str());
													float z1 = (float)atof(surfaceResidueCoordinates[i].substr(46, 8).c_str());

													float x2 = (float)atof(surfaceResidueCoordinates[j].substr(30, 8).c_str());
													float y2 = (float)atof(surfaceResidueCoordinates[j].substr(38, 8).c_str());
													float z2 = (float)atof(surfaceResidueCoordinates[j].substr(46, 8).c_str());

													float x3 = (float)atof(surfaceResidueCoordinates[k].substr(30, 8).c_str());
													float y3 = (float)atof(surfaceResidueCoordinates[k].substr(38, 8).c_str());
													float z3 = (float)atof(surfaceResidueCoordinates[k].substr(46, 8).c_str());

													TranslationParameter selectedRS;

													selectedRS = CalculateGeoTranslation(x1, y1, z1, x2, y2, z2, x3, y3, z3, binSize);

													AddToHashTable(fragmentHashTable, selectedRS.p2, i + 1, j + 1, k + 1);
													AddToHashTable(fragmentHashTable, selectedRS.p3, i + 1, j + 1, k + 1);

													for (int f = lineIndex1; f < fragmentLength + lineIndex1; f++)
													{
														if (f != i && f != j && f != k)
														{
															Point p;
															p.x = (float)atof(surfaceResidueCoordinates[f].substr(30, 8).c_str());
															p.y = (float)atof(surfaceResidueCoordinates[f].substr(38, 8).c_str());
															p.z = (float)atof(surfaceResidueCoordinates[f].substr(46, 8).c_str());
															p = CalculateNewPoint(selectedRS, p, binSize);
															AddToHashTable(fragmentHashTable, p, i + 1, j + 1, k + 1);
														}
													}


													/*cout << "\n-------------- Printing Fragment Hash Table ---------------" << endl << endl;
													for (auto it1 = fragmentHashTable.begin(); it1 != fragmentHashTable.end(); ++it1)
													{
													cout << it1->first << ":";
													for (auto it2 = fragmentHashTable[it1->first].begin(); it2 != fragmentHashTable[it1->first].end(); ++it2)
													{
													cout << "   " << it2->first << ":";
													for (auto it3 = fragmentHashTable[it1->first][it2->first].begin(); it3 != fragmentHashTable[it1->first][it2->first].end(); ++it3)
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

													getchar();*/


													//compares two hash tables
													matched = CompareTwoHashTables(clusterHashTable, fragmentHashTable, expectedMatchedPoints);
													if (matched) // if the cluster hash table and the new one are similar
													{
														//cout << "\nMatching output is: " << matched;
														break;
													}
												}
											if (matched) break;
										}// for j
										if (matched) break;
									}// for i
									if (matched) break;
								}//for c

								if (matched)
									proteinDescriptorVector[scannedClusters - 1]++;

								scannedClusters++;
							}
							clustersFile.clear();
							clustersFile.seekg(0, ios::beg);
							getline(clustersFile, line);	//To pass the line having cluster number
						}
					}

					for (int i = 0; i < numberOfClusters; i++) //Write the protein descriptor vector into file
					{
						proteinDescriptorVectorFile << proteinDescriptorVector[i] << "\t";
					}
					proteinDescriptorVectorFile << "\n";
					proteinDescriptorVectorFile.close();

				}
				else
					cout << "\nError creating the protein surface descriptor vector file " + proteinName;
			}

			proteinFile.close();
			proteinDescriptorVectorFile.close();
		}
		else
			cout << "\nError reading the protein surface file " + proteinName;
	}
}


//****************************************************************************************
/*
This process compares the interface descriptor vectors with protein surface descriptor vectors.
The vectors are sparse.
*/
void CompareProteinWithInterfaces(int fragmentLength, int overlappingResidues, int binSize, int expectedMatchedPoints, int numberOfClusters)
{

	int *proteinDescriptorVector;
	proteinDescriptorVector = new int[numberOfClusters];
	int *interfaceDescriptorVector;
	interfaceDescriptorVector = new int[numberOfClusters];
	int *descriptorVectorDifferences;
	descriptorVectorDifferences = new int[numberOfClusters];
	string line, proteinFileName;

	ifstream proteinList("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Piface\\PifaceMixedList.txt"); //list of template protein surfaces
	if (proteinList)
	{
		ofstream interfacesSimilarToProteins("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\PifaceInterfacesSimilarToProteins.txt");
		if (interfacesSimilarToProteins)
		{

			while (getline(proteinList, proteinFileName))	//each line has a protein name
			{

				interfacesSimilarToProteins << "Interfaces similar to " << proteinFileName << ":\n";

				ifstream proteinDescriptorVectorFile("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\Piface_Interfaces_Descriptor_Vectors\\" + proteinFileName + "_DescriptorVector.txt");
				if (proteinDescriptorVectorFile)
				{
					getline(proteinDescriptorVectorFile, line);
					istringstream iss1(line);
					vector <string> lineTokens1{ istream_iterator<string>{iss1}, istream_iterator<string>{} };

					// Read the line and put it in an array
					for (int i = 0; i < numberOfClusters; i++)
						proteinDescriptorVector[i] = stoi(lineTokens1[i]);

					//for (int i = 0; i < numberOfClusters; i++)
					//cout << proteinDescriptorVector[i] << "\t";

					bool fragmentByFragmentSimilar = true;

					ifstream interfaceDescriptorVectorLibrary("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\PifaceRepresentativeInterfacesDescriptorVectorLibray.txt"); //list of all PRISM interfaces
					if (interfaceDescriptorVectorLibrary)
					{
						while (getline(interfaceDescriptorVectorLibrary, line)) //Read each line (interface descriptor vector)
						{
							istringstream iss2(line);
							vector <string> lineTokens2{ istream_iterator<string>{iss2}, istream_iterator<string>{} };
							int countSimilarFrequencies = 0;

							for (int i = 0; i < numberOfClusters; i++)
							{
								interfaceDescriptorVector[i] = stoi(lineTokens2[i + 1]);
								//cout << interfaceDescriptorVector[i] << "\t";
								descriptorVectorDifferences[i] = abs(proteinDescriptorVector[i] - interfaceDescriptorVector[i]);

								if (descriptorVectorDifferences[i] <= 3) //count how many elements are the same
									countSimilarFrequencies++;

								if (fragmentByFragmentSimilar && descriptorVectorDifferences[i] > 3)
									fragmentByFragmentSimilar = false;

							}
							//cout << endl;

							if (fragmentByFragmentSimilar || (countSimilarFrequencies / numberOfClusters) >= 0.7)
								interfacesSimilarToProteins << lineTokens2[0] << "\n";

						}
					}
					else
						cout << "\nError reading the interface descriptor vector library file.";
				}
				else
					cout << "\nError reading the protein descriptor vector file " + proteinFileName;
			}
			interfacesSimilarToProteins.close();
		}
		else
			cout << "\nError creating the output file.";
	}
	else
		cout << "\nError reading the protein list file.";

}


//****************************************************************************
/* This process creates the descriptor vectors for the template interface partners of PRISM */

int CreateInterfaceDescriptors(int fragmentLength, int overlappingResidues, int binSize, int expectedMatchedPoints, int numberOfClusters)
{

	int *interfaceDescriptorVector;
	interfaceDescriptorVector = new int[numberOfClusters];

	string line, interfaceName = " ";
	int incrementValue = fragmentLength - overlappingResidues;
	int lineInFragment = -1;
	bool continuousResidues;
	string *fragment;
	fragment = new string[fragmentLength];

	ifstream interfacesFile, clustersFile;
	ofstream interfaceVectorLibraryFile;


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
		cout << "\nOther exception thrown." << endl;
		return -1;
	}


	try
	{
		interfaceVectorLibraryFile.open("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\InterfaceDescriptorVectorLibray.txt");

		if (interfaceVectorLibraryFile.fail())
		{
			cout << "\nError creating the interface vector library file.";
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
		cout << "\nOther exception thrown" << endl;
		return -1;
	}


	try
	{
		clustersFile.open("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\Clusters.txt");

		if (clustersFile.fail())
		{
			cout << "\nError reading the cluster hash tables file.";
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
	getline(clustersFile, line);	//To pass the line having cluster number



	while (getline(interfacesFile, line))	//each line has an interface name or residue information
	{

		//check if it is a line showing an interface name or not
		if (line.find("Interface Name") == 0)
		{

			if (interfaceName != " ")
			{
				interfaceVectorLibraryFile << interfaceName << "\t";
				for (int i = 0; i < numberOfClusters; i++) //Write the interface descriptor vector into file
				{
					interfaceVectorLibraryFile << interfaceDescriptorVector[i] << "\t";
				}
				interfaceVectorLibraryFile << "\n";
			}

			interfaceName = line.substr(16);
			lineInFragment = 0;
			continuousResidues = true;
			for (int i = 0; i < numberOfClusters; i++) //initialize the interface descriptor vector 
			{
				interfaceDescriptorVector[i] = 0;
			}
		}
		else
		{

			fragment[lineInFragment] = line;

			if (lineInFragment == 0)
				lineInFragment++;
			else if (lineInFragment > 0)
			{

				int previousResidueNumber = atoi(fragment[lineInFragment - 1].substr(22, 5).c_str());
				int nextResidueNumber = atoi(fragment[lineInFragment].substr(22, 5).c_str());

				if (previousResidueNumber < nextResidueNumber - 1)
				{
					continuousResidues = false;
					fragment[0] = fragment[lineInFragment];
					lineInFragment = 1;
				}
				else
					lineInFragment++;
			}

			if (lineInFragment == fragmentLength)	//if it is a continuous fragment
			{

				//for (int h = 0; h < fragmentLength; h++)
				//cout << endl << fragment[h];

				unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> fragmentHashTable;

				//fragmentHashTable.clear();
				unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> clusterHashTable;
				bool matched = false;
				int scannedClusters = 1;

				while (!matched && scannedClusters <= numberOfClusters) //till we find which cluster hash table is similar to the fragment hashtable or none of them are similar
				{
					//cout << "\n" << scannedClusters;

					//clusterHashTable.clear();
					LoadHashTableFromFile(clustersFile, scannedClusters, clusterHashTable);

					//
					////
					///////

					for (int c = 0; c < (fragmentLength - 2); c++)	// we select three consecutive points as Reference Set (RS)
					{
						for (int i = c; i < c + 3; i++)
						{
							for (int j = c; j < c + 3; j++)
							{
								for (int k = c; k < c + 3; k++)
								{
									if (i != j && i != k && j != k)	//three different points to create RS
									{

										fragmentHashTable.clear();

										float x1 = (float)atof(fragment[i].substr(30, 8).c_str());
										float y1 = (float)atof(fragment[i].substr(38, 8).c_str());
										float z1 = (float)atof(fragment[i].substr(46, 8).c_str());

										float x2 = (float)atof(fragment[j].substr(30, 8).c_str());
										float y2 = (float)atof(fragment[j].substr(38, 8).c_str());
										float z2 = (float)atof(fragment[j].substr(46, 8).c_str());

										float x3 = (float)atof(fragment[k].substr(30, 8).c_str());
										float y3 = (float)atof(fragment[k].substr(38, 8).c_str());
										float z3 = (float)atof(fragment[k].substr(46, 8).c_str());

										TranslationParameter selectedRS;
										//cout << "GeoTranslation started.\n";
										selectedRS = CalculateGeoTranslation(x1, y1, z1, x2, y2, z2, x3, y3, z3, binSize);
										//cout << "GeoTranslation success.\n";

										AddToHashTable(fragmentHashTable, selectedRS.p2, i + 1, j + 1, k + 1);
										AddToHashTable(fragmentHashTable, selectedRS.p3, i + 1, j + 1, k + 1);

										for (int f = 0; f < fragmentLength; f++)
										{
											if (f != i && f != j && f != k)
											{
												Point p;
												p.x = (float)atof(fragment[f].substr(30, 8).c_str());
												p.y = (float)atof(fragment[f].substr(38, 8).c_str());
												p.z = (float)atof(fragment[f].substr(46, 8).c_str());
												p = CalculateNewPoint(selectedRS, p, binSize);
												AddToHashTable(fragmentHashTable, p, i + 1, j + 1, k + 1);
											}
										}




										/*cout << "\n-------------- Printing Fragment Hash Table ---------------" << endl << endl;
										for (auto it1 = fragmentHashTable.begin(); it1 != fragmentHashTable.end(); ++it1)
										{
										cout << it1->first << ":";
										for (auto it2 = fragmentHashTable[it1->first].begin(); it2 != fragmentHashTable[it1->first].end(); ++it2)
										{
										cout << "   " << it2->first << ":";
										for (auto it3 = fragmentHashTable[it1->first][it2->first].begin(); it3 != fragmentHashTable[it1->first][it2->first].end(); ++it3)
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

										getchar();*/





										//compares two hash tables
										matched = CompareTwoHashTables(clusterHashTable, fragmentHashTable, expectedMatchedPoints);
										if (matched) // if the cluster hash table and the new one are similar
										{
											//cout << "\nMatching output is: " << matched;
											break;
										}
									}
								}// for k
								if (matched) break;
							}// for j
							if (matched) break;
						}// for i
						if (matched) break;
					}//for c

					if (matched)
					{
						interfaceDescriptorVector[scannedClusters - 1]++;
						//cout << "\nMatched";
					}

					scannedClusters++;
				}

				clustersFile.clear();
				clustersFile.seekg(0, ios::beg);
				getline(clustersFile, line);	//To pass the line having cluster number

												/*lineInFragment = fragmentLength - 1;
												for (int i = 1; i < fragmentLength; i++)
												fragment[i - 1] = fragment[i];*/
				lineInFragment = overlappingResidues;
				for (int i = 0; i < overlappingResidues; i++)
					fragment[i] = fragment[fragmentLength - overlappingResidues + i];

			}
		}
	}

	interfaceVectorLibraryFile << interfaceName << "\t";
	for (int i = 0; i < numberOfClusters; i++) //Write the last interface descriptor vector into the file
	{
		interfaceVectorLibraryFile << interfaceDescriptorVector[i] << "\t";
	}
	interfaceVectorLibraryFile << "\n";

	interfaceVectorLibraryFile.close();
	clustersFile.close();
	interfacesFile.close();
}


//****************************************************************************
/*
This process creates the descriptor vectors for the template representative interfaces of PIFACE
PIFACE interfaces have two sides in the same file
*/

int CreateInterfaceDescriptors_v1(int fragmentLength, int overlappingResidues, int binSize, int expectedMatchedPoints, int numberOfClusters)
{

	int *interfaceDescriptorVector;
	interfaceDescriptorVector = new int[numberOfClusters];

	string line, interfaceFileName;
	int incrementValue = fragmentLength - overlappingResidues;
	//bool continuousResidues;

	ifstream interfacesList, clustersFile;
	ofstream interfaceDescriptorVectorLibraryFile;


	//Reads the interfaces list
	try
	{
		interfacesList.open("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Piface\\tempInterfaceRepresentativesList.txt"); //All PIFACE interface representatives

		if (interfacesList.fail())
		{
			cout << "\nError reading the interface list file.";
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


	//Creates the interface descriptor vectors file
	try
	{
		interfaceDescriptorVectorLibraryFile.open("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\PifaceRepresentativeInterfacesDescriptorVectorLibray.txt");

		if (interfaceDescriptorVectorLibraryFile.fail())
		{
			cout << "\nError creating the interface vector library file.";
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
		cout << "\nOther exception thrown";
		return -1;
	}


	//Reads the fragment clusters file
	try
	{
		clustersFile.open("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\Clusters.txt");

		if (clustersFile.fail())
		{
			cout << "\nError reading the cluster hash tables file.";
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
	getline(clustersFile, line); //To pass the line having cluster number



	while (getline(interfacesList, interfaceFileName)) //each line has an interface partner name
	{

		//interfaceFileName = interfaceFileName.substr(0, interfaceFileName.size() - 1);
		cout << endl << interfaceFileName;

		ifstream interfaceFile(("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Piface\\Interface Structures\\" + interfaceFileName + ".pdb").c_str());	//read the interface file

		if (interfaceFile)
		{

			int counter = 0;
			while (getline(interfaceFile, line)) //count the number of lines in the interface file
			{
				counter++;
			}
			interfaceFile.close();
			//cout << counter << "\n";

			string *interfaceLines;
			interfaceLines = new string[counter];

			ifstream interfaceFile(("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Piface\\Interface Structures\\" + interfaceFileName + ".pdb").c_str());	//read the interface file
			int lineCounter = 0;
			while (getline(interfaceFile, line)) //read the lines of the interface file
			{
				interfaceLines[lineCounter++] = line;
			}
			interfaceFile.close();

			for (int i = 0; i < numberOfClusters; i++) //initialize the interface descriptor vector 
			{
				interfaceDescriptorVector[i] = 0;
			}

			unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> fragmentHashTable;

			for (int lineIndex1 = 0; lineIndex1 <= counter - (fragmentLength); lineIndex1 = lineIndex1 + incrementValue)	//extracts continuous fragments from interfaces
			{
				//cout << "\nincrement value: " << incrementValue;
				//cout << "\nlineIndex: " << lineIndex1;
				int previousResidueNumber = atoi(interfaceLines[lineIndex1].substr(22, 5).c_str());
				char previousResidueChainId = interfaceLines[lineIndex1][21];
				bool continuousResidues = true;

				for (int lineIndex2 = (lineIndex1 + 1); lineIndex2 < (lineIndex1 + fragmentLength); lineIndex2++)
				{

					int nextResidueNumber = atoi(interfaceLines[lineIndex2].substr(22, 5).c_str());
					char nextResidueChainId = interfaceLines[lineIndex2][21];

					if (previousResidueNumber < nextResidueNumber - 1 || previousResidueChainId != nextResidueChainId)
					{
						//cout << "\nnon continuous residues\t" << lineIndex2;
						continuousResidues = false;
						lineIndex1 = lineIndex2 - incrementValue;/////////
						break;
					}
					previousResidueNumber = nextResidueNumber;
					previousResidueChainId = nextResidueChainId;
				}


				if (continuousResidues == true)	//We made sure that it is a continuous fragment. Now we want to find which fragment in our fragment library is similar to this fragment.
				{

					//for (int tmp = lineIndex1; tmp < (lineIndex1 + fragmentLength); tmp++)
					//cout << "\n" << interfaceLines[tmp];


					unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> clusterHashTable;
					bool matched = false;
					int scannedClusters = 1;

					while (!matched && scannedClusters <= numberOfClusters) //till we find which cluster hash table is similar to the fragment hashtable or none of them are similar
					{
						//cout << "\n" << scannedClusters;

						LoadHashTableFromFile(clustersFile, scannedClusters, clusterHashTable);


						for (int c = lineIndex1; c < (lineIndex1 + fragmentLength - 2); c++)	// we select three consecutive points as Reference Set (RS)
						{
							for (int i = c; i < c + 3; i++)
							{
								for (int j = c; j < c + 3; j++)
								{
									for (int k = c; k < c + 3; k++)
										if (i != j && i != k && j != k)	//three different points to create RS
										{

											fragmentHashTable.clear();

											float x1 = (float)atof(interfaceLines[i].substr(30, 8).c_str());
											float y1 = (float)atof(interfaceLines[i].substr(38, 8).c_str());
											float z1 = (float)atof(interfaceLines[i].substr(46, 8).c_str());

											float x2 = (float)atof(interfaceLines[j].substr(30, 8).c_str());
											float y2 = (float)atof(interfaceLines[j].substr(38, 8).c_str());
											float z2 = (float)atof(interfaceLines[j].substr(46, 8).c_str());

											float x3 = (float)atof(interfaceLines[k].substr(30, 8).c_str());
											float y3 = (float)atof(interfaceLines[k].substr(38, 8).c_str());
											float z3 = (float)atof(interfaceLines[k].substr(46, 8).c_str());

											TranslationParameter selectedRS;
											//cout << "GeoTranslation started.\n";
											selectedRS = CalculateGeoTranslation(x1, y1, z1, x2, y2, z2, x3, y3, z3, binSize);
											//cout << "GeoTranslation success.\n";

											AddToHashTable(fragmentHashTable, selectedRS.p2, i + 1, j + 1, k + 1);
											AddToHashTable(fragmentHashTable, selectedRS.p3, i + 1, j + 1, k + 1);

											for (int f = lineIndex1; f < fragmentLength + lineIndex1; f++)
											{
												if (f != i && f != j && f != k)
												{

													Point p;
													p.x = (float)atof(interfaceLines[f].substr(30, 8).c_str());
													p.y = (float)atof(interfaceLines[f].substr(38, 8).c_str());
													p.z = (float)atof(interfaceLines[f].substr(46, 8).c_str());
													p = CalculateNewPoint(selectedRS, p, binSize);
													AddToHashTable(fragmentHashTable, p, i + 1, j + 1, k + 1);
												}
											}


											matched = CompareTwoHashTables(clusterHashTable, fragmentHashTable, expectedMatchedPoints);	//compares two hash tables
											if (matched) //if the cluster hash table and the new one are similar
												break;
										}
									if (matched) break;
								}// for j
								if (matched) break;
							}// for i
							if (matched) break;
						}//for c

						if (matched)
							interfaceDescriptorVector[scannedClusters - 1]++;

						scannedClusters++;
					}

					clustersFile.clear();
					clustersFile.seekg(0, ios::beg);
					getline(clustersFile, line);	//To pass the line having cluster number

				}
			}

			interfaceDescriptorVectorLibraryFile << interfaceFileName << "\t";
			for (int i = 0; i < numberOfClusters; i++) //Write the last interface descriptor vector into the file
			{
				interfaceDescriptorVectorLibraryFile << interfaceDescriptorVector[i] << "\t";
			}
			interfaceDescriptorVectorLibraryFile << "\n";
		}
		else
			cout << "\nError reading the interface file " + interfaceFileName;
	}

	interfaceDescriptorVectorLibraryFile.close();
	clustersFile.close();
	interfacesList.close();
}


//****************************************************************************
/* This process extracts the overlapping fragments from all the interfaces */

int ExtractFragments(int fragmentLength, int overlappingResidues)
{

	int numberOfFragments = 0;
	string line;
	//int incrementValue = fragmentLength - overlappingResidues;
	string *fragment;
	fragment = new string[fragmentLength];
	int lineInFragment = -1;
	ifstream interfacesFile;
	ofstream fragmentsFile;
	bool continuousResidues;


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
		fragmentsFile.open(("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Protein_Fragments\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "\\fragments.txt").c_str());
		if (fragmentsFile.fail())
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
				fragmentsFile << "Fragment_" << numberOfFragments << "\n";
				for (int i = 0; i < fragmentLength; i++)
					fragmentsFile << fragment[i] << "\n";


				numberOfFragments++;

				lineInFragment = overlappingResidues;
				for (int i = 0; i < overlappingResidues; i++)
					fragment[i] = fragment[fragmentLength - overlappingResidues + i];

			}
		}
	}
	interfacesFile.close();
	fragmentsFile.close();
	cout << "\nNumber of fragments: " << numberOfFragments;
	return numberOfFragments;
}


//****************************************************************************
/*
This process randomly selects a seed fragment as a new cluster and puts the similar fragments in that cluster.
It continues till all the fragments are assigned to a cluster.
*/

int ClusterFragments(int fragmentLength, int numberOfFragments, int overlappingResidues, int expectedMatchedPoints, int binSize, int minNumberOfClusterMembers)
{

	string line;
	string *lines;
	lines = new string[fragmentLength];
	int ln;
	bool isSeed = true;
	int numberOfClusters = 0, numberOfSimilarFragments = 0;
	string similarFragments;
	int seedFragmentNumber;


	unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> seedHashTable;
	unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> nonSeedHashTable;
	ofstream clusteringFile;
	ifstream fragmentsFile;

	vector <int> unclusteredFragments(numberOfFragments);	//Create a vector of fragments which are all unclustered now
	for (int i = 0; i < numberOfFragments; i++)
	{
		unclusteredFragments[i] = i;
	}


	//Creates a file to put the clustering information there
	try
	{
		clusteringFile.open(("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\ClusteringLog.txt").c_str());
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


	//read the fragments file
	try
	{
		fragmentsFile.open(("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Protein_Fragments\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "\\fragments.txt").c_str());
		if (fragmentsFile.fail())
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
	catch (...)
	{
		cout << "\nOther exception thrown." << endl;
		return -1;
	}



	//Clustering the fragments
	srand((unsigned int)time(NULL));
	while (unclusteredFragments.size() > 0)
	{
		//Randomly select one fragment to be the seed
		int seed = rand() % unclusteredFragments.size();
		seedFragmentNumber = unclusteredFragments[seed];
		cout << "\nSeed fragment: " << seedFragmentNumber << "\n";
		//cout << "\nunclusteredFragments.size: " << unclusteredFragments.size() << "\n";
		similarFragments = " ";
		numberOfSimilarFragments = 0;

		unclusteredFragments.erase(unclusteredFragments.begin() + seed);

		ln = 0;
		int positionInFile = seed * (fragmentLength + 1) + 1;
		//cout << "\nPosition in File: " << positionInFile;
		int index = 0;
		while (getline(fragmentsFile, line))	//read the fragment lines from the fragments file
		{
			index++;
			if (index > positionInFile)
				lines[ln++] = line;
			if (index == positionInFile + fragmentLength)
				break;
		}

		fragmentsFile.clear();
		fragmentsFile.seekg(0, ios::beg);

		seedHashTable.clear();

		/* Select three consecutive points as Reference Set (RS) and create the hashtable*/
		for (int c = 0; c < fragmentLength - 2; c++)
			for (int i = c; i < c + 3; i++)
				for (int j = c; j < c + 3; j++)
					for (int k = c; k < c + 3; k++)
						if (i != j && i != k && j != k)	//three different points to create RS
						{

							float x1 = (float)atof(lines[i].substr(30, 8).c_str());
							float y1 = (float)atof(lines[i].substr(38, 8).c_str());
							float z1 = (float)atof(lines[i].substr(46, 8).c_str());
							//cout << "x1: " << x1 << "\ty1: " << y1 << "\tz1: " << z1 << "\n";

							float x2 = (float)atof(lines[j].substr(30, 8).c_str());
							float y2 = (float)atof(lines[j].substr(38, 8).c_str());
							float z2 = (float)atof(lines[j].substr(46, 8).c_str());
							//cout << "x2: " << x2 << "\ty2 " << y2 << "\tz2: " << z2 << "\n";

							float x3 = (float)atof(lines[k].substr(30, 8).c_str());
							float y3 = (float)atof(lines[k].substr(38, 8).c_str());
							float z3 = (float)atof(lines[k].substr(46, 8).c_str());
							//cout << "x3: " << x3 << "\ty3 " << y3 << "\tz3: " << z3 << "\n";

							TranslationParameter selectedRS;
							selectedRS = CalculateGeoTranslation(x1, y1, z1, x2, y2, z2, x3, y3, z3, binSize);


							if (selectedRS.Rx != 0 || selectedRS.Ry != 0 || selectedRS.Rz != 0)
							{
								/* Add the second and third points of RS to the hash table */
								AddToHashTable(seedHashTable, selectedRS.p2, i + 1, j + 1, k + 1);
								AddToHashTable(seedHashTable, selectedRS.p3, i + 1, j + 1, k + 1);

								/* Add non-RS points of the fragment to the hash table */
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

		//print the hash table
		/*cout << "-------------- Printing Hash Table ---------------" << endl << endl;
		for (auto it1 = seedHashTable.begin(); it1 != seedHashTable.end(); ++it1)
		{
		cout << it1->first << ":";
		for (auto it2 = seedHashTable[it1->first].begin(); it2 != seedHashTable[it1->first].end(); ++it2)
		{
		cout << "   " << it2->first << ":";
		for (auto it3 = seedHashTable[it1->first][it2->first].begin(); it3 != seedHashTable[it1->first][it2->first].end(); ++it3)
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
		}*/








		/* Compare the seed hashtable with all other unclustered fragments hashtables */
		int fragmentNo = 0;
		while (fragmentNo < (int)unclusteredFragments.size())
		{
			//cout << "unclusteredFragments.size: " << unclusteredFragments.size() << "\n" << "fragmentNo" << fragmentNo << "\n";

			ln = 0;
			index = 0;
			positionInFile = unclusteredFragments[fragmentNo] * (fragmentLength + 1) + 1;
			//cout << endl << "Fragment no: " << unclusteredFragments[fragmentNo] << "\tPosition in file: " << positionInFile;
			while (getline(fragmentsFile, line))	//read the fragment lines from the fragments file
			{
				index++;
				if (index > positionInFile)
				{
					lines[ln++] = line;
					//cout << endl << lines[ln - 1];
				}
				if (index == positionInFile + fragmentLength)
					break;
			}

			fragmentsFile.clear();
			fragmentsFile.seekg(0, ios::beg);


			bool matched = false;	//if the fragment hashtable and seed hash table matches
			for (int c = 0; c < fragmentLength - 2; c++)	// we select three consecutive points as Reference Set (RS)
			{
				for (int i = c; i < c + 3; i++)
				{
					for (int j = c; j < c + 3; j++)
					{
						for (int k = c; k < c + 3; k++)
							if (i != j && i != k && j != k)	//three different points to create RS
							{

								nonSeedHashTable.clear();

								float x1 = (float)atof(lines[i].substr(30, 8).c_str());
								float y1 = (float)atof(lines[i].substr(38, 8).c_str());
								float z1 = (float)atof(lines[i].substr(46, 8).c_str());
								//cout << "x1: " << x1 << "\ty1: " << y1 << "\tz1: " << z1 << "\n";

								float x2 = (float)atof(lines[j].substr(30, 8).c_str());
								float y2 = (float)atof(lines[j].substr(38, 8).c_str());
								float z2 = (float)atof(lines[j].substr(46, 8).c_str());
								//cout << "x2: " << x2 << "\ty2 " << y2 << "\tz2: " << z2 << "\n";

								float x3 = (float)atof(lines[k].substr(30, 8).c_str());
								float y3 = (float)atof(lines[k].substr(38, 8).c_str());
								float z3 = (float)atof(lines[k].substr(46, 8).c_str());
								//cout << "x3: " << x3 << "\ty3 " << y3 << "\tz3: " << z3 << "\n";

								TranslationParameter selectedRS;
								selectedRS = CalculateGeoTranslation(x1, y1, z1, x2, y2, z2, x3, y3, z3, binSize);

								if (selectedRS.Rx != 0 || selectedRS.Ry != 0 || selectedRS.Rz != 0)
								{
									AddToHashTable(nonSeedHashTable, selectedRS.p2, i + 1, j + 1, k + 1);
									AddToHashTable(nonSeedHashTable, selectedRS.p3, i + 1, j + 1, k + 1);

									for (int f = 0; f < fragmentLength; f++)
									{
										if (f != i && f != j && f != k)
										{
											Point p;
											p.x = (float)atof(lines[f].substr(30, 8).c_str());
											p.y = (float)atof(lines[f].substr(38, 8).c_str());
											p.z = (float)atof(lines[f].substr(46, 8).c_str());
											p = CalculateNewPoint(selectedRS, p, binSize);
											AddToHashTable(nonSeedHashTable, p, i + 1, j + 1, k + 1);
										}
									}

									//compares two hash tables
									matched = CompareTwoHashTables(seedHashTable, nonSeedHashTable, expectedMatchedPoints);
									if (matched) // if the seed fragment and the new one are similar
										break;
								}
							}
						if (matched) break;
					}// for j
					if (matched) break;
				}// for i
				if (matched) break;
			}//for c

			 /*cout << "-------------- Printing Second Hash Table ---------------" << endl << endl;
			 for (auto it1 = nonCentroidHashTable.begin(); it1 != nonCentroidHashTable.end(); ++it1)
			 {
			 cout << it1->first << ":";
			 for (auto it2 = nonCentroidHashTable[it1->first].begin(); it2 != nonCentroidHashTable[it1->first].end(); ++it2)
			 {
			 cout << "   " << it2->first << ":";
			 for (auto it3 = nonCentroidHashTable[it1->first][it2->first].begin(); it3 != nonCentroidHashTable[it1->first][it2->first].end(); ++it3)
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
			 }*/

			if (matched)
			{
				numberOfSimilarFragments++;
				similarFragments += to_string(unclusteredFragments[fragmentNo]) + "\t";
				//cout << "\n" << similarFragments;
				unclusteredFragments.erase(unclusteredFragments.begin() + fragmentNo);
			}
			else
			{
				fragmentNo++;
			}

		}
		cout << "\nNumber of cluster members: " << numberOfSimilarFragments;
		if (numberOfSimilarFragments > minNumberOfClusterMembers)
		{
			clusteringFile << "\nCluster_" << ++numberOfClusters << "\t" << seedFragmentNumber << "\t" << similarFragments;
			int res = WriteHashTableToFile(("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\Clusters.txt").c_str(), numberOfClusters, seedHashTable);
			if (res == -1)
				break;
		}
	}
	clusteringFile.close();
	fragmentsFile.close();
	cout << "\nNumber of clusters: " << numberOfClusters;
	return numberOfClusters;

}


//****************************************************************************
/*
This process randomly selects a seed fragment as a new cluster and puts
the similar fragments in that cluster. It continues till all the fragments
are assigned to a cluster.
This version reads 85 characters lines from the file. It randomly
accesses the file so it is significantly faster than the orijinal version.
*/

int ClusterFragments_v1(int fragmentLength, int numberOfFragments, int overlappingResidues, int expectedMatchedPoints, int binSize, int minNumberOfClusterMembers)
{

	int numberOfClusters = 0;
	int numberOfSimilarFragments = 0;
	string *lines;
	lines = new string[fragmentLength];
	string similarFragments, line;
	int seedFragmentNumber;
	char rd_buffer[85];


	unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> seedHashTable;
	unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> nonSeedHashTable;
	vector <int> unclusteredFragments;	//A vector of fragments which are all unclustered yet

	FILE * fragmentsFile;
	ofstream clusteringFile;
	fstream unclusteredFragmentsFile;


	//Reads the unclustered fragments file or creates it if doesn't exist
	try
	{
		unclusteredFragmentsFile.open(("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\UnclusteredFragments.txt").c_str(), fstream::in | fstream::out | fstream::app);
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
		fragmentsFile = fopen(("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Protein_Fragments\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "\\fragments.txt").c_str(), "r");
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
	catch (...)
	{
		cout << "\nOther exception thrown." << endl;
		return -1;
	}



	//Reads the clustering information file or creates it if doesn't exist
	try
	{
		clusteringFile.open(("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\ClusteringLog.txt").c_str(), ios::out | ios::app);
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

		//unclusteredFragmentsFile.seekg(0, unclusteredFragmentsFile.beg);
		for (int i = 0; i < numberOfFragments; i++)
		{
			//unclusteredFragments[i] = i;
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


	//Clustering the fragments
	while (unclusteredFragments.size() > 0)
	{
		int seed = rand() % unclusteredFragments.size();	//Randomly select one fragment to be the seed


		seed = 743643;
		seedFragmentNumber = unclusteredFragments[seed];



		cout << "\nSeed fragment: " << seedFragmentNumber << "\t";
		//cout << "\nunclusteredFragments.size: " << unclusteredFragments.size() << "\n";
		similarFragments = " ";
		numberOfSimilarFragments = 0;

		unclusteredFragments.erase(unclusteredFragments.begin() + seed);

		int positionInFile = seed * (fragmentLength * 85);
		//cout << "\nPosition in File: " << positionInFile;
		fseek(fragmentsFile, positionInFile, SEEK_SET);

		for (int index = 0; index < fragmentLength; index++)
		{
			if (fread(rd_buffer, sizeof(char), 85, fragmentsFile) < 85)
			{
				cout << "\nError reading the fragments file.";
				return -1;
			}
			//memcpy(&lines[index], rd_buffer, strlen(rd_buffer));
			lines[index].assign(&rd_buffer[0], &rd_buffer[0] + 85);
			//cout << endl << lines[index];
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
							//cout << "x1: " << x1 << "\ty1: " << y1 << "\tz1: " << z1 << "\n";

							float x2 = (float)atof(lines[j].substr(30, 8).c_str());
							float y2 = (float)atof(lines[j].substr(38, 8).c_str());
							float z2 = (float)atof(lines[j].substr(46, 8).c_str());
							//cout << "x2: " << x2 << "\ty2 " << y2 << "\tz2: " << z2 << "\n";

							float x3 = (float)atof(lines[k].substr(30, 8).c_str());
							float y3 = (float)atof(lines[k].substr(38, 8).c_str());
							float z3 = (float)atof(lines[k].substr(46, 8).c_str());
							//cout << "x3: " << x3 << "\ty3 " << y3 << "\tz3: " << z3 << "\n";

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

		//print the hash table
		/*cout << "-------------- Printing Hash Table ---------------" << endl << endl;
		for (auto it1 = seedHashTable.begin(); it1 != seedHashTable.end(); ++it1)
		{
		cout << it1->first << ":";
		for (auto it2 = seedHashTable[it1->first].begin(); it2 != seedHashTable[it1->first].end(); ++it2)
		{
		cout << "   " << it2->first << ":";
		for (auto it3 = seedHashTable[it1->first][it2->first].begin(); it3 != seedHashTable[it1->first][it2->first].end(); ++it3)
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
		}*/




		//Compare the seed hashtable with all other unclustered fragments hashtables
		int fragmentNo = 0;
		while (fragmentNo < (int)unclusteredFragments.size())
		{
			//cout << "unclusteredFragments.size: " << unclusteredFragments.size() << "\n" << "fragmentNo" << fragmentNo << "\n";

			positionInFile = unclusteredFragments[fragmentNo] * (fragmentLength * 85);
			//cout << endl << "Fragment no: " << unclusteredFragments[fragmentNo] << "\tPosition in file: " << positionInFile;
			fseek(fragmentsFile, positionInFile, SEEK_SET);

			for (int index = 0; index < fragmentLength; index++)
			{
				if (fread(rd_buffer, sizeof(char), 85, fragmentsFile) < 85)
				{
					cout << "\nError reading the fragments file.";
					return -1;
				}
				//memcpy(&lines[index], rd_buffer, strlen(rd_buffer));
				lines[index].assign(&rd_buffer[0], &rd_buffer[0] + 85);
			}


			bool matched = false;	//if the fragment hashtable and seed hash table matches
			for (int c = 0; c < fragmentLength - 2; c++)	// we select three consecutive points as Reference Set (RS)
			{
				for (int i = c; i < c + 3; i++)
				{
					for (int j = c; j < c + 3; j++)
					{
						for (int k = c; k < c + 3; k++)
							if (i != j && i != k && j != k)	//three different points to create RS
							{

								nonSeedHashTable.clear();

								float x1 = (float)atof(lines[i].substr(30, 8).c_str());
								float y1 = (float)atof(lines[i].substr(38, 8).c_str());
								float z1 = (float)atof(lines[i].substr(46, 8).c_str());
								//cout << "x1: " << x1 << "\ty1: " << y1 << "\tz1: " << z1 << "\n";

								float x2 = (float)atof(lines[j].substr(30, 8).c_str());
								float y2 = (float)atof(lines[j].substr(38, 8).c_str());
								float z2 = (float)atof(lines[j].substr(46, 8).c_str());
								//cout << "x2: " << x2 << "\ty2 " << y2 << "\tz2: " << z2 << "\n";

								float x3 = (float)atof(lines[k].substr(30, 8).c_str());
								float y3 = (float)atof(lines[k].substr(38, 8).c_str());
								float z3 = (float)atof(lines[k].substr(46, 8).c_str());
								//cout << "x3: " << x3 << "\ty3 " << y3 << "\tz3: " << z3 << "\n";

								TranslationParameter selectedRS;
								selectedRS = CalculateGeoTranslation(x1, y1, z1, x2, y2, z2, x3, y3, z3, binSize);

								if (selectedRS.Rx != 0 || selectedRS.Ry != 0 || selectedRS.Rz != 0)
								{
									AddToHashTable(nonSeedHashTable, selectedRS.p2, i + 1, j + 1, k + 1);
									AddToHashTable(nonSeedHashTable, selectedRS.p3, i + 1, j + 1, k + 1);

									for (int f = 0; f < fragmentLength; f++)
									{
										if (f != i && f != j && f != k)
										{
											Point p;
											p.x = (float)atof(lines[f].substr(30, 8).c_str());
											p.y = (float)atof(lines[f].substr(38, 8).c_str());
											p.z = (float)atof(lines[f].substr(46, 8).c_str());
											p = CalculateNewPoint(selectedRS, p, binSize);
											AddToHashTable(nonSeedHashTable, p, i + 1, j + 1, k + 1);
										}
									}

									//compares two hash tables
									matched = CompareTwoHashTables(seedHashTable, nonSeedHashTable, expectedMatchedPoints);
									if (matched) // if the seed fragment and the new one are similar
										break;
								}
							}
						if (matched) break;
					}// for j
					if (matched) break;
				}// for i
				if (matched) break;
			}//for c

			 /*cout << "-------------- Printing Second Hash Table ---------------" << endl << endl;
			 for (auto it1 = nonCentroidHashTable.begin(); it1 != nonCentroidHashTable.end(); ++it1)
			 {
			 cout << it1->first << ":";
			 for (auto it2 = nonCentroidHashTable[it1->first].begin(); it2 != nonCentroidHashTable[it1->first].end(); ++it2)
			 {
			 cout << "   " << it2->first << ":";
			 for (auto it3 = nonCentroidHashTable[it1->first][it2->first].begin(); it3 != nonCentroidHashTable[it1->first][it2->first].end(); ++it3)
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
			 }*/

			if (matched)
			{
				numberOfSimilarFragments++;
				similarFragments += to_string(unclusteredFragments[fragmentNo]) + "\t";
				//cout << "\n" << similarFragments;
				unclusteredFragments.erase(unclusteredFragments.begin() + fragmentNo);
			}
			else
			{
				fragmentNo++;
			}

		}

		cout << "Number of cluster members: " << numberOfSimilarFragments + 1;


		unclusteredFragmentsFile.open(("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\UnclusteredFragments.txt").c_str(), fstream::out | fstream::app);
		unclusteredFragmentsFile << "SeedFragment:" << seedFragmentNumber << "\t NumberOfUnclusteredFragments:" << unclusteredFragments.size();
		for (int i = 0; i < (int)unclusteredFragments.size(); i++)	//Writes the unclustered fragments to file
			unclusteredFragmentsFile << "\t" << unclusteredFragments[i];
		unclusteredFragmentsFile << "\n";
		unclusteredFragmentsFile.close();


		if (numberOfSimilarFragments + 1 >= minNumberOfClusterMembers)
		{
			clusteringFile << "\nCluster_" << ++numberOfClusters << "\t" << seedFragmentNumber << "\t" << similarFragments;
			int res = WriteHashTableToFile(("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\Clusters.txt").c_str(), numberOfClusters, seedHashTable);
			if (res == -1)
				break;
		}
	}

	clusteringFile.close();
	fclose(fragmentsFile);
	cout << "\nNumber of clusters: " << numberOfClusters;
	return numberOfClusters;

}


//****************************************************************************
/*
This process creates the descriptor vectors for the template interface partners of PRISM as compact dictionaries
*/

int CreateInterfaceDescriptors_v2(int fragmentLength, int overlappingResidues, int binSize, int expectedMatchedPoints, int numberOfClusters)
{

	std::map<int, int> interfaceDescriptor;

	string line, interfaceFileName;
	int incrementValue = fragmentLength - overlappingResidues;
	//bool continuousResidues;	//To see if the extracted fragment is continuous or not

	ifstream interfacesList, clustersFile;
	ofstream interfaceDescriptorsLibraryFile;


	//Reads the interfaces list
	try
	{
		interfacesList.open("D:\\PhD\\prism\\prism_standalone\\template\\InterfaceList.txt");

		if (interfacesList.fail())
		{
			cout << "\nError reading the interface list file.";
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


	//Creates the interface descriptors library file
	try
	{
		interfaceDescriptorsLibraryFile.open("D:\\PhD\\My_Thesis\\Second_Step\\Data\\PRISM_Interface_Descriptors\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\PrismInterfaceDescriptors.txt");

		if (interfaceDescriptorsLibraryFile.fail())
		{
			cout << "\nError creating the interface vector library file.";
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
		cout << "\nOther exception thrown";
		return -1;
	}


	//Reads the fragment clusters file
	try
	{
		clustersFile.open("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Clustering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\Clusters.txt");

		if (clustersFile.fail())
		{
			cout << "\nError reading the cluster hash tables file.";
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


	getline(clustersFile, line); //To pass the line having the cluster number


	while (getline(interfacesList, interfaceFileName)) //each line has an interface partner name
	{

		//interfaceFileName = interfaceFileName.substr(0, interfaceFileName.size() - 1);
		cout << endl << interfaceFileName;

		ifstream interfaceFile(("D:\\PhD\\prism\\prism_standalone\\template\\interfaces\\" + interfaceFileName).c_str());	//read the interface file

		if (interfaceFile)
		{

			vector <string> interfaceResidueCoordinates;
			while (getline(interfaceFile, line))	//read the lines of the interface partner file
				if (line.find("ATOM") == 0)
					interfaceResidueCoordinates.push_back(line);

			interfaceFile.close();

			interfaceDescriptorsLibraryFile << interfaceFileName;

			if ((int)interfaceResidueCoordinates.size() >= fragmentLength)	//if the interface file has at least the size of fragmentLength 
			{
				interfaceDescriptor.clear();

				unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> fragmentHashTable;

				for (int lineIndex1 = 0; lineIndex1 <= (int)interfaceResidueCoordinates.size() - (fragmentLength); lineIndex1 = lineIndex1 + incrementValue)	//extracts continuous fragments from interfaces
				{

					int previousResidueNumber = atoi(interfaceResidueCoordinates[lineIndex1].substr(22, 5).c_str());
					bool continuousResidues = true;

					for (int lineIndex2 = (lineIndex1 + 1); lineIndex2 < (lineIndex1 + fragmentLength); lineIndex2++)
					{

						int nextResidueNumber = atoi(interfaceResidueCoordinates[lineIndex2].substr(22, 5).c_str());

						if (previousResidueNumber < nextResidueNumber - 1)
						{
							continuousResidues = false;
							lineIndex1 = lineIndex2 - incrementValue;///////// As icrement value is added in the first loop
							break;
						}
						previousResidueNumber = nextResidueNumber;
					}


					if (continuousResidues == true)	//We made sure that it is a continuous fragment. Now we want to find which fragment in our fragment library is similar to this fragment
					{

						unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> clusterHashTable;
						bool matched = false;
						int scannedClusters = 1;

						while (!matched && scannedClusters <= numberOfClusters) //till we find which cluster hash table is similar to the fragment hashtable or none of them are similar
						{
							//cout << "\n" << scannedClusters;

							LoadHashTableFromFile(clustersFile, scannedClusters, clusterHashTable);


							for (int c = lineIndex1; c < (lineIndex1 + fragmentLength - 2); c++)	// we select three consecutive points as Reference Set (RS)
							{
								for (int i = c; i < c + 3; i++)
								{
									for (int j = c; j < c + 3; j++)
									{
										for (int k = c; k < c + 3; k++)
											if (i != j && i != k && j != k)	//three different points to create RS
											{

												fragmentHashTable.clear();

												float x1 = (float)atof(interfaceResidueCoordinates[i].substr(30, 8).c_str());
												float y1 = (float)atof(interfaceResidueCoordinates[i].substr(38, 8).c_str());
												float z1 = (float)atof(interfaceResidueCoordinates[i].substr(46, 8).c_str());

												float x2 = (float)atof(interfaceResidueCoordinates[j].substr(30, 8).c_str());
												float y2 = (float)atof(interfaceResidueCoordinates[j].substr(38, 8).c_str());
												float z2 = (float)atof(interfaceResidueCoordinates[j].substr(46, 8).c_str());

												float x3 = (float)atof(interfaceResidueCoordinates[k].substr(30, 8).c_str());
												float y3 = (float)atof(interfaceResidueCoordinates[k].substr(38, 8).c_str());
												float z3 = (float)atof(interfaceResidueCoordinates[k].substr(46, 8).c_str());

												TranslationParameter selectedRS;
												selectedRS = CalculateGeoTranslation(x1, y1, z1, x2, y2, z2, x3, y3, z3, binSize);

												AddToHashTable(fragmentHashTable, selectedRS.p2, i + 1, j + 1, k + 1);
												AddToHashTable(fragmentHashTable, selectedRS.p3, i + 1, j + 1, k + 1);

												for (int f = lineIndex1; f < fragmentLength + lineIndex1; f++)
												{
													if (f != i && f != j && f != k)
													{

														Point p;
														p.x = (float)atof(interfaceResidueCoordinates[f].substr(30, 8).c_str());
														p.y = (float)atof(interfaceResidueCoordinates[f].substr(38, 8).c_str());
														p.z = (float)atof(interfaceResidueCoordinates[f].substr(46, 8).c_str());
														p = CalculateNewPoint(selectedRS, p, binSize);
														AddToHashTable(fragmentHashTable, p, i + 1, j + 1, k + 1);
													}
												}


												matched = CompareTwoHashTables(clusterHashTable, fragmentHashTable, expectedMatchedPoints);	//compares two hash tables
												if (matched) //if the cluster hash table and the new one are similar
													break;
											}
										if (matched) break;
									}// for j
									if (matched) break;
								}// for i
								if (matched) break;
							}//for c

							if (matched)
								if (interfaceDescriptor.find(scannedClusters) != interfaceDescriptor.end())	// found
									interfaceDescriptor[scannedClusters]++;
								else // not found
									interfaceDescriptor[scannedClusters] = 1;

							scannedClusters++;
						}

						clustersFile.clear();
						clustersFile.seekg(0, ios::beg);
						getline(clustersFile, line);	//To pass the line having cluster number
					}
				}

				interfaceDescriptorsLibraryFile << "\t";
				for (auto it1 = interfaceDescriptor.begin(); it1 != interfaceDescriptor.end(); ++it1)
					interfaceDescriptorsLibraryFile << it1->first << ":" << it1->second << ",";
			}
			else
				cout << "\n" << interfaceFileName << "interface file is shorter than fragment length.";

			interfaceDescriptorsLibraryFile << "\n";
		}
		else
			cout << "\nError reading the interface file " + interfaceFileName;
	}

	interfaceDescriptorsLibraryFile.close();
	clustersFile.close();
	interfacesList.close();
}



//****************************************************************************
/* This process loads a hash table from a file */

void LoadHashTableFromFile(ifstream &clustersFile, int clusterNumber, unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> &hashTable)
{

	string line;
	//bool clusterFound = false;

	hashTable.clear();
	//clustersFile.clear();
	//clustersFile.seekg(0, ios::beg);


	while (getline(clustersFile, line))
	{
		//if (!clusterFound && line.find("Cluster") == 0 && stoi(line.substr(8)) == clusterNumber)
		//clusterFound = true;

		if (line.find("Cluster") != 0)
		{
			istringstream iss(line);
			vector <string> lineTokens{ istream_iterator<string>{iss}, istream_iterator<string>{} };

			Point p;
			p.x = (float)stoi(lineTokens[0]);
			p.y = (float)stoi(lineTokens[1]);
			p.z = (float)stoi(lineTokens[2]);

			//if (hashTable.find(p.x) != hashTable.end() && hashTable[p.x].find(p.y) != hashTable[p.x].end() && hashTable[p.x][p.y].find(p.z) != hashTable[p.x][p.y].end())	// found
			//hashTable[p.x][p.y][p.z].insert(hashTable[p.x][p.y][p.z].end(), lineTokens[3]);
			//hashTable[p.x][p.y][p.z].push_back(lineTokens[3]);
			//else // not found
			//hashTable[p.x][p.y][p.z] = { lineTokens[3] };

			hashTable[p.x][p.y][p.z].push_back(lineTokens[3]);

			for (int i = 4; i < (int)lineTokens.size(); i++)
			{
				//hashTable[p.x][p.y][p.z].insert(hashTable[p.x][p.y][p.z].end(), lineTokens[i]);
				hashTable[p.x][p.y][p.z].push_back(lineTokens[i]);
			}
		}
		//else if (clusterFound && line.find("Cluster") == 0)
		else
			break;
	}




	/*	ofstream fp("D:\\PhD\\My_Thesis\\Second_Step\\Code_Outputs\\Clustering_Results\\FrLn4_OvRd3_MtPn3_BnSz3\\Cluster_16_New.txt");

	if (!fp)
	{
	cout << "\nError creating the hash table file.";
	}
	else
	{

	for (auto it1 = hashTable.begin(); it1 != hashTable.end(); ++it1)
	for (auto it2 = hashTable[it1->first].begin(); it2 != hashTable[it1->first].end(); ++it2)
	for (auto it3 = hashTable[it1->first][it2->first].begin(); it3 != hashTable[it1->first][it2->first].end(); ++it3)
	{
	fp << it1->first << "\t" << it2->first << "\t" << it3->first;
	for (int it4 = 0; it4 < it3->second.size(); it4++)
	fp << "\t" << it3->second[it4];
	fp << "\n";
	}

	fp.close();
	}*/
}