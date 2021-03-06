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


//*****************************************************************************************
/*
This function creates the descriptor vectors of PRISM interfaces as compact vectors showing
the frequency of each fragment in a paralle manner.
This version has less rotations i.e. (1,2,3) and (3,2,1)
*/

int CreateInterfaceDescriptors_v2_LessRotations(int fragmentLength, int overlappingResidues, int binSize, int expectedMatchedPoints, int numberOfClusters)
{

	string line;
	int incrementValue = fragmentLength - overlappingResidues;
	bool continuousResidues;	//To see if the extracted fragment is continuous or not

	ifstream interfacesListFile;
	FILE * clustersFile;
	char * clustersBuffer;
	int clustersFileSize;
	ofstream interfaceDescriptorsLibraryFile;
	vector <string> clusters;

	list<string> interfacesList;
	vector<unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>>> fragmentHashTables;
	unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> fHashTable;


	//Reads the interfaces list
	try
	{
		interfacesListFile.open(rootDir + "InterfaceList.txt");

		if (interfacesListFile.fail())
		{
			cout << "\nError reading the interface list file.";
			return -1;
		}
		else
		{
			while (getline(interfacesListFile, line))
			{
				interfacesList.push_back(line);
				//interfacesList.push_back(line.substr(0, line.size() - 1));	//linux version
			}
		}
		interfacesListFile.close();
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
		interfaceDescriptorsLibraryFile.open(rootDir + "Less_Rotations_PRISM_Interface_Descriptors" + slash + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + slash + "PrismInterfaceDescriptors.txt");

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
		//clustersFile.open(rootFolder + "Clustering_Results"+dirSepreator+"FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + dirSepreator+"Clusters.txt");
		clustersFile = fopen((rootDir + "Less_Rotations_Clustering_Results" + slash + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + slash + "Clusters.txt").c_str(), "rb");


		fseek(clustersFile, 0, SEEK_END);
		clustersFileSize = ftell(clustersFile);
		/*if (!GetFileSizeEx(clustersFile, cFileSize)) {
		/* Handle error
		}
		clustersFileSize=(int)cFileSize;*/
		clustersBuffer = (char *)malloc(clustersFileSize + 1);
		clustersBuffer[clustersFileSize] = 0;
		fseek(clustersFile, 0, SEEK_SET);
		int freadsize = 0;
		freadsize += fread(clustersBuffer + freadsize, 1, clustersFileSize - freadsize, clustersFile);
		if (clustersFileSize != freadsize)
		{
			cout << "\nError reading all of the cluster hash tables file.";
			return -1;
		}
		fclose(clustersFile);

		std::stringstream ss(clustersBuffer);
		std::string cluterLine;

		while (std::getline(ss, cluterLine, '\n'))
		{
			if (cluterLine.find("Cluster") == 0)
			{
				if (fHashTable.size()>0)
				{
					fragmentHashTables.push_back(fHashTable);
					fHashTable.clear();
				}
			}
			else
			{
				istringstream iss(cluterLine);
				vector <string> lineTokens{ istream_iterator<string>{iss}, istream_iterator<string>{} };

				Point p;
				p.x = (float)stoi(lineTokens[0]);
				p.y = (float)stoi(lineTokens[1]);
				p.z = (float)stoi(lineTokens[2]);


				fHashTable[p.x][p.y][p.z].push_back(lineTokens[3]);

				for (int i = 4; i < lineTokens.size(); i++)
				{
					fHashTable[p.x][p.y][p.z].push_back(lineTokens[i]);
				}
			}
		}
		if (fHashTable.size()>0)	//To get the last cluster
		{
			fragmentHashTables.push_back(fHashTable);
		}
		cout << "fragmentHashTables size : " << fragmentHashTables.size();
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



#pragma omp parallel
	{
		string interfaceFileName;
		string interfaceline;
		//string interfaceDescriptorsLibrary;
		map<int, int> interfaceDescriptor;
		int finish = 0;
		while (true) //each line has an interface partner name
		{
#pragma omp critical
			{
				if (interfacesList.empty())
				{
					finish = 1;
				}
				else
				{
					interfaceFileName = interfacesList.front();
					cout << endl << interfaceFileName;
					interfacesList.pop_front();
				}
			}

			if (finish)
			{
				break;
			}


			//interfaceDescriptorsLibrary = "";
			ifstream interfaceFile((rootDir + "interfaces" + slash + interfaceFileName).c_str());	//read the interface file

			if (interfaceFile)
			{

				vector <string> interfaceResidueCoordinates;
				while (getline(interfaceFile, interfaceline))	//read the lines of the interface partner file
					if (interfaceline.find("ATOM") == 0)
						interfaceResidueCoordinates.push_back(interfaceline);

				interfaceFile.close();

				//need buffer for parallel
				//interfaceDescriptorsLibraryFile << interfaceFileName;
				//interfaceDescriptorsLibrary += interfaceFileName;

				if (interfaceResidueCoordinates.size() >= fragmentLength)	//if the interface file has at least the size of frag,entLength 
				{
					interfaceDescriptor.clear();

					unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> fragmentHashTable;

					for (int lineIndex1 = 0; lineIndex1 <= interfaceResidueCoordinates.size() - (fragmentLength); lineIndex1 = lineIndex1 + incrementValue)	//extracts continuous fragments from interfaces
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

							while (!matched && scannedClusters <= fragmentHashTables.size()) //till we find which cluster hash table is similar to the fragment hashtable or none of them are similar
							{
								//cout << "\n" << scannedClusters;

								//LoadHashTableFromFile(clustersFile, scannedClusters, clusterHashTable);
								//LoadHashTableFromFile(clustersFilebuffer, scannedClusters, clusterHashTable);
								clusterHashTable = fragmentHashTables[scannedClusters - 1];
								int i, j, k;
								for (int c = lineIndex1; c < (lineIndex1 + fragmentLength - 2); c++)	// we select three consecutive points as Reference Set (RS)
								{

									i = c;
									j = c + 1;
									k = c + 2;
									fragmentHashTable.clear();

									float x1 = atof(interfaceResidueCoordinates[i].substr(30, 8).c_str());
									float y1 = atof(interfaceResidueCoordinates[i].substr(38, 8).c_str());
									float z1 = atof(interfaceResidueCoordinates[i].substr(46, 8).c_str());

									float x2 = atof(interfaceResidueCoordinates[j].substr(30, 8).c_str());
									float y2 = atof(interfaceResidueCoordinates[j].substr(38, 8).c_str());
									float z2 = atof(interfaceResidueCoordinates[j].substr(46, 8).c_str());

									float x3 = atof(interfaceResidueCoordinates[k].substr(30, 8).c_str());
									float y3 = atof(interfaceResidueCoordinates[k].substr(38, 8).c_str());
									float z3 = atof(interfaceResidueCoordinates[k].substr(46, 8).c_str());

									TranslationParameter selectedRS;
									selectedRS = CalculateGeoTranslation(x1, y1, z1, x2, y2, z2, x3, y3, z3, binSize);

									if (selectedRS.Rx != 0 || selectedRS.Ry != 0 || selectedRS.Rz != 0)
									{
										AddToHashTable(fragmentHashTable, selectedRS.p2, i + 1, j + 1, k + 1);
										AddToHashTable(fragmentHashTable, selectedRS.p3, i + 1, j + 1, k + 1);

										for (int f = lineIndex1; f < fragmentLength + lineIndex1; f++)
										{
											if (f != i && f != j && f != k)
											{

												Point p;
												p.x = atof(interfaceResidueCoordinates[f].substr(30, 8).c_str());
												p.y = atof(interfaceResidueCoordinates[f].substr(38, 8).c_str());
												p.z = atof(interfaceResidueCoordinates[f].substr(46, 8).c_str());
												p = CalculateNewPoint(selectedRS, p, binSize);
												AddToHashTable(fragmentHashTable, p, i + 1, j + 1, k + 1);
											}
										}


										matched = CompareTwoHashTables(clusterHashTable, fragmentHashTable, expectedMatchedPoints);	//compares two hash tables
										if (matched) //if the cluster hash table and the new one are similar
											break;
									}

									k = c;
									i = c + 2;

									fragmentHashTable.clear();

									x1 = atof(interfaceResidueCoordinates[i].substr(30, 8).c_str());
									y1 = atof(interfaceResidueCoordinates[i].substr(38, 8).c_str());
									z1 = atof(interfaceResidueCoordinates[i].substr(46, 8).c_str());

									x2 = atof(interfaceResidueCoordinates[j].substr(30, 8).c_str());
									y2 = atof(interfaceResidueCoordinates[j].substr(38, 8).c_str());
									z2 = atof(interfaceResidueCoordinates[j].substr(46, 8).c_str());

									x3 = atof(interfaceResidueCoordinates[k].substr(30, 8).c_str());
									y3 = atof(interfaceResidueCoordinates[k].substr(38, 8).c_str());
									z3 = atof(interfaceResidueCoordinates[k].substr(46, 8).c_str());

									//TranslationParameter selectedRS;
									selectedRS = CalculateGeoTranslation(x1, y1, z1, x2, y2, z2, x3, y3, z3, binSize);

									if (selectedRS.Rx != 0 || selectedRS.Ry != 0 || selectedRS.Rz != 0)
									{

										AddToHashTable(fragmentHashTable, selectedRS.p2, i + 1, j + 1, k + 1);
										AddToHashTable(fragmentHashTable, selectedRS.p3, i + 1, j + 1, k + 1);

										for (int f = lineIndex1; f < fragmentLength + lineIndex1; f++)
										{
											if (f != i && f != j && f != k)
											{

												Point p;
												p.x = atof(interfaceResidueCoordinates[f].substr(30, 8).c_str());
												p.y = atof(interfaceResidueCoordinates[f].substr(38, 8).c_str());
												p.z = atof(interfaceResidueCoordinates[f].substr(46, 8).c_str());
												p = CalculateNewPoint(selectedRS, p, binSize);
												AddToHashTable(fragmentHashTable, p, i + 1, j + 1, k + 1);
											}
										}


										matched = CompareTwoHashTables(clusterHashTable, fragmentHashTable, expectedMatchedPoints);	//compares two hash tables
										if (matched) //if the cluster hash table and the new one are similar
											break;
									}

								}//for c

								if (matched)
									if (interfaceDescriptor.find(scannedClusters) != interfaceDescriptor.end())	// found
										interfaceDescriptor[scannedClusters]++;
									else // not found
										interfaceDescriptor[scannedClusters] = 1;

								scannedClusters++;
							}
						}
					}

#pragma omp critical
					{
						interfaceDescriptorsLibraryFile << interfaceFileName;
						interfaceDescriptorsLibraryFile << "\t";
						for (auto it1 = interfaceDescriptor.begin(); it1 != interfaceDescriptor.end(); ++it1)
							interfaceDescriptorsLibraryFile << it1->first << ":" << it1->second << ",";
						interfaceDescriptorsLibraryFile << "\n";
					}
				}
				else
					cout << "\n" << interfaceFileName << "interface file is shorter than fragment length.";
			}
			else
				cout << "\nError reading the interface file " + interfaceFileName;


		}
	}

	interfaceDescriptorsLibraryFile.close();
}



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

/*void LoadHashTableFromFile(ifstream &clustersFile, int clusterNumber, unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> &hashTable)
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
	}
}
*/

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
		unclusteredFragmentsFile.open((rootDir + dataPath + "Less_Rotations_Clustering_Results" + slash + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + slash + "UnclusteredFragments.txt").c_str(), fstream::in | fstream::out | fstream::app);
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
		fragmentsFile = fopen((rootDir + "Protein_Fragments" + slash + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + slash + "fragments.txt").c_str(), "r");
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
		clusteringFile.open((rootDir + "Clustering_Results" + slash + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + slash + "ClusteringLog.txt"), ios::out | ios::app);
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

		unclusteredFragmentsFile.open((rootDir + "Clustering_Results" + slash + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + slash + "UnclusteredFragments.txt").c_str(), fstream::out | fstream::app);
		unclusteredFragmentsFile << "SeedFragment:" << seedFragmentNumber << "\t NumberOfUnclusteredFragments:" << unclusteredFragments.size();
		for (int i = 0; i < (int)unclusteredFragments.size(); i++)	//Writes the unclustered fragments to file
			unclusteredFragmentsFile << "\t" << unclusteredFragments[i];
		unclusteredFragmentsFile << "\n";
		unclusteredFragmentsFile.close();


		if (numberOfSimilarFragments + 1 >= minNumberOfClusterMembers)
		{
			clusteringFile << "\nCluster_" << ++numberOfClusters << "\t" << seedFragmentNumber << "\t" << similarFragments;
			int res = WriteHashTableToFile((rootDir + "Clustering_Results" + slash + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + slash + "Clusters.txt").c_str(), numberOfClusters, seedHashTable);
			if (res == -1)
				break;
		}
	}

	clusteringFile.close();
	cout << "\nNumber of clusters: " << numberOfClusters;
	return numberOfClusters;

}

//****************************************************************************************
/*
This process compares the descriptors of proteins in PRISM pair_list with interface descriptor vectors.
The descriptors are compressed. Serial results can be processed.
*/

int CompareProteinWithInterfaces_v1(int fragmentLength, int overlappingResidues, int binSize, int expectedMatchedPoints, int includeRatio)
{

	std::map<int, int> proteinSurfaceDescriptor_1, proteinSurfaceDescriptor_2;
	std::map<int, int> interfaceSideDescriptor_1, interfaceSideDescriptor_2;

	int key, value;

	string line, proteinName_1, proteinName_2, interfaceSideName_1 = "noname", interfaceSideName_2 = "noname";
	ifstream pairList, interfaceDescriptorsLibrary;
	ofstream interfacesSimilarToProteins;
	vector <string> proteinList;
	vector <string> potentialSimilarInterfaces;

	//Read the target proteins pair list
	try
	{
		pairList.open("D:\\PhD\\prism\\prism_standalone\\pair_list");

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


	//Create the output file to put interfaces similar to the proteins there	
	try
	{
		interfacesSimilarToProteins.open("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Filtering_Results\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\" + to_string(includeRatio) + "%InterfacesSimilarToProteinPairs.txt");

		if (interfacesSimilarToProteins.fail())
		{
			cout << "\nError creating the output file.";
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


	//Reads the interface descriptors library file
	try
	{
		interfaceDescriptorsLibrary.open("D:\\PhD\\My_Thesis\\Second_Step\\Data\\PRISM_Interface_Descriptors\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\PrismInterfaceDescriptors.txt"); //All PRISM interface descriptors

		if (interfaceDescriptorsLibrary.fail())
		{
			cout << "\nError reading the interface descriptors library file.";
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



	while (getline(pairList, line))	//each line has two protein names
	{
		int spac = line.find(" ", 0);

		proteinName_1 = line.substr(0, spac);
		cout << "\n" << proteinName_1;

		proteinName_2 = line.substr(spac + 1, line.length());
		cout << " - " << proteinName_2;


		ifstream proteinDescriptorFile_1("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Protein_Descriptor_Vectors\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\" + proteinName_1 + "DescriptorVector.txt");
		ifstream proteinDescriptorFile_2("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Protein_Descriptor_Vectors\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\" + proteinName_2 + "DescriptorVector.txt");
		if (proteinDescriptorFile_1 && proteinDescriptorFile_2)
		{

			interfacesSimilarToProteins << "Interfaces similar to " << proteinName_1 << "-" << proteinName_2 << ":\n";

			proteinSurfaceDescriptor_1.clear();
			proteinSurfaceDescriptor_2.clear();

			//Load protein_1 surface descriptor from file
			getline(proteinDescriptorFile_1, line);
			int length = line.length();
			int commaPos = -1;

			while (commaPos < length - 1)
			{
				key = stoi(line.substr(commaPos + 1, line.find(":", commaPos + 1)));
				value = stoi(line.substr(line.find(":", commaPos + 1) + 1, line.find(",", commaPos + 1)));
				//cout << "\nKey: " << key << "\tValue: " << value;
				proteinSurfaceDescriptor_1[key] = value;
				commaPos = line.find(",", commaPos + 1);
			}


			//Load protein_2 surface descriptor from file
			getline(proteinDescriptorFile_2, line);
			length = line.length();
			commaPos = -1;

			while (commaPos < length - 1)
			{
				key = stoi(line.substr(commaPos + 1, line.find(":", commaPos + 1)));
				value = stoi(line.substr(line.find(":", commaPos + 1) + 1, line.find(",", commaPos + 1)));
				//cout << "\nKey: " << key << "\tValue: " << value;
				proteinSurfaceDescriptor_2[key] = value;
				commaPos = line.find(",", commaPos + 1);
			}



			//Compare each protein pair with all interfaces
			int lineNumber = 0;
			bool side1Existence = false;

			while (getline(interfaceDescriptorsLibrary, line)) //each line has an interface descriptor
			{

				//lineNumber++;

				//if (lineNumber % 2 != 0)	//Read side 1 of the interface
				if (interfaceSideName_1 == "noname")
				{

					length = line.length();

					if (length > 13)
					{

						side1Existence = true;
						interfaceSideDescriptor_1.clear();
						interfaceSideName_1 = line.substr(0, 12);
						//cout << "\n" << "interfaceSideName_1: " << interfaceSideName_1;
						//getchar();

						commaPos = 12;

						while (commaPos < length - 1)	//Load interface descriptor from file
						{
							key = stoi(line.substr(commaPos + 1, line.find(":", commaPos + 1)));
							value = stoi(line.substr(line.find(":", commaPos + 1) + 1, line.find(",", commaPos + 1)));
							//cout << "\nKey: " << key << "\tValue: " << value;
							interfaceSideDescriptor_1[key] = value;
							commaPos = line.find(",", commaPos + 1);
						}

					}
					else
						side1Existence = false;

					//getchar();
				}

				else if (interfaceSideName_2 == "noname") //Read side 2 of the interface and do the comparisons
				{

					length = line.length();

					if (length > 13 && side1Existence)
					{

						interfaceSideName_2 = line.substr(0, 12);
						//cout << "\n" << "interfaceSideName_2: " << interfaceSideName_2;
						//getchar();


						if (interfaceSideName_1.substr(0, 6) == interfaceSideName_2.substr(0, 6))
						{

							interfaceSideDescriptor_2.clear();
							commaPos = 12;

							while (commaPos < length - 1)	//Load interface descriptor from file
							{
								key = stoi(line.substr(commaPos + 1, line.find(":", commaPos + 1)));
								value = stoi(line.substr(line.find(":", commaPos + 1) + 1, line.find(",", commaPos + 1)));
								//cout << "\nKey: " << key << "\tValue: " << value;
								interfaceSideDescriptor_2[key] = value;
								commaPos = line.find(",", commaPos + 1);
							}
							//getchar();


							bool interfaceSurfaceSimilarity = CompareProteinSurfaceWithInterfaceSide(proteinSurfaceDescriptor_1, interfaceSideDescriptor_1, (float)includeRatio);
							//cout << endl << interfaceSurfaceSimilarity;
							if (interfaceSurfaceSimilarity)
							{
								//cout << endl << "Prtein A is similar to interface side 1.";
								interfaceSurfaceSimilarity = CompareProteinSurfaceWithInterfaceSide(proteinSurfaceDescriptor_2, interfaceSideDescriptor_2, (float)includeRatio);
								if (interfaceSurfaceSimilarity)
								{
									interfacesSimilarToProteins << interfaceSideName_1.substr(0, 6) << "\n";
									//cout << endl << "Prtein B is similar to interface side 2.";
								}
							}

							//else
							if (!interfaceSurfaceSimilarity)
							{
								interfaceSurfaceSimilarity = CompareProteinSurfaceWithInterfaceSide(proteinSurfaceDescriptor_1, interfaceSideDescriptor_2, (float)includeRatio);
								if (interfaceSurfaceSimilarity)
								{
									//cout << endl << "Prtein A is similar to interface side 2.";
									interfaceSurfaceSimilarity = CompareProteinSurfaceWithInterfaceSide(proteinSurfaceDescriptor_2, interfaceSideDescriptor_1, (float)includeRatio);
									if (interfaceSurfaceSimilarity)
									{
										interfacesSimilarToProteins << interfaceSideName_1.substr(0, 6) << "\n";
										//cout << endl << "Prtein B is similar to interface side 1.";
									}
								}
							}
							//cout << endl << interfaceSurfaceSimilarity;
							//getchar();
							interfaceSideName_1 = "noname";
						}
						else
							interfaceSideName_1 = interfaceSideName_2;

						interfaceSideName_2 = "noname";

					}

				}
			}

			interfaceDescriptorsLibrary.clear();
			interfaceDescriptorsLibrary.seekg(0, ios::beg);

		}
		else
			cout << "\nError reading the protein surface descriptor file for protein: " + proteinName_1 + "or " + proteinName_2;
	}

	interfacesSimilarToProteins.close();
	pairList.close();
	interfaceDescriptorsLibrary.close();
	return 0;
}

//****************************************************************************
/*
This process creates descriptor vectors for protein in PRISM pair_list from which the surfaces are extracted.
It creates the descriptor vector as a map just showing the frequency of the existent fragments, not as a sparce vector.
*/

int CreateProteinDescriptors_v2(string jobName, int fragmentLength, int overlappingResidues, int binSize, int expectedMatchedPoints, int numberOfClusters)
{

	map<int, int> proteinSurfaceDescriptor;

	string line, proteinName;
	string pairListPath = "D:/PhD/prism/prism_standalone/pair_list_new";	//Linux is sensitive to upper cases
	int incrementValue = fragmentLength - overlappingResidues;
	vector <string> proteinList;

	ifstream pairListFile, clustersFile;

	//Read the pair_list file
	try
	{
		pairListFile.open(pairListPath); //list of template proteins

		if (pairListFile.fail())
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
		clustersFile.open(rootDir + "Clustering_Results/FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "/Clusters.txt");
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



	while (getline(pairListFile, line))	//each line has two protein names
	{

		int spac = line.find(" ", 0);

		proteinName = line.substr(0, spac);
		if (!any_of(proteinList.begin(), proteinList.end(), bind2nd(std::equal_to<std::string>(), proteinName)))	// not found
			proteinList.push_back(proteinName);

		proteinName = line.substr(spac + 1, line.length());
		proteinName.erase(remove_if(proteinName.begin(),
			proteinName.end(),
			[](unsigned char x) {return isspace(x); }),
			proteinName.end());

		if (!any_of(proteinList.begin(), proteinList.end(), bind2nd(std::equal_to<std::string>(), proteinName)))	// not found
			proteinList.push_back(proteinName);
	}

	pairListFile.close();


	for (int proteinIndex = 0; proteinIndex < (int)proteinList.size(); ++proteinIndex)
	{

		proteinName = proteinList[proteinIndex];

		ifstream proteinFile(("D:/PhD/prism/prism_standalone/jobs/" + jobName + "/surfaceExtract/" + proteinName + ".asa.pdb").c_str());	//read the protein surface file
		if (proteinFile)
		{
			ifstream proteinDescriptorFile(rootDir + "Protein_Descriptor_Vectors/FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "/" + proteinName + "DescriptorVector.txt");
			if (!proteinDescriptorFile)	//if the descriptor file of this protein doesn't exist
			{

				vector <string> surfaceResidueCoordinates;
				int lineCounter = 0;
				while (getline(proteinFile, line))	//read the lines of the protein file
					if (line.find("ATOM") == 0)
						surfaceResidueCoordinates.push_back(line);

				proteinFile.close();

				if ((int)surfaceResidueCoordinates.size() >= fragmentLength)
				{
					ofstream proteinDescriptorFile(rootDir + "Protein_Descriptor_Vectors/FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "/" + proteinName + "DescriptorVector.txt");

					if (proteinDescriptorFile)
					{
						proteinSurfaceDescriptor.clear();

						unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> fragmentHashTable;

						for (int lineIndex1 = 0; lineIndex1 <= (int)surfaceResidueCoordinates.size() - (fragmentLength); lineIndex1 = lineIndex1 + incrementValue)	//extracts continuous fragments from protein surfaces
						{

							int previousResidueNumber = atoi(surfaceResidueCoordinates[lineIndex1].substr(22, 5).c_str());
							bool continuousResidues = true;

							for (int lineIndex2 = (lineIndex1 + 1); lineIndex2 < (lineIndex1 + fragmentLength); lineIndex2++)
							{

								int nextResidueNumber = atoi(surfaceResidueCoordinates[lineIndex2].substr(22, 5).c_str());

								if (previousResidueNumber < nextResidueNumber - 1)
								{
									//cout << "\nnon continuous residues\t" << lineIndex2;
									continuousResidues = false;
									lineIndex1 = lineIndex2 - incrementValue;
									break;
								}
								previousResidueNumber = nextResidueNumber;
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

														if (selectedRS.Rx != 1000)
														{

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


															matched = CompareTwoHashTables(clusterHashTable, fragmentHashTable, expectedMatchedPoints);	//compares two hash tables
															if (matched) break;	// if the cluster hash table and the new one are similar
														}
													}
												if (matched) break;
											}// for j
											if (matched) break;
										}// for i
										if (matched) break;
									}//for c

									if (matched)
										if (proteinSurfaceDescriptor.find(scannedClusters) != proteinSurfaceDescriptor.end())	//found
											proteinSurfaceDescriptor[scannedClusters]++;

										else //not found
											proteinSurfaceDescriptor[scannedClusters] = 1;

									scannedClusters++;
								}

								clustersFile.clear();
								clustersFile.seekg(0, ios::beg);
								getline(clustersFile, line);	//To pass the line having cluster number

							}
						}

						for (auto it1 = proteinSurfaceDescriptor.begin(); it1 != proteinSurfaceDescriptor.end(); ++it1)
							proteinDescriptorFile << it1->first << ":" << it1->second << ",";

						proteinDescriptorFile << "\n";
						proteinDescriptorFile.close();

					}
					else
						cout << "\nError creating the protein surface descriptor vector file " + proteinName;
				}
				else
					cout << "\nProtein surface " + proteinName + "  is smaller than the fragment length.";
			}

			proteinFile.close();
			proteinDescriptorFile.close();
		}
		else
			cout << "\nError reading the protein surface file " + proteinName;
	}
	return 0;
}


//****************************************************************************
/* This process loads a hash table from a file */

void LoadHashTableFromFile(ifstream &clustersFile, int clusterNumber, unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> &hashTable)
{

	string line;

	hashTable.clear();

	while (getline(clustersFile, line))
	{

		if (line.find("Cluster") != 0)
		{
			istringstream iss(line);
			vector <string> lineTokens{ istream_iterator<string>{iss}, istream_iterator<string>{} };

			Point p;
			p.x = (float)stoi(lineTokens[0]);
			p.y = (float)stoi(lineTokens[1]);
			p.z = (float)stoi(lineTokens[2]);

			hashTable[p.x][p.y][p.z].push_back(lineTokens[3]);

			for (int i = 4; i < (int)lineTokens.size(); i++)
				hashTable[p.x][p.y][p.z].push_back(lineTokens[i]);
		}
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

//****************************************************************************************
/*
This process compares the descriptors of proteins in PRISM pair_list with interface descriptor vectors.
The descriptors are compressed. Parallel results can be processed.
*/

int CompareProteinWithInterfaces_v2(int fragmentLength, int overlappingResidues, int binSize, int expectedMatchedPoints, int includeRatio)
{

	map<int, int> proteinSurfaceDescriptor_1, proteinSurfaceDescriptor_2;
	map<int, int> interfaceSideDescriptorMap, interfaceSideDescriptor_1, interfaceSideDescriptor_2;
	//map<string, string> PrismInterfaceDescriptors;
	map<string, map<int, int>>PrismInterfaceDescriptors;
	int key, value;

	string line, proteinName_1, proteinName_2, interfaceSideName_1 = "noname", interfaceSideName_2 = "noname", interfaceSideName, interfaceDescriptorString;
	ifstream pairList, interfaceDescriptorsLibrary, interfacesListFile;
	vector<string> interfacesList;
	ofstream interfacesSimilarToProteins;
	vector <string> proteinList;
	vector <string> potentialSimilarInterfaces;

	//Read the target proteins pair list
	try
	{
		pairList.open("D:/PhD/prism/prism_standalone/Temp_pair_list_new");

		if (pairList.fail())
		{
			cout << "\nError reading the protein pair list file.";
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


	//Reads the interfaces list
	try
	{
		interfacesListFile.open(rootDir + "InterfaceList.txt");

		if (interfacesListFile.fail())
		{
			cout << "\nError reading the interface list file.";
			return -1;
		}
		else
		{
			while (getline(interfacesListFile, line))
			{
				//solves linux "\n" problem
				line.erase(std::remove_if(line.begin(),
					line.end(),
					[](unsigned char x) {return isspace(x); }),
					line.end());

				interfacesList.push_back(line);
				//interfacesList.push_back(line.substr(0, line.size() - 1));	//linux version
			}
		}
		interfacesListFile.close();
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


	//Create the output file to put interfaces similar to the proteins there	
	try
	{

		interfacesSimilarToProteins.open(rootDir + "Filtering_Results/FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "/" + to_string(includeRatio) + "%InterfacesSimilarToProteinPairs.txt");

		if (interfacesSimilarToProteins.fail())
		{
			cout << "\nError creating the output file.";
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


	//Reads the interface descriptors library file to PrismInterfaceDescriptors
	try
	{
		interfaceDescriptorsLibrary.open(rootDir + "PRISM_Interface_Descriptors/FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "/TempPrismInterfaceDescriptors.txt"); //All PRISM interface descriptors

		if (interfaceDescriptorsLibrary.fail())
		{
			cout << "\nError reading the interface descriptors library file.";
			return -1;
		}
		else
		{
			while (getline(interfaceDescriptorsLibrary, line))
			{


				//int tabPos = line.find('\t');
				interfaceSideName = line.substr(0, 12);
				//cout << endl << interfaceSideName;
				interfaceDescriptorString = line.substr(12, line.length() - 12);
				//cout << endl << interfaceDescriptor;

				//solves linux "\n" problem
				interfaceDescriptorString.erase(std::remove_if(interfaceDescriptorString.begin(),
					interfaceDescriptorString.end(),
					[](unsigned char x) {return isspace(x); }),
					interfaceDescriptorString.end());

				//cout << endl << interfaceDescriptor;
				interfaceSideDescriptorMap.clear();

				if ((int)interfaceDescriptorString.length() > 2)
				{
					int commaPos = -1;
					while (commaPos < (int)interfaceDescriptorString.length() - 1)	//Load interface descriptor from file
					{
						key = stoi(interfaceDescriptorString.substr(commaPos + 1, interfaceDescriptorString.find(":", commaPos + 1)));
						value = stoi(interfaceDescriptorString.substr(interfaceDescriptorString.find(":", commaPos + 1) + 1, interfaceDescriptorString.find(",", commaPos + 1)));
						//cout << "\nKey: " << key << "\tValue: " << value;
						interfaceSideDescriptorMap[key] = value;
						commaPos = interfaceDescriptorString.find(",", commaPos + 1);

					}

					//cout << interfaceSideDescriptor_1.size();
					if ((int)interfaceSideDescriptorMap.size() > 0)
						PrismInterfaceDescriptors[interfaceSideName] = interfaceSideDescriptorMap;
				}
				//cout << endl << line;
			}
		}
		interfaceDescriptorsLibrary.close();
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



	/*for (auto it1 = PrismInterfaceDescriptors.cbegin(); it1 != PrismInterfaceDescriptors.cend(); ++it1)
	{
	cout << endl << it1->first << " ";

	for (auto it2 = it1->second.cbegin(); it2 != it1->second.cend(); ++it2)
	cout << it2->first << " " << it2->second << " ";

	//cout << endl;
	getchar();
	}*/



	while (getline(pairList, line))	//each line has two protein names
	{
		int spac = line.find(" ", 0);

		proteinName_1 = line.substr(0, spac);
		cout << "\n" << proteinName_1;

		proteinName_2 = line.substr(spac + 1, line.length());
		proteinName_2.erase(remove_if(proteinName_2.begin(),
			proteinName_2.end(),
			[](unsigned char x) {return isspace(x); }),
			proteinName_2.end());

		cout << " - " << proteinName_2;


		ifstream proteinDescriptorFile_1(rootDir + "Protein_Descriptor_Vectors/FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "/" + proteinName_1 + "DescriptorVector.txt");
		ifstream proteinDescriptorFile_2(rootDir + "Protein_Descriptor_Vectors/FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "/" + proteinName_2 + "DescriptorVector.txt");
		if (proteinDescriptorFile_1 && proteinDescriptorFile_2)
		{

			interfacesSimilarToProteins << "Interfaces similar to " << proteinName_1 << "-" << proteinName_2 << ":\n";

			proteinSurfaceDescriptor_1.clear();
			proteinSurfaceDescriptor_2.clear();

			//Load protein_1 surface descriptor from file
			//getline(proteinDescriptorFile_1, line);
			//int length = line.length();
			//int commaPos = -1;

			while (getline(proteinDescriptorFile_1, line, ','))
			{
				spac = line.find(":", 0);
				if (spac>0)
				{
					key = stoi(line.substr(0, spac));
					value = stoi(line.substr(spac + 1, line.length()));
					proteinSurfaceDescriptor_1[key] = value;
				}


			}
			/*while (commaPos < length - 1)
			{
			key = stoi(line.substr(commaPos + 1, line.find(":", commaPos + 1)));
			value = stoi(line.substr(line.find(":", commaPos + 1) + 1, line.find(",", commaPos + 1)));
			//cout << "\nKey: " << key << "\tValue: " << value;
			proteinSurfaceDescriptor_1[key] = value;
			commaPos = line.find(",", commaPos + 1);
			}*/

			/*cout << endl;
			for (auto it1 = proteinSurfaceDescriptor_1.cbegin(); it1 != proteinSurfaceDescriptor_1.cend(); ++it1)
			{
			cout << it1->first << " " << it1->second << " ";
			}
			getchar();*/

			//Load protein_2 surface descriptor from file
			//getline(proteinDescriptorFile_2, line);
			//length = line.length();
			//commaPos = -1;

			while (getline(proteinDescriptorFile_2, line, ','))
			{
				spac = line.find(":", 0);
				if (spac > 0)
				{
					key = stoi(line.substr(0, spac));
					value = stoi(line.substr(spac + 1, line.length()));
					proteinSurfaceDescriptor_2[key] = value;
				}
				//cout << line<<endl;
			}

			/*while (commaPos < length - 1)
			{
			key = stoi(line.substr(commaPos + 1, line.find(":", commaPos + 1)));
			value = stoi(line.substr(line.find(":", commaPos + 1) + 1, line.find(",", commaPos + 1)));
			//cout << "\nKey: " << key << "\tValue: " << value;
			proteinSurfaceDescriptor_2[key] = value;
			commaPos = line.find(",", commaPos + 1);
			}*/

			/*cout << endl;
			for (auto it1 = proteinSurfaceDescriptor_2.cbegin(); it1 != proteinSurfaceDescriptor_2.cend(); ++it1)
			{
			cout << it1->first << " " << it1->second << " ";
			}
			getchar();*/




			int length, commaPos;
			bool side1Existence = false;


			for (int i = 0; i < interfacesList.size(); )
			{

				interfaceSideName_1 = interfacesList[i];

				if (i + 1 < interfacesList.size())
				{
					interfaceSideName_2 = interfacesList[i + 1];
				}

				interfaceSideName_1 = interfaceSideName_1.substr(0, interfaceSideName_1.find('_'));
				interfaceSideName_2 = interfaceSideName_2.substr(0, interfaceSideName_2.find('_'));
				//cout << endl << interfaceSideName_1;
				//cout << endl << interfaceSideName_2;

				if (interfaceSideName_1 == interfaceSideName_2) //&& interfaceSideName_1!=""
				{
					interfaceSideName_1 = interfacesList[i];
					interfaceSideName_2 = interfacesList[i + 1];
					//cout << endl << interfaceSideName_1;
					//cout << endl << interfaceSideName_2;

					if ((int)PrismInterfaceDescriptors[interfaceSideName_1].size() > 0 && (int)PrismInterfaceDescriptors[interfaceSideName_2].size() > 0)
					{
						//cout << endl << interfaceSideName_1;
						//cout << endl << interfaceSideName_2;
						//cout << endl << "They will be compared.";

						if (CompareProteinSurfaceWithInterfaceSide(proteinSurfaceDescriptor_1, PrismInterfaceDescriptors[interfaceSideName_1], (float)includeRatio) && CompareProteinSurfaceWithInterfaceSide(proteinSurfaceDescriptor_2, PrismInterfaceDescriptors[interfaceSideName_2], (float)includeRatio))
						{
							interfacesSimilarToProteins << interfaceSideName_1.substr(0, interfaceSideName_1.find('_')) << " \n";
							//cout << endl << "Similar";
						}
						else if (CompareProteinSurfaceWithInterfaceSide(proteinSurfaceDescriptor_1, PrismInterfaceDescriptors[interfaceSideName_2], (float)includeRatio) && CompareProteinSurfaceWithInterfaceSide(proteinSurfaceDescriptor_2, PrismInterfaceDescriptors[interfaceSideName_1], (float)includeRatio))
						{
							interfacesSimilarToProteins << interfaceSideName_1.substr(0, interfaceSideName_1.find('_')) << " \n";
							//cout << endl << "Similar";
						}
						//getchar();
					}

					i += 2;
				}
				else
				{
					i++;
				}

			}



			//-------------------------------------------------------------------------------------------
			/*for (auto& kv : PrismInterfaceDescriptors) {
			std::cout << kv.first << " has value " << kv.second << std::endl;
			if (interfaceSideName_1 == "")
			{

			length = line.length();

			if (length > 13)
			{

			side1Existence = true;
			interfaceSideDescriptor_1.clear();
			interfaceSideName_1 = line.substr(0, 12);
			//cout << "\n" << "interfaceSideName_1: " << interfaceSideName_1;
			//getchar();



			while (commaPos < length - 1)	//Load interface descriptor from file
			{
			key = stoi(line.substr(commaPos + 1, line.find(":", commaPos + 1)));
			value = stoi(line.substr(line.find(":", commaPos + 1) + 1, line.find(",", commaPos + 1)));
			//cout << "\nKey: " << key << "\tValue: " << value;
			interfaceSideDescriptor_1[key] = value;
			commaPos = line.find(",", commaPos + 1);
			}

			}
			else
			side1Existence = false;

			//getchar();
			}

			}/////







			//Compare each protein pair with all interfaces
			/*int lineNumber = 0;


			//while (getline(interfaceDescriptorsLibrary, line)) //each line has an interface descriptor
			while (getline(interfaceDescriptorsLibrary, line)) //each line has an interface descriptor
			{

			//lineNumber++;

			//if (lineNumber % 2 != 0)	//Read side 1 of the interface
			if (interfaceSideName_1 == "noname")
			{

			length = line.length();

			if (length > 13)
			{

			side1Existence = true;
			interfaceSideDescriptor_1.clear();
			interfaceSideName_1 = line.substr(0, 12);
			//cout << "\n" << "interfaceSideName_1: " << interfaceSideName_1;
			//getchar();

			commaPos = 12;

			while (commaPos < length - 1)	//Load interface descriptor from file
			{
			key = stoi(line.substr(commaPos + 1, line.find(":", commaPos + 1)));
			value = stoi(line.substr(line.find(":", commaPos + 1) + 1, line.find(",", commaPos + 1)));
			//cout << "\nKey: " << key << "\tValue: " << value;
			interfaceSideDescriptor_1[key] = value;
			commaPos = line.find(",", commaPos + 1);
			}

			}
			else
			side1Existence = false;

			//getchar();
			}

			else if (interfaceSideName_2 == "noname") //Read side 2 of the interface and do the comparisons
			{

			length = line.length();

			if (length > 13 && side1Existence)
			{

			interfaceSideName_2 = line.substr(0, 12);
			//cout << "\n" << "interfaceSideName_2: " << interfaceSideName_2;
			//getchar();


			if (interfaceSideName_1.substr(0, 6) == interfaceSideName_2.substr(0, 6))
			{

			interfaceSideDescriptor_2.clear();
			commaPos = 12;

			while (commaPos < length - 1)	//Load interface descriptor from file
			{
			key = stoi(line.substr(commaPos + 1, line.find(":", commaPos + 1)));
			value = stoi(line.substr(line.find(":", commaPos + 1) + 1, line.find(",", commaPos + 1)));
			//cout << "\nKey: " << key << "\tValue: " << value;
			interfaceSideDescriptor_2[key] = value;
			commaPos = line.find(",", commaPos + 1);
			}
			//getchar();


			bool interfaceSurfaceSimilarity = CompareProteinSurfaceWithInterfaceSide(proteinSurfaceDescriptor_1, interfaceSideDescriptor_1, (float)includeRatio);
			//cout << endl << interfaceSurfaceSimilarity;
			if (interfaceSurfaceSimilarity)
			{
			//cout << endl << "Prtein A is similar to interface side 1.";
			interfaceSurfaceSimilarity = CompareProteinSurfaceWithInterfaceSide(proteinSurfaceDescriptor_2, interfaceSideDescriptor_2, (float)includeRatio);
			if (interfaceSurfaceSimilarity)
			{
			interfacesSimilarToProteins << interfaceSideName_1.substr(0, 6) << "\n";
			//cout << endl << "Prtein B is similar to interface side 2.";
			}
			}

			//else
			if (!interfaceSurfaceSimilarity)
			{
			interfaceSurfaceSimilarity = CompareProteinSurfaceWithInterfaceSide(proteinSurfaceDescriptor_1, interfaceSideDescriptor_2, (float)includeRatio);
			if (interfaceSurfaceSimilarity)
			{
			//cout << endl << "Prtein A is similar to interface side 2.";
			interfaceSurfaceSimilarity = CompareProteinSurfaceWithInterfaceSide(proteinSurfaceDescriptor_2, interfaceSideDescriptor_1, (float)includeRatio);
			if (interfaceSurfaceSimilarity)
			{
			interfacesSimilarToProteins << interfaceSideName_1.substr(0, 6) << "\n";
			//cout << endl << "Prtein B is similar to interface side 1.";
			}
			}
			}
			//cout << endl << interfaceSurfaceSimilarity;
			//getchar();
			interfaceSideName_1 = "noname";
			}
			else
			interfaceSideName_1 = interfaceSideName_2;

			interfaceSideName_2 = "noname";

			}

			}
			}////

			//interfaceDescriptorsLibrary.clear();
			//interfaceDescriptorsLibrary.seekg(0, ios::beg);*/
			//-------------------------------------------------------------------------------------------
		}
		else
			cout << "\nError reading the protein surface descriptor file for protein: " + proteinName_1 + "or " + proteinName_2;
	}


	interfacesSimilarToProteins.close();
	pairList.close();

	return 0;
}

//****************************************************************************************
bool CompareProteinSurfaceWithInterfaceSide(map<int, int> proteinSurfaceDescriptor, std::map<int, int> interfaceSideDescriptor, float includeRatio)
{




	cout << endl;
	for (auto it1 = proteinSurfaceDescriptor.cbegin(); it1 != proteinSurfaceDescriptor.cend(); ++it1)
	{
		cout << it1->first << " " << it1->second << " ";
	}
	//getchar();

	cout << endl;
	for (auto it1 = interfaceSideDescriptor.cbegin(); it1 != interfaceSideDescriptor.cend(); ++it1)
	{
		cout << it1->first << " " << it1->second << " ";
	}
	getchar();



	int includedFragmentsCount = 0;
	int interfaceFragmentsCount = 0;
	std::map<int, int>::iterator proteinIterator = proteinSurfaceDescriptor.begin();

	for (auto interfaceIterator = interfaceSideDescriptor.begin(); interfaceIterator != interfaceSideDescriptor.end(); ++interfaceIterator)
	{
		//cout << endl << interfaceIterator->first;
		//getchar();
		interfaceFragmentsCount += interfaceIterator->second;

		while (proteinIterator != proteinSurfaceDescriptor.end() && proteinIterator->first < interfaceIterator->first)
		{

			//cout << endl << proteinIterator->first;
			//getchar();
			proteinIterator++;
		}

		if (proteinIterator != proteinSurfaceDescriptor.end())
		{
			if (proteinIterator->first == interfaceIterator->first)
			{
				if (proteinIterator->second >= interfaceIterator->second)
					includedFragmentsCount += interfaceIterator->second;
				else
					includedFragmentsCount += interfaceIterator->second - proteinIterator->second;

				proteinIterator++;
			}
		}
	}

	//cout << endl << interfaceFragmentsCount;
	//cout << endl << includedFragmentsCount;

	float interfaceSurfaceSimilarityRate = float(includedFragmentsCount) / float(interfaceFragmentsCount);
	//cout << endl << interfaceSurfaceSimilarityRate;

	if (interfaceSurfaceSimilarityRate >= includeRatio / 100)
		return true;
	else
		return false;
	//	potentialSimilarInterfaces.push_back((interfaceName.substr(4, 5) == interfaceName.substr(7, 8)) ? interfaceName.substr(5, 6) : interfaceName.substr(4, 5));
	//add this interface to potential similar interfaces list to be compared to protein_2
	//}
}