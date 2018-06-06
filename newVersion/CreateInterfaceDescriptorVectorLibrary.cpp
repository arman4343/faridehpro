/*
This code creates the library of protein interface descriptors by using the protein small
fragments library. It actually converts the template interface files of PRISM to their descriptors.
Written by Farideh Halakou
*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <map>
#include <unordered_map>
#include <list>
#include "My_Functions.h"


#define PI 3.14159265

using namespace std;


//****************************************************************************************
/*
This function creates the descriptor vectors of PRISM interfaces as compact vectors showing
the frequency of each fragment in a paralle manner.
*/

int CreateInterfaceDescriptors_v2_Parallel(int fragmentLength, int overlappingResidues, int binSize, int expectedMatchedPoints, int numberOfClusters)
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
		else {
			while (getline(interfacesListFile, line)) {
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
		interfaceDescriptorsLibraryFile.open(rootDir + "Parallel_PRISM_Interface_Descriptors" + slash + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + slash + "PrismInterfaceDescriptors.txt");

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
		clustersFile = fopen((rootDir + "Parallel_Clustering_Results" + slash + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + slash + "Clusters.txt").c_str(), "rb");


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

		stringstream ss(clustersBuffer);
		string clusterLine;

		while (std::getline(ss, clusterLine, '\n'))
		{
			if (clusterLine.find("Cluster") == 0)	//if this is a new cluster
			{
				if (fHashTable.size()>0)
				{
					fragmentHashTables.push_back(fHashTable);
					fHashTable.clear();
				}
			}
			else 
			{
				//Read cluster hash table from file
				istringstream iss(clusterLine);
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
		cout << "fragmentHashTables size: " << fragmentHashTables.size();
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
			//interfaceFileName = interfaceFileName.substr(0, interfaceFileName.size() - 1);

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

							//clustersFile.clear();
							//clustersFile.seekg(0, ios::beg);
							//getline(clustersFile, line);	//To pass the line having cluster number
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







