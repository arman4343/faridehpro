/*
This code creates the protein surface descriptor vectors by extracting the small protein fragments and finding their
similarities to the fragment clusters. Then it compares the proteins descriptor vectors with interface descriptor vector
library to find which interface can be on these proteins surfaces.

Written by Farideh Halakou
*/

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
#include <ctype.h>
#include "My_Functions.h"


#define PI 3.14159265
#pragma warning(disable : 4503) 
using namespace std;


//****************************************************************************
/*
This process creates descriptor vectors for protein in PRISM pair_list from which the surfaces are extracted.
It creates the descriptor vector as a map just showing the frequency of the existent fragments, not as a sparce vector.
*/

int CreateProteinDescriptors_v2(string jobName, int fragmentLength, int overlappingResidues, int binSize, int expectedMatchedPoints, int numberOfClusters)
{

	std::map<int, int> proteinSurfaceDescriptor;

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
	
	pairList.close();


	for (int proteinIndex = 0; proteinIndex < (int)proteinList.size(); ++proteinIndex)
	{

		proteinName = proteinList[proteinIndex];

		ifstream proteinFile(("D:\\PhD\\prism\\prism_standalone\\jobs\\" + jobName + "\\surfaceExtract\\" + proteinName + ".asa.pdb").c_str());	//read the protein surface file
		if (proteinFile)
		{
			ifstream proteinDescriptorFile("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Protein_Descriptor_Vectors\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\" + proteinName + "DescriptorVector.txt");
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
					ofstream proteinDescriptorFile("D:\\PhD\\My_Thesis\\Second_Step\\Data\\Protein_Descriptor_Vectors\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\" + proteinName + "DescriptorVector.txt");

					if (proteinDescriptorFile)
					{
						proteinSurfaceDescriptor.clear();

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

														if (selectedRS.Rx != 0 || selectedRS.Ry != 0 || selectedRS.Rz != 0)
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


//****************************************************************************************
/*
This process compares the descriptors of proteins in PRISM pair_list with interface descriptor vectors.
The descriptors are compressed.
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

//****************************************************************************************
/*
This process compares the descriptors of proteins in PRISM pair_list with interface descriptor vectors.
The descriptors are compressed.
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
		pairList.open(rootDir + "pair_list");

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

		interfacesSimilarToProteins.open(rootDir + "Parallel_Filtering_Results" + slash + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + slash + to_string(includeRatio) + "%InterfacesSimilarToProteinPairs.txt");

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
		interfaceDescriptorsLibrary.open(rootDir + "Parallel_PRISM_Interface_Descriptors" + slash + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + slash + "PrismInterfaceDescriptors.txt"); //All PRISM interface descriptors

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
				interfaceDescriptorString = line.substr(12, line.length()-12);
				//cout << endl << interfaceDescriptor;

				//solves linux "\n" problem
				interfaceDescriptorString.erase(std::remove_if(interfaceDescriptorString.begin(),
					interfaceDescriptorString.end(),
					[](unsigned char x) {return isspace(x); }),
					interfaceDescriptorString.end());

				//cout << endl << interfaceDescriptor;
				interfaceSideDescriptor_1.clear();
				
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
		cout << " - " << proteinName_2;


		ifstream proteinDescriptorFile_1(rootDir + "Parallel_Protein_Descriptor_Vectors" + slash + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + slash + proteinName_1 + "DescriptorVector.txt");
		ifstream proteinDescriptorFile_2(rootDir + "Parallel_Protein_Descriptor_Vectors" + slash + "FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + slash + proteinName_2 + "DescriptorVector.txt");
		if (proteinDescriptorFile_1 && proteinDescriptorFile_2)
		{

			interfacesSimilarToProteins << "Interfaces similar to " << proteinName_1 << "-" << proteinName_2 << ":\n";

			proteinSurfaceDescriptor_1.clear();
			proteinSurfaceDescriptor_2.clear();

			//Load protein_1 surface descriptor from file
			//getline(proteinDescriptorFile_1, line);
			//int length = line.length();
			//int commaPos = -1;

			while (getline(proteinDescriptorFile_1, line,','))
			{
				spac = line.find(":", 0);
				if (spac>0)
				{
					key= stoi(line.substr(0, spac));
					value = stoi(line.substr( spac + 1, line.length() ));
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
			getchar();
			*/



			int length, commaPos;
			bool side1Existence = false;

			
			for (int i = 0; i < interfacesList.size(); ) 
			{

				interfaceSideName_1 = interfacesList[i];

				if (i + 1 < interfacesList.size())
				{
					interfaceSideName_2 = interfacesList[i+1];
				}

				interfaceSideName_1 = interfaceSideName_1.substr(0, interfaceSideName_1.find('_'));
				interfaceSideName_2 = interfaceSideName_2.substr(0, interfaceSideName_2.find('_'));
				//cout << endl << interfaceSideName_1;
				//cout << endl << interfaceSideName_2;

				if (interfaceSideName_1== interfaceSideName_2) //&& interfaceSideName_1!=""
				{
					interfaceSideName_1 = interfacesList[i];
					interfaceSideName_2 = interfacesList[i + 1];
					//cout << endl << interfaceSideName_1;
					//cout << endl << interfaceSideName_2;

					if ((int)PrismInterfaceDescriptors[interfaceSideName_1].size() > 0 && (int)PrismInterfaceDescriptors[interfaceSideName_2].size() > 0 )
					{
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

//----------------------------------------------------------------
bool CompareProteinSurfaceWithInterfaceSide(std::map<int, int> proteinSurfaceDescriptor, std::map<int, int> interfaceSideDescriptor, float includeRatio)
{

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

	float interfaceSurfaceSimilarityRate = float (includedFragmentsCount) / float(interfaceFragmentsCount);
	//cout << endl << interfaceSurfaceSimilarityRate;

	if (interfaceSurfaceSimilarityRate >= includeRatio / 100)
		return true;
	else
		return false;
	//	potentialSimilarInterfaces.push_back((interfaceName.substr(4, 5) == interfaceName.substr(7, 8)) ? interfaceName.substr(5, 6) : interfaceName.substr(4, 5));
	//add this interface to potential similar interfaces list to be compared to protein_2
	//}
}