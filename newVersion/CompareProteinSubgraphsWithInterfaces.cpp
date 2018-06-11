/*
This code ...

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
This process ...
*/

int CompareProteinPairsWithInterfaces(string jobName, int numberOfClusters, int fragmentLength, int overlappingResidues, int expectedMatchedPoints, int binSize)
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

