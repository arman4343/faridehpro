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
This process creates a graph of protein surface for all the proteins in the target list. Then finds the 
subgraphs in each protein graph and creates the descriptors of each subgraph. Then it compares the
subgraphs' descriptors with interface descriptirs.

*/

int CompareProteinsWithInterfaces_v3(string jobName, int fragmentLength, int overlappingResidues, int expectedMatchedPoints, int binSize)
{

	vector <FragmentInfo> proteinSurfaceDescriptor;

	string line, proteinName, interfaceSideName, interfaceDescriptorString;
	map <int, int> interfaceSideDescriptorMap;
	map <string, map<int, int>> PrismInterfaceDescriptors;

	ifstream proteinsListFile, interfacesListFile, interfaceDescriptorsLibraryFile;
	ofstream interfacesSimilarToProteinsFile;
	vector <string> proteinsList, interfacesList;
	string prismPath = "D:/PhD/prism/prism_standalone/";

	FILE * clustersFile;
	char * clustersBuffer;
	int clustersFileSize, key, value;

	vector<unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>>> fragmentHashTables;
	unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> fHashTable;


	//Read the target proteins list file
	try
	{
		proteinsListFile.open(rootDir + "Target_Proteins_List.txt"); //list of target proteins

		if (proteinsListFile.fail())
		{
			cout << "\nError reading the pair list file.";
			return -1;
		}
		else
		{
			while (getline(proteinsListFile, line))
			{
				//solves linux "\n" problem
				line.erase(std::remove_if(line.begin(),
					line.end(),
					[](unsigned char x) {return isspace(x); }),
					line.end());

				proteinsList.push_back(line);
				//interfacesList.push_back(line.substr(0, line.size() - 1));	//linux version
			}
		}
		proteinsListFile.close();

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


	//Reads the fragment clusters file
	try
	{
		clustersFile = fopen((rootDir + "Clustering_Results/FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "/Clusters.txt").c_str(), "rb");

		fseek(clustersFile, 0, SEEK_END);
		clustersFileSize = ftell(clustersFile);
		clustersBuffer = (char *)malloc(clustersFileSize + 1);
		clustersBuffer[clustersFileSize] = 0;
		fseek(clustersFile, 0, SEEK_SET);
		int freadsize = 0;
		//freadsize += fread(clustersBuffer + freadsize, 1, clustersFileSize - freadsize, clustersFile);
		freadsize = fread(clustersBuffer, 1, clustersFileSize, clustersFile);
		if (clustersFileSize != freadsize)
		{
			cout << "\nError reading all of the cluster hash tables file.";
			return -1;
		}
		fclose(clustersFile);

		stringstream ss(clustersBuffer);
		string clusterLine;

		while (getline(ss, clusterLine, '\n'))
		{
			if (clusterLine.find("Cluster") == 0)	//if this is a new cluster
			{
				if (fHashTable.size()>0)
				{
					//Add the cluster hash table to fragmentHashTables vector
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
		//cout << "\nfragmentHashTables size: " << fragmentHashTables.size();
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


	//Reads the interfaces list
	/*try
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
	*/

	//Create the output file to put interfaces similar to the proteins there	
	try
	{

		interfacesSimilarToProteinsFile.open(rootDir + "Filtering_Results/FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "/" + "InterfacesSimilarToProteins.txt");

		if (interfacesSimilarToProteinsFile.fail())
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
	/*try
	{
		interfaceDescriptorsLibraryFile.open(rootDir + "PRISM_Interface_Descriptors/FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "/PrismInterfaceDescriptors.txt");

		if (interfaceDescriptorsLibraryFile.fail())
		{
			cout << "\nError reading the interface descriptors library file.";
			return -1;
		}
		else
		{
			while (getline(interfaceDescriptorsLibraryFile, line))
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
		interfaceDescriptorsLibraryFile.close();
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
	*/



	for (int proteinIndex = 0; proteinIndex < (int)proteinsList.size(); ++proteinIndex)
	{

		proteinName = proteinsList[proteinIndex];
		cout << endl << proteinName;

		ifstream proteinFile((prismPath + "jobs/" + jobName + "/surfaceExtract/" + proteinName + ".asa.pdb").c_str());	//read the protein surface file
		if (proteinFile)
		{
			
			vector <string> surfaceResidueCoordinates;
			int lineCounter = 0;
			while (getline(proteinFile, line))	//read the lines of the protein file
				if (line.find("ATOM") == 0)
					surfaceResidueCoordinates.push_back(line);

			proteinFile.close();

			if ((int)surfaceResidueCoordinates.size() >= fragmentLength)
			{
				ofstream proteinDescriptorFile(rootDir + "Protein_Descriptor_Vectors/FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "/" + proteinName + "_DescriptorVector.txt");

				if (proteinDescriptorFile)
				{

					CreateProteinSurfaceGraph(proteinSurfaceDescriptor, surfaceResidueCoordinates, fragmentHashTables, fragmentLength, overlappingResidues, expectedMatchedPoints, binSize);
					for (auto it1 = 0; it1 < proteinSurfaceDescriptor.size(); ++it1)
						proteinDescriptorFile << proteinSurfaceDescriptor[it1].fragmentNo << ", " << proteinSurfaceDescriptor[it1].clusterNo << ", [" << proteinSurfaceDescriptor[it1].fragmentCenter.x << "," << proteinSurfaceDescriptor[it1].fragmentCenter.y << "," << proteinSurfaceDescriptor[it1].fragmentCenter.z << "]\n";
		
					proteinDescriptorFile.close();

				}
				else
					cout << "\nError creating the protein surface descriptor vector file " + proteinName;
			}
			else
				cout << "\nProtein surface " + proteinName + "  is smaller than the fragment length.";

			proteinFile.close();
		}
		else
			cout << "\nError reading the protein surface file " + proteinName;
	}
	return 0;
}



//****************************************************************************************
/*
This process extracts the continuous fragments from the protein surface and create a graph of 
fragments based on the fagment center distances.
*/

void CreateProteinSurfaceGraph(vector <FragmentInfo>& proteinSurfaceDescriptor, vector <string>& surfaceResidueCoordinates, vector<unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>>>& fragmentHashTables, int fragmentLength, int overlappingResidues, int expectedMatchedPoints, int binSize)
{

	proteinSurfaceDescriptor.clear();
	int fragmentNumber = 0;	//Fragment number on protein surface
	int incrementValue = fragmentLength - overlappingResidues;

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


			//for (int temp = lineIndex1; temp < (lineIndex1 + fragmentLength); temp++)
			//	cout << endl << surfaceResidueCoordinates[temp];
			//getchar();

			unordered_map<float, unordered_map<float, unordered_map<float, vector <string>>>> clusterHashTable;
			bool matched = false;
			int scannedClusters = 1;

			while (!matched && scannedClusters <= fragmentHashTables.size()) //till we find which cluster hash table is similar to the fragment hashtable or none of them are similar
			{
				
				clusterHashTable.clear();
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
				{
					
					FragmentInfo fragInfo;

					fragInfo.fragmentNo = fragmentNumber;
					fragInfo.clusterNo = scannedClusters;

					//Calculate fragment's center
					float centerX = 0, centerY = 0, centerZ = 0;
					for (int fragmentIterator = lineIndex1; fragmentIterator < (lineIndex1 + fragmentLength); fragmentIterator++)
					{
						centerX += (float)atof(surfaceResidueCoordinates[fragmentIterator].substr(30, 8).c_str());
						centerY += (float)atof(surfaceResidueCoordinates[fragmentIterator].substr(38, 8).c_str());
						centerZ += (float)atof(surfaceResidueCoordinates[fragmentIterator].substr(46, 8).c_str());
					}

					fragInfo.fragmentCenter.x = (float)centerX / fragmentLength;
					fragInfo.fragmentCenter.y = (float)centerY / fragmentLength;
					fragInfo.fragmentCenter.z = (float)centerZ / fragmentLength;

					proteinSurfaceDescriptor.push_back(fragInfo);
				}

				scannedClusters++;
			}

			fragmentNumber++;
		}
	}

}