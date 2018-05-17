#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
//#include <#cmath>

using namespace std;

//****************************************************************************
/*
This process calculates the cosine distance between the PRISM interface side descriptors.
*/

int CalculateCosineDistance(int fragmentLength, int overlappingResidues, int binSize, int expectedMatchedPoints, int numberOfInterfaces, int numberOfClusters)
{

	
	ifstream descriptorsFile;
	ofstream descriptorsDistances;
	string line, interfaceName;
	int numberOfDescriptors=0, lineLen = 0;
	double averageDistance, sumOfAverageDistances = 0;
	vector <string> fileContent, interfaceNames;
	vector<vector<int>> interfaceDescriptors;
	string dataPath = "D:\\PhD\\My_Thesis\\Second_Step\\Data\\";

	//Read the PRISM interface descriptors file
	try
	{
		descriptorsFile.open(dataPath + "PRISM_Interface_Descriptors\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\PrismInterfaceDescriptors.txt");

		if (descriptorsFile.fail())
		{
			cout << "\nError reading the interface descriptors file.";
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


	//Create the output file to put the interface descriptors average distances	
	try
	{
		descriptorsDistances.open(dataPath + "PRISM_Interface_Descriptors\\FrLn" + to_string(fragmentLength) + "_OvRd" + to_string(overlappingResidues) + "_MtPn" + to_string(expectedMatchedPoints) + "_BnSz" + to_string(binSize) + "\\DescriptorsCosineDistances.txt");

		if (descriptorsDistances.fail())
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


	descriptorsDistances << "Interface_Name\tAverage_Distance\tMost_Similar_Interface_Name\tMinimum_Distance\tLeast_Similar_Interface_Name\tMaximum_Distance\n";


	while (getline(descriptorsFile, line))
		fileContent.push_back(line);
	descriptorsFile.close();


	for (int i=0; i < numberOfInterfaces; i++)
	{
		line = fileContent[i];
		bool descriptorExists = false;
		interfaceName = line.substr(0, 12);

		if (line.length() > 13)	//if the descriptor exists
		{

			descriptorExists = true;
			interfaceNames.push_back(interfaceName);
			cout << endl << interfaceName;
			vector<int> interfaceDescriptor;
			for (int t = 0; t < numberOfClusters; t++)
				interfaceDescriptor.push_back(0);


			//Load interface descriptor
			line = line.substr(13, line.length()-13);
			//cout << endl << line;
			int vectorIndex = 0;
			int commaPos = -1;
			lineLen = line.length();
			//cout << endl << lineLen;
			//getchar();
			while (commaPos < lineLen - 1)
			{
				int key = stoi(line.substr(commaPos + 1, line.find(":", commaPos + 1)));
				int value = stoi(line.substr(line.find(":", commaPos + 1) + 1, line.find(",", commaPos + 1)));
				//cout << "\nKey: " << key << "\tValue: " << value;
				//getchar();
				while (vectorIndex < key - 1)
					vectorIndex++;

				interfaceDescriptor[vectorIndex] = value;
				vectorIndex++;
				commaPos = line.find(",", commaPos + 1);
			}

			interfaceDescriptors.push_back(interfaceDescriptor);

			numberOfDescriptors += 1;
			double maxDistance = -1000;
			double minDistance = 1000;
			string mostSimilar = "-";
			string leastSimilar = "-";
			double sumOfDistances = 0.0;

			for (int j = i + 1; j < numberOfInterfaces; j++)
			{
				line = fileContent[j];
				descriptorExists = false;
				interfaceName = line.substr(0, 12);

				if (line.length() > 13)	//if the descriptor exists
				{

					descriptorExists = true;
					interfaceNames.push_back(interfaceName);
					//cout << endl << interfaceName;
					interfaceDescriptor.clear();
					for (int t = 0; t < numberOfClusters; t++)
						interfaceDescriptor.push_back(0);


					//Load interface descriptor
					line = line.substr(13, line.length() - 13);
					//cout << endl << line;
					int vectorIndex = 0;
					int commaPos = -1;
					lineLen = line.length();
					//cout << endl << lineLen;
					//getchar();
					while (commaPos < lineLen - 1)
					{
						int key = stoi(line.substr(commaPos + 1, line.find(":", commaPos + 1)));
						int value = stoi(line.substr(line.find(":", commaPos + 1) + 1, line.find(",", commaPos + 1)));
						//cout << "\nKey: " << key << "\tValue: " << value;
						//getchar();
						while (vectorIndex < key - 1)
							vectorIndex++;

						interfaceDescriptor[vectorIndex] = value;
						vectorIndex++;
						commaPos = line.find(",", commaPos + 1);
					}

					//cout << endl;
					//for (int e = 0; e < interfaceDescriptor.size(); e++)
					//	cout << interfaceDescriptor[e] << "\t";

					interfaceDescriptors.push_back(interfaceDescriptor);
					numberOfDescriptors += 1;
			

					//Calculate cosine distance
					double mul = 0.0, d_a = 0.0, d_b = 0.0;
					for (int k = 0; k < numberOfClusters; ++k)
					{
						mul += interfaceDescriptors[i][k] * interfaceDescriptor[k];
						d_a += interfaceDescriptors[i][k] * interfaceDescriptors[i][k];
						d_b += interfaceDescriptor[k] * interfaceDescriptor[k];
					}

					double cosineDistance = 1 - (mul / (sqrt(d_a) * sqrt(d_b)));
					sumOfDistances += cosineDistance;
					//cout << endl << sumOfDistances;
					//getchar();

					if (cosineDistance > maxDistance)
					{
						maxDistance = cosineDistance;
						leastSimilar = interfaceName;
					}
					else if (cosineDistance < minDistance)
					{
						minDistance = cosineDistance;
						mostSimilar = interfaceName;
					}
				}
			}

			cout << endl << sumOfDistances << endl << numberOfDescriptors;
			averageDistance = sumOfDistances / double(numberOfDescriptors - 1);
			descriptorsDistances << interfaceNames[i] << "\t" << averageDistance << "\t" << mostSimilar << "\t" << minDistance << "\t" << leastSimilar << "\t" << maxDistance << "\n";
			sumOfAverageDistances += averageDistance;
			break;
		}
	}

	//----------------------------------------------------------------------

	for (int i=1; i < numberOfDescriptors; i++)
	{
		
		double maxDistance = -1000;
		double minDistance = 1000;
		string mostSimilar = "-";
		string leastSimilar = "-";
		double sumOfDistances = 0;

		for (int j = 0; j < numberOfDescriptors; j++)
		{
			if (i != j)
			{

				//Calculate cosine distance

				double mul = 0.0, d_a = 0.0, d_b = 0.0;

				for (int k = 0; k < numberOfClusters; ++k)
				{
					mul += interfaceDescriptors[i][k] * interfaceDescriptors[j][k];
					d_a += interfaceDescriptors[i][k] * interfaceDescriptors[i][k];
					d_b += interfaceDescriptors[j][k] * interfaceDescriptors[j][k];
				}

				double cosineDistance = 1 - (mul / (sqrt(d_a) * sqrt(d_b)));
				sumOfDistances += cosineDistance;

				if (cosineDistance > maxDistance)
				{
					maxDistance = cosineDistance;
					leastSimilar = interfaceNames[j];
				}
				else if (cosineDistance < minDistance)
				{
					minDistance = cosineDistance;
					mostSimilar = interfaceNames[j];
				}

			}
		}

		averageDistance = sumOfDistances / double(numberOfDescriptors - 1);
		sumOfAverageDistances += averageDistance;
		descriptorsDistances << interfaceNames[i] << "\t" << averageDistance << "\t" << mostSimilar << "\t" << minDistance << "\t" << leastSimilar << "\t" << maxDistance << "\n";
	}
	
	
	descriptorsDistances << "Average distance between descriptors is: " << sumOfAverageDistances / double(numberOfDescriptors);
	descriptorsDistances.close();
	return 0;
}
