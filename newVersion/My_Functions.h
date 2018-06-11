#ifndef _MY_FUNCTIONS_
#define _MY_FUNCTIONS_

#include <map>
#include <list>

struct Point
{
	float x = 0, y = 0, z = 0;
};

struct TranslationParameter {
	Point Transfer;
	float Ry = 0, Rz = 0, Rx = 0;
	Point p2, p3;
};

extern std::string rootDir;
extern std::string slash;


int main(void);
int ExtractFragments_v1(int, int);
int Extract_NonContinuous_Fragments(int, int);
int ReadNumberOfFragments(int, int);
int ClusterFragments_v1_Parallel(int, int, int, int, int, int);
TranslationParameter CalculateGeoTranslation(float, float, float, float, float, float, float, float, float, int);
Point CalculateNewPoint(TranslationParameter, Point, int);
void AddToHashTable(std::unordered_map<float, std::unordered_map<float, std::unordered_map<float, std::vector <std::string>>>>&, Point, int, int, int);
bool CompareTwoHashTables(std::unordered_map<float, std::unordered_map<float, std::unordered_map<float, std::vector <std::string>>>>, std::unordered_map<float, std::unordered_map<float, std::unordered_map<float, std::vector <std::string>>>>, int);
int WriteHashTableToFile(std::string, int, std::unordered_map<float, std::unordered_map<float, std::unordered_map<float, std::vector <std::string>>>>);
int CreateInterfaceDescriptors_v2_Parallel(int, int, int, int);
int CompareProteinPairsWithInterfaces(std::string, int, int, int, int, int);



int CreateProteinDescriptors(std::string, int, int, int, int, int);
int CreateProteinDescriptors_v2(std::string, int, int, int, int, int);
int CreateProteinDescriptors_v1(int, int, int, int, int);
void CompareProteinWithInterfaces(int, int, int, int, int);
int CompareProteinWithInterfaces_v1(int, int, int, int, int);
int CompareProteinWithInterfaces_v2(int, int, int, int, int);
bool CompareProteinSurfaceWithInterfaceSide(std::map<int, int>, std::map<int, int>, float);



//Debugging
void TestGeometricCalculations(void);
void TestComparingHashTables(void);
void PrintHashTable(std::unordered_map<float, std::unordered_map<float, std::unordered_map<float, std::vector <std::string>>>>&);
int ExtractFragmentsInClusters(std::list<int>, int, int);
int CalculateCosineDistance(int, int, int, int, int, int);


//Old versions
int ExtractFragments(int, int);
int ClusterFragments(int, int, int, int, int, int);
int ClusterFragments_v1(int, int, int, int, int, int);
int ClusterFragments_v2_Parallel(int, int, int, int, int, int);
int CreateInterfaceDescriptors(int, int, int, int, int);
int CreateInterfaceDescriptors_v1(int, int, int, int, int);
int CreateInterfaceDescriptors_v2(int, int, int, int, int);
int CreateInterfaceDescriptors_v2_LessRotations(int, int, int, int, int);
void LoadHashTableFromFile(std::ifstream &, int, std::unordered_map<float, std::unordered_map<float, std::unordered_map<float, std::vector <std::string>>>>&);
int ReadNumberOfClusters(int, int, int, int);

#endif
