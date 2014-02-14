#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctype.h>
#include <list>
#include <vector>
#include <sys/time.h>
#include <strings.h>
#include <time.h>
//#include <tr1/unordered_map>
//#include <tr1/unordered_map>

#include <ext/hash_map>
using namespace __gnu_cxx;
#include<string>
using namespace std;
//using namespace std::tr1;
//#include <hash_map.h>

typedef struct sum
{
	int laneNumber;
	int clusterRaw;
	int clusterFilter;
}Summary;

struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) == 0;
  }
};



vector<Summary *> summaryVector;


hash_map<char*, vector<char*>, hash<char*>, eqstr> readsMapped, hashTable;
//hash_map<const char*, vector<const char*>, hash<const char*>, eqstr> readsUnique;
hash_map<const char*, vector< vector<int> >, hash<const char*>, eqstr> junctions;
hash_map<char*, vector<char* >, hash<char*> , eqstr>::iterator mappedIterator, hashTableIterator, tempIter;
hash_map<const char*, const char*, hash<const char*>, eqstr> chrmReadsMapped;
vector<char *> toBeDeleted;

/*hash_map<const char*, const char*, hash<const char*>, eqstr> readsMapped;
hash_map<const char*, vector<const char*>, hash<const char*>, eqstr> readsUnique;
hash_map<const char*, vector< vector<int> >, hash<const char*>, eqstr> junctions;
hash_map<const char*, const char*, hash<const char*>, eqstr>::iterator mappedIterator;
hash_map<const char*, vector<const char*>, hash<const char*>, eqstr> intermediate;*/


int lineNum = 0;
int gencodeTotalLineNum = 0;
int genomeTotalLineNum = 0;
int gencodeCurrentLineNum = 0;
int rounds = 0;
int totalNum = 0;
int numMapped = 0;
int numMappedGencode = 0;
int numMappedGenome = 0;
int uniqueMapped = 0;
int multipleMapped = 0;
int uniqueMappedGencode = 0;
int uniqueMappedGenome = 0;
int readLength = 0;
double actb = 0.0;
char machineName[100];
char folderName[100];
FILE *gencodeFp, *genomeFp, *multipleFp, *uniqueFp, *genomeFp2;
char *flagTable, *gencodeName, *genomeName;

void Convert(char *thousands, int a);
void ParseHTM(char *fileName, int laneNumber);
void ParseXML(char *fileName, int laneNumber);
void ProcessGencode();
void ProcessGenome();
void ACTBExpression(char *fileName);
void ProcessChrM(char *chrmGencodeName, char *chrmGenomeName);
void CreateHashTable();
