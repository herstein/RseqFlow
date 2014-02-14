#include "Merge_genome_transcriptome.h"


//read in junctionIndex files
void ProcessJunction(char* fileName)
{
	FILE *fp;
	fp = fopen(fileName, "r");
	printf("Begin to process %s!\n", fileName);
	if (fp == NULL)
	{
		printf("The junction file doesn't exist!\n");
		exit(1);
	}
	char oneLine[50000];
	char *start, *end, *comma;
	int a, b, c, d;
	char *isophone;
	while (fgets(oneLine, 50000, fp) != NULL)
	{
		comma = strchr(oneLine , ',');
		while (comma != NULL)
		{
			*comma = ' ';
			comma = strchr(++comma, ',');
		}
		if (strlen(oneLine) < 10)
		{
			break;
		}
		isophone = new char[20];
		sscanf(oneLine, "%s", isophone);
		start = strchr(oneLine, '(');

		while (start != NULL)
		{
			start++;
			end = strchr(start, ')');
			if (end == NULL)
			{
				printf("Error: right bracket missing!\n");
				exit(-1);
			}
			*end = ' ';
			sscanf(start, "%d %d %d %d", &a, &b, &c, &d);

			vector<int> temp;
			temp.push_back(a);
			temp.push_back(b);
			temp.push_back(c);
			temp.push_back(d);

			const char* isophone2 = (const char*)isophone;

			if (junctions.find(isophone2) == junctions.end())
			{
				vector< vector<int> > temp2;
				temp2.push_back(temp);
				junctions[isophone2] = temp2;
			}
			else
				junctions[isophone2].push_back(temp);
			start = strchr(end, '(');	

		}

		
	}
	fclose(fp);


}

void GencodePreprocess()
{
	FILE *fp = fopen(gencodeName, "r");
	printf("Begin gencode preprocessing!\n");
	if (fp == NULL)
	{
		printf("The gencode file doesn't exist!\n");
		exit(1);
	}
	char oneLine[500];
	while (fgets(oneLine, 500, fp) != NULL)
	{
		gencodeTotalLineNum++;
	}
	flagTable = new char[(int)gencodeTotalLineNum + 1];
	/*for (int i = 1; i < gencodeTotalLineNum + 1 ; i++ )
	{
		flagTable[i] = 'N'; 
	}*/

	//flagTable = new char[(int)1000000 + 1];


}

void CreateHashTable()
{

	char oneLine[500];
	char *readID, *chrIDLocation;
	char chrID[10];
	char location[10];
	rounds++;
	char previousReadID[80];
	char currentReadID[80];

	printf("Round %d: CreateHashTable!\n", rounds);
	while (fgets(oneLine, 500, genomeFp) != NULL)
	{

		lineNum++;
		//printf("%d\n", lineNum);
		/*if (lineNum > 2252310)
		{
			printf("Line %d: %s\n", lineNum, oneLine);
		}*/
		
		//sleep(1);
		genomeTotalLineNum++;
		readID = new char[80];
		chrIDLocation = new char[50];
		//printf("start to read..\n");
		sscanf(oneLine, "%s %*s %s %s", readID, chrID, location);
		if (chrID[0]=='*')
                {
			delete readID;
                        delete chrIDLocation;
                        continue;
                }
		//printf("sscanf complete %s\n", chrID);
		sprintf(chrIDLocation, "%s:%s", chrID, location);
		//printf("sprintf complete\n");
		strcpy(currentReadID, readID);
		//printf("copy readID %s\n", readID);
		if (hashTable.find(readID) == hashTable.end())
		{
			hashTable[readID].push_back(chrIDLocation);

			toBeDeleted.push_back(readID);
			toBeDeleted.push_back(chrIDLocation);

		}
		else
		{
			bool found = false;
			for (unsigned int i = 0; i < hashTable[readID].size() ; i++ )
			{
				if (strcmp(hashTable[readID][i], chrIDLocation) == 0)
				{
					found = true;
					break;
				}
			}
			if (found == false)
			{
				hashTable[readID].push_back(chrIDLocation);
				toBeDeleted.push_back(readID);
				toBeDeleted.push_back(chrIDLocation);
			}
			else 
			{
				delete readID;
				delete chrIDLocation;
			}
		}
		if ((double)hashTable.size() >= 3000000 && strcmp(previousReadID, currentReadID) != 0)
		{
			break;
		}
		strcpy(previousReadID, currentReadID);


	}

	if (lineNum == 0)
	{
		return;
	}
	ProcessGencode();

}

void ProcessGencode()
{

	vector<char*> temporary;
	
	printf("Round %d: Process gencode.sam\n", rounds);
	FILE *fp = fopen(gencodeName, "r");
	

	char oneLine[500];
	gencodeCurrentLineNum = 0;
	char previousReadID[80];		



	while (fgets(oneLine, 500, fp) != NULL)
	{
		gencodeCurrentLineNum++;
		//printf("%d\n", gencodeCurrentLineNum);
		if (flagTable[gencodeCurrentLineNum] == 'P')
		{
			continue;
		}
		/*char *pch = strchr(oneLine, '_');
		pch = strchr(pch+1, '_');
		pch = strchr(pch+1, '_');
		pch = strchr(pch+1, '_');
		*pch = ' ';
		pch = strchr(pch + 1, '_');
		*pch = ' ';
		pch = strchr(pch + 1, '=');
		*pch = ' ';*/
		char readID[80];

		int secondColumn = 0;
		char chr_isoID[30];
		char chrID[10];
		int startInGenome = 0;
		char thirdColumn[100];
		int fifthColumn = 0;
		int firstHalf = 0;
		int middle = 0;
		int secondHalf = 0;
		char CIGAR[20];
		char flag;
		int eighthColumn = 0;
		int ninthColumn = 0;
		char series1[80];
		char series2[80];
		char twelfthColumn[20];
		//twelfthColumn[0] = '\0';
		char thirteenthColumn[20];
		thirteenthColumn[0] = '\0';
		char lastColumn[20];
		lastColumn[0] = '\0';
		//char geneID[15];
		char isophone[20];
		int startInGencode = 0;

		sscanf(oneLine, "%s %d %s %d %d %s %c %d %d %s %s %s %s %s %*s", readID, &secondColumn, thirdColumn, &startInGencode, &fifthColumn, CIGAR, &flag, &eighthColumn, &ninthColumn, series1, series2, twelfthColumn, thirteenthColumn, lastColumn);
		if (thirdColumn[0]=='*')
                {
                        continue;
                }
		//sscanf(oneLine, "%s %*s", readID);
		//printf("%s\n", readID);	
		/*sscanf(oneLine, "%s %d %*s %*s %s  %s   %d  %*s %c %d %d %s %s %s %s %s",
			readID, &secondColumn, isophone, chr_isoID, &fifthColumn, &flag, &eighthColumn, &ninthColumn, series1, series2, twelfthColumn, thirteenthColumn, lastColumn);*/
		sscanf(CIGAR, "%d[^M]M", &readLength);
		//printf("%d\n", readLength);
		char *pch = strchr(thirdColumn, '_');
		pch = strchr(pch + 1, '_');
		*pch = ' ';
		pch = strchr(pch + 1, '=');
		*pch = ' ';
		sscanf(thirdColumn, "%*s %s %s", isophone, chr_isoID);		

		pch = strchr(chr_isoID, ':');
		if (pch == NULL)
		{
			printf("Format error: missing ':' between chromID & isophone!\n");
			printf("%s\n",chr_isoID);
			exit(-1);
		}
		*pch = ' ';
		sscanf(chr_isoID, "%s", chrID);


		
		if (strcmp(previousReadID, readID) != 0)
		{
			if (hashTable.find(previousReadID) != hashTable.end())
			{
				for (unsigned int i = 0; i < temporary.size() ; i++)
				{
					fprintf(multipleFp, "%s", temporary[i]);
					delete temporary[i];
					numMapped++;
				}
			
				temporary.clear();

			}
		}
		
		strcpy(previousReadID, readID);

		if (hashTable.find(readID) != hashTable.end())
		{
			bool found2 = false;
			if (junctions.find((const char*)isophone) != junctions.end())
			{
				
				for (unsigned int i = 0; i < junctions[(const char*)isophone].size() ; i++ )
				{
					if ((startInGencode + readLength - 1 > junctions[(const char*)isophone][i][0]) && (startInGencode <= junctions[(const char*)isophone][i][0]))
					{
						found2 = true;
						firstHalf = junctions[(const char*)isophone][i][0] - startInGencode + 1;
						middle = junctions[(const char*)isophone][i][3] - junctions[(const char*)isophone][i][2];
						secondHalf = readLength - firstHalf;
						startInGenome = junctions[(const char*)isophone][i][2] - firstHalf + 1; 
						break;
					}
				}
			}


			if (found2 == true)
			{
			
				char *chrIDLocation = new char[50];
				sprintf(chrIDLocation, "%s:%d", chrID, startInGenome);
				bool found = false;
				for (unsigned int i = 0; i < hashTable[readID].size(); i++ )
				{
					if (strcmp(hashTable[readID][i], chrIDLocation) == 0)
					{
						found = true;
						delete chrIDLocation;
						break;
					}
				}
				if (found == false)
				{
					hashTable[readID].push_back(chrIDLocation);
					char *oneEntry = new char[300];
					sprintf(oneEntry, "%s\t%d\t%s\t%d\t%d\t%dM%dN%dM\t%c\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n", readID, secondColumn, chrID, startInGenome, fifthColumn, firstHalf, middle, secondHalf, flag, eighthColumn, ninthColumn, series1, series2, twelfthColumn, thirteenthColumn, lastColumn);
					toBeDeleted.push_back(chrIDLocation);
					temporary.push_back(oneEntry);
				}

			}
			flagTable[gencodeCurrentLineNum] = 'P';
		}
				
		//delete readID;	
	}
	
	
	ProcessGenome();
}

void ProcessGenome()
{

	printf("Round %d: Process genome.sam \n",rounds);
	char oneLine[500];

	char readID[80];
	char chrID[10];
	char location[10];


	while (lineNum != 0)
	{
		if (fgets(oneLine, 500, genomeFp2) == NULL)
		{
			printf("Error: can't read lines in second genome processing\n");
			exit(-1);
		}
		lineNum--;
		
		sscanf(oneLine, "%s %*s %s %s", readID, chrID, location);
		if (chrID[0]=='*')
		{
			continue;
		}
		if (hashTable.find(readID) == hashTable.end())
		{
			printf("Error: can't find readID in second genome processing\n");
			exit(-1);

		}
		else
		{
			if (hashTable[readID].size() == 1)
			{
				fprintf(uniqueFp, "%s", oneLine);
				uniqueMapped++;
				numMapped++;
			}
			else 
			{
				fprintf(multipleFp, "%s", oneLine);
				numMapped++;
			}
		}
		


	}
	


	hashTable.clear();
	while (toBeDeleted.size() > 0)
	{
		//printf("%s\n", toBeDeleted[i]);
		char *temp = toBeDeleted.back();
		toBeDeleted.pop_back();
		delete temp;
		//printf("After: %s\n", toBeDeleted[i]);
	}
	toBeDeleted.clear();
	

	CreateHashTable();

}

void FinalGencodeProcess()
{
	printf("Start final gencode processing!\n");
	vector<char*> temporary;
	FILE *fp;
	fp = fopen(gencodeName, "r");
	int lineNumTemp = 0;
	if (fp == NULL)
	{
		printf("The %s doesn't exist!\n", gencodeName);
		exit(-1);
	}
	char oneLine[500];
	char previousReadID[80];		

	while (fgets(oneLine, 500, fp) != NULL)
	{
		lineNumTemp++;
		if (flagTable[lineNumTemp] == 'P')
		{
			continue;
		}
		/*char *pch = strchr(oneLine, '_');
		pch = strchr(pch+1, '_');
		pch = strchr(pch+1, '_');
		pch = strchr(pch+1, '_');
		*pch = ' ';
		pch = strchr(pch + 1, '_');
		*pch = ' ';
		pch = strchr(pch + 1, '=');
		*pch = ' ';*/
		char *readID = new char[80];
		int secondColumn = 0;
		char thirdColumn[100];
		char chr_isoID[30];
		char chrID[10];
		int startInGenome = 0;
		int fifthColumn = 0;
		int firstHalf = 0;
		int middle = 0;
		int secondHalf = 0;
		char flag;
		int eighthColumn = 0;
		int ninthColumn = 0;
		char series1[80];
		char series2[80];
		char twelfthColumn[20];
		//twelfthColumn[0] = '\0';
		char thirteenthColumn[20];
		thirteenthColumn[0] = '\0';
		char lastColumn[20];
		lastColumn[0] = '\0';
		//char geneID[15];
		char isophone[20];
		int startInGencode = 0;
		
		sscanf(oneLine, "%s %d %s %d %d %*s %c %d %d %s %s %s %s %s", readID, &secondColumn, thirdColumn, &startInGencode, &fifthColumn, &flag, &eighthColumn, &ninthColumn, series1, series2, twelfthColumn, thirteenthColumn, lastColumn);
		/*sscanf(oneLine, "%s %d %*s %*s %s  %s   %d %d %*s %c %d %d %s %s %s",
			readID, &secondColumn, isophone, chr_isoID, &startInGencode, &fifthColumn, &flag, &eighthColumn, &ninthColumn, series1, series2, lastColumn);*/
		if (thirdColumn[0]=='*')
                {
                        continue;
                }
		char *pch = strchr(thirdColumn, '_');
		pch = strchr(pch + 1, '_');
		*pch = ' ';
		pch = strchr(pch + 1, '=');
		*pch = ' ';
		sscanf(thirdColumn, "%*s %s %s", isophone, chr_isoID);		

		pch = strchr(chr_isoID, ':');
		if (pch == NULL)
		{
			printf("Format error: missing ':' between chromID & isophone!\n");
			printf("%s\n",chr_isoID);
			exit(-1);
		}
		*pch = ' ';
		sscanf(chr_isoID, "%s", chrID);		


		
		if (strcmp(previousReadID, readID) != 0)
		{
			char *toBeDeleted1;
			hashTableIterator = hashTable.find(previousReadID);
			if (hashTableIterator != hashTable.end())
			{
				for (unsigned int i = 0; i < hashTable[previousReadID].size() ; i++ )
				{
					toBeDeleted1 = hashTable[previousReadID][i];
					delete toBeDeleted1;
				}
				toBeDeleted1 = hashTableIterator->first;
				delete toBeDeleted1;
				hashTable.clear();
				if (temporary.size() > 1)
				{
					for (unsigned int i = 0; i < temporary.size() ; i++)
					{
						fprintf(multipleFp, "%s", temporary[i]);
						delete temporary[i];
						numMapped++;
					}
				}
				else if (temporary.size() == 1)
				{
					fprintf(uniqueFp, "%s", temporary[0]);
					delete temporary[0];
					numMapped++;
					uniqueMapped++;
				}
				temporary.clear();
				hashTable.clear();

			}
		}
		
		strcpy(previousReadID, readID);

		if (hashTable.find(readID) != hashTable.end())
		{

			//fprintf(fp2, "%s\t%d\t%s\t%d\t%d\t%dM%dN%dM\t%c\t%d\t%d\t%s\t%s\t%s\n", readID, secondColumn, chrID, startInGenome, fifthColumn, firstHalf, middle, secondHalf, flag, eighthColumn, ninthColumn, series1, series2, lastColumn);
			char *chrIDLocation = new char[50];
			sprintf(chrIDLocation, "%s:%d", chrID, startInGenome);
			bool found = false;
			for (unsigned int i = 0; i < hashTable[readID].size(); i++ )
			{
				if (strcmp(hashTable[readID][i], chrIDLocation) == 0)
				{
					found = true;
					delete chrIDLocation;
					break;
				}
			}
			if (found == false)
			{
				hashTable[readID].push_back(chrIDLocation);
				char *oneEntry = new char[300];
				sprintf(oneEntry, "%s\t%d\t%s\t%d\t%d\t%dM%dN%dM\t%c\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n", readID, secondColumn, chrID, startInGenome, fifthColumn, firstHalf, middle, secondHalf, flag, eighthColumn, ninthColumn, series1, series2, twelfthColumn, thirteenthColumn, lastColumn);
				hashTable[readID].push_back(chrIDLocation);
				temporary.push_back(oneEntry);
			}

			flagTable[gencodeCurrentLineNum] = 'P';
		}
		else 
		{
			delete readID;				
			
		}

	}

}



int main(int argc, char *argv[])
{
	argc--, argv++;

	char **argvTemp = argv;
	int argcTemp = argc;

	char *junctionName = new char[50];
	gencodeName = new char[50];
	genomeName = new char[50];
	char *uniqueName = new char[50];
	char *multipleName = new char[50];


	//extract 
	int counter = 1;
	for(;argcTemp>0;argcTemp--, argvTemp++, counter++)
	{
		if (counter == 1)
		{
			junctionName = *argvTemp;
		}
		else if (counter == 2)
		{
			gencodeName = *argvTemp;
		}
		else if (counter == 3)
		{
			genomeName = *argvTemp;
		}
		else if (counter == 4)
		{
			uniqueName = *argvTemp;
		}
		else if (counter == 5)
		{
			multipleName = *argvTemp;
		}

		
	}

	if (counter != 6)
	{
		printf("Invalid number of arguments! \n");
		exit(-1);
	}
	
	ProcessJunction(junctionName);
	genomeFp = fopen(genomeName, "r");
	if (genomeFp == NULL)
	{
		printf("The genome file doesn't exist!\n");
		exit(1);
	}
	genomeFp2 = fopen(genomeName, "r");
	if (genomeFp2 == NULL)
	{
		printf("The genome file doesn't exist!\n");
		exit(1);
	}
	gencodeFp = fopen(gencodeName, "r");
	if (gencodeFp == NULL)
	{
		printf("The gencode file doesn't exist!\n");
		exit(1);
	}
	multipleFp = fopen(multipleName,"w");
	uniqueFp = fopen(uniqueName, "w");
	GencodePreprocess();
	CreateHashTable();
	FinalGencodeProcess();
	fclose(genomeFp);
	fclose(gencodeFp);
	fclose(multipleFp);
	fclose(uniqueFp);
	printf("Total number of lines read from %s: %d\n", gencodeName, gencodeTotalLineNum );
	printf("Total number of lines read from %s: %d\n", genomeName, genomeTotalLineNum);
	

	return 0;
}
