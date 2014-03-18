
int mm = 3;
int go = 11;
int ge = 4;
float MaxDiff = 0.99;
int MaxScore;

struct SamData {
    char * name;
    char * chr;
    char * cigar;
    int pos;
    int miss;
    int opens;
    int extension;

    char * muttype;
    int * mutlen;

    float AS;

    struct SamData * next;
};

struct SamData * Shead = NULL;
struct SamData * Scurr = NULL;
struct SamData * Sprev = NULL;

struct SamData * CreateStruct(void)
{
    struct SamData * ptr = (struct SamData *)calloc(1, sizeof(struct SamData));

    ptr->next = NULL;
    ptr->name = NULL;
    ptr->chr = NULL;
    ptr->cigar = NULL;
    ptr->muttype = NULL;
    ptr->mutlen = NULL;

    return ptr;
}

char * ConcatenateStrings(char * s1, char * s2)
{
    char * destination = (char *)calloc(strlen(s1) + strlen(s2) + 1, sizeof(char));

    strcpy(destination, s1);
    strcat(destination, s2);

    return destination;
}

int CompareStrings(char * s1, char * s2)
{
    if(strlen(s1) == strlen(s2))
    {
        if(strncmp(s1, s2, strlen(s1)) == 0)
            return 1;
    }

    return 0;
}

/*char * CopyString(char * source, int length)
{
    char * destination = (char *)calloc(length + 1, sizeof(char));
    strncpy(destination, source, length);

    return destination;
}*/

void CalculateAlignmentScore(struct SamData * psr)
{
    float tempScore = 0;

    if(psr->miss > 0)
        tempScore += (float)psr->miss * (float)mm;

    if(psr->opens > 0)
        tempScore += (float)psr->opens * (float)go;

    if(psr->extension > 0)
        tempScore += (float)psr->extension * (float)ge;

    psr->AS = ((float)MaxScore - tempScore) / ((float)MaxScore);
}

void ThrowStruct(struct SamData * tempo)
{
    if(tempo->cigar)
        free(tempo->cigar);

    if(tempo->name)
        free(tempo->name);

    if(tempo->chr)
        free(tempo->chr);

    if (tempo) free(tempo);
}

void SortStruct(struct SamData * tempStruct)
{
    //printf("%.3f - %.3f\n", Shead->AS, tempStruct->AS);

    if(tempStruct->AS > Shead->AS)
    {
        Sprev = Shead;
        Shead = tempStruct;
        tempStruct->next = Sprev;

        return;
    }

    if(tempStruct->AS > Shead->AS * MaxDiff)
    {
        Scurr->next = tempStruct;
        Scurr = Scurr->next;
        return;
    }

    ThrowStruct(tempStruct);
}

void parseCigar(struct SamData * tempPoint)
{
    int count = 0;
    int toks = 0;
    int length = strlen(tempPoint->cigar);
    int i;
    int gops = 0;
    int gext = 0;

    for(i = 0; i < length; i++)
    {
       if(!isdigit(tempPoint->cigar[i]))
        toks++;
    }

    tempPoint->mutlen = (int *)calloc(toks + 1, sizeof(int));
    tempPoint->muttype = (char *)calloc(toks + 1, sizeof(char));

    toks = 0;

    tempPoint->mutlen[toks] = atoi(tempPoint->cigar);

    while(count < length)
    {
        if(!isdigit(tempPoint->cigar[count]))
        {
            tempPoint->muttype[toks] = tempPoint->cigar[count];
            toks++;
            count++;

            if(count < length)
            {
                tempPoint->cigar += count;
                tempPoint->mutlen[toks] = atoi(tempPoint->cigar);
                tempPoint->cigar -= count;
            }
        }

        count++;
    }

    //printf("%s - %d\n", tempPoint->cigar, toks);

    for(i = 0; i < toks; i++)
    {
        //printf("%c\n", tempPoint->muttype[i]);
        if(tempPoint->muttype[i] == 'D' || tempPoint->muttype[i] == 'I')
        {
            gops++;
            gext += tempPoint->mutlen[i] - 1;
        }
    }

    toks = tempPoint->miss;
    tempPoint->miss = toks - gops - gext;
    tempPoint->opens = gops;
    tempPoint->extension = gext;

    CalculateAlignmentScore(tempPoint);

    /*if(tempPoint->extension > 0)
    {
        printf("AS: %.3f (mis: %d   gapo: %d   gape: %d)\n", tempPoint->AS, tempPoint->miss, tempPoint->opens, tempPoint->extension);
        printf("\n");
    }*/

    if (tempPoint->mutlen) free(tempPoint->mutlen);
    if (tempPoint->muttype) free(tempPoint->muttype);

    SortStruct(tempPoint);
}

void getRestAlignments(char * source)
{
    int toks = 0;
    char * psr = strtok(source, ",");
    char * tpoint;
    struct SamData * pst;

    while(psr != NULL)
    {
        if(toks % 3 == 0)
        {
            if(toks == 0)
            {
                pst = CreateStruct();
                pst->name = CopyString(Shead->name, strlen(Shead->name));
                //printf("\nchr: %s\n", psr);
            }


            else
            {
                //printf("dist: %d\n", atoi(psr));
                pst->miss = atoi(psr);
                parseCigar(pst);
                tpoint = strstr(psr, ";");

                tpoint++;

                if(strlen(tpoint) > 0)
                {
                    tpoint += 1;
                    pst = CreateStruct();
                    pst->name = CopyString(Shead->name, strlen(Shead->name));
                   // printf("\nchr: %s\n", tpoint);
                }

            }
        }


        if(toks % 3 == 1)
        {
            psr += 1;
            pst->pos = atoi(psr);
           // printf("pos %d\n", atoi(psr));
            psr -= 1;
        }

        if(toks % 3 == 2)
        {
            pst->cigar = CopyString(psr, strlen(psr));
            //printf("cigar: %s\n", psr);
        }


        toks++;
        psr = strtok(NULL, ",");
    }

    toks++;
}

void TokenizeRest(char * source)
{
    if(strncmp(source, "XM:i:", 5) == 0)
    {
        source += 5;
        Shead->miss = atoi(source);
        //printf("mastermiss: %d\n", Shead->miss);
        source -= 5;
    }

    if(strncmp(source, "XO:i:", 5) == 0)
    {
        source += 5;
        Shead->opens = atoi(source);
        source -= 5;
    }

    if(strncmp(source, "XG:i:", 5) == 0)
    {
        source += 5;
        Shead->extension = atoi(source);
        source -= 5;
    }
}

void FreeStruct(void)
{
    while(Shead != NULL)
    {
        Scurr = Shead;
        Shead = Shead->next;

        if(Scurr->name)
            free(Scurr->name);

        if(Scurr->chr)
            free(Scurr->chr);

        if(Scurr->cigar)
            free(Scurr->cigar);

        if (Scurr) free(Scurr);
    }

    Shead = NULL;
    Scurr = NULL;
    Sprev = NULL;
}

/*void ReadSAM(char * input, FILE * out)
{
    char * ptr;
    char * temp = NULL;
    int pos;
    int token = 0;

    //printf("%s\n", input);

    if(Shead != NULL)
        FreeStruct();

    Shead = CreateStruct();
    Scurr = Shead;

    ptr = strtok(input, "\t");

    while(ptr != NULL)
    {
        if(token == 0)
            Shead->name = CopyString(ptr, strlen(ptr));

        if(token == 2)
        {
            if(CompareStrings(ptr, "*") == 0)
               Shead->chr = CopyString(ptr, strlen(ptr));
        }

        if(token == 3)
            Shead->pos = atoi(ptr);

        if(token == 3)
            Shead->cigar = CopyString(ptr, strlen(ptr));

        if(token == 9)
            MaxScore = strlen(ptr) * 3;


        if(token > 10)
        {
            if(strncmp(ptr, "XA:Z:", 4) == 0)
            {
                Shead->extension = Shead->extension - Shead->opens;
                CalculateAlignmentScore(Shead);

                ptr += 4;
                temp = CopyString(ptr, strlen(ptr) - 1);
                ptr -= 4;

                //getRestAlignments(temp);
            }

            else
                TokenizeRest(ptr);
        }

        ptr = strtok(NULL, "\t");
        token++;
    }

    if(temp)
    {
        getRestAlignments(temp);

        if(temp)
            free(temp);
    }
}*/

/*void ReadFile(char * input)
{
    FILE * handle = fopen(input, "r+");
    char * outName = ConcatenateStrings(input, ".aln");
    FILE * output = fopen(outName, "w+");
    char line[1000000];

    while(fgets(line, sizeof(line), handle))
    {
        if(line[0] != '@')
            ReadSAM(line, output);
    }

    if(Shead != NULL)
        FreeStruct();

    fclose(handle);
    fclose(output);
    free(outName);
}

void InputFiles(void)
{
    char filename[] = "ion.sam";

    ReadFile(filename);
}*/

/*int main()
{
    InputFiles();
    printf("Hello world!\n");
    return 0;
}*/

