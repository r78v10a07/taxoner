int maxNode;

struct nodes {
    int id;
    int nextId;
    int ok;
};


struct nodes * n = NULL;

void FreeData(void)
{
    if(n != NULL)
        free(n);
}


void AllocateNode(void)
{
    int i;

    if(maxNode < 1)
    {
        printf("No taxons exist\ncheck ncbi nodes.dmp file\n");
        exit(EXIT_FAILURE);
    }

    n = (struct nodes *)calloc(maxNode, sizeof(struct nodes));

    for(i = 0; i < maxNode; i++)
    {
        n[i].id = 0;
        n[i].nextId = 0;
        n[i].ok = 0;
    }
}

void FillNodes(char * infile)
{
    FILE * handle = fopen(infile, "r+");
    char line[BUFSIZ];
    int tempnode;
    int token = 0;
    char * ptr;

    while(fgets(line, sizeof(line), handle))
    {
        ptr = strtok(line, "|");
        token = 0;

        while(ptr != NULL)
        {
            if(token == 0)
                tempnode = atoi(ptr);

            if(token == 1)
            {
                n[tempnode].id = tempnode;
                n[tempnode].nextId = atoi(ptr);
            }

            ptr = strtok(NULL, "|");
            token++;
        }
    }

    fclose(handle);
}

void GetMaxNodes(char * infile)
{
    FILE * handle = fopen(infile, "r+");
    char line[BUFSIZ];
    int tempnode = 0;

    while(fgets(line, sizeof(line), handle))
    {
        tempnode = atoi(line);

        if(tempnode > maxNode)
            maxNode = tempnode + 1;
    }

    printf("Highest taxon (node): %d\n", maxNode - 1);fflush(NULL);
    maxTax = maxNode;
    fclose(handle);

    AllocateNode();
    FillNodes(infile);
}
