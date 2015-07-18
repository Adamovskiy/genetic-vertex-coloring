#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

const int GENERATIONS_COUNT = 1000000;
const float MUTATION_RISK =0.0001;
const int GENERATION_SIZE = 10000; // must be even

const float ZERO_THRESHOLD = 0.00000001; // if number is less that this, consider as zero

int verticesCount = -1;
int edgesCount = -1;

typedef struct {
	int adjacentsCount;
	int * adjacents; // array of adjacent vertex indexes
} vertex;

vertex * vertices;

void createVertices() {
	vertices = (vertex *)malloc(sizeof(vertex)*verticesCount);
	for (int i = 0; i < verticesCount; i++) {
		vertices[i].adjacentsCount = 0;
		vertices[i].adjacents = (int *)malloc(sizeof(int)*verticesCount); // TODO optimize size
	}
}

void createEdgeOneDirection(int a, int b) {
	for (int i = 0; i<vertices[a].adjacentsCount; i++) {
		if (vertices[a].adjacents[i] == b)
			return;
	}
	vertices[a].adjacents[vertices[a].adjacentsCount] = b;
	vertices[a].adjacentsCount++;
}

// bi-directional
void createEdge(int a, int b) {
	if (a==b) {
        printf("Bad graph structure: vertex %i is connected with itself.\n", a);
		exit(1);
	}
	createEdgeOneDirection(a,b);
	createEdgeOneDirection(b,a);
}

/*
 * Read graph from given file. Format discription and data files can be found at
 * http://www.info.univ-angers.fr/pub/porumbel/graphs/.
 */
void readFile(const char* filename) {
	FILE * file = fopen( filename, "r" );
	if (file == NULL) {
        printf("Input file not found: %s\n", filename);
		exit(1);
	}
	char lineDescriptor = 'x';
	char edge[5];
	while (!feof(file)) {
		fscanf(file, "%c", &lineDescriptor);
		int e1,e2;
		switch (lineDescriptor) {
		case 'c':
			do {
				if (fgets(edge, 5, file) == NULL)
                    break;
			}while (NULL != edge && edge[strlen(edge)-1] != '\n' );
		break;
		case 'p':
			fscanf(file, " %s %i %i\n", edge, &verticesCount, &edgesCount);
			createVertices();
		break;
		case 'e':
			fscanf(file, " %i %i\n", &e1, &e2);
			createEdge(e1-1, e2-1); // numeration starts from 0
		break;
		}
	}
	fclose( file );
}


typedef struct {
	int * colors; // colors of vertices. 0 means uncolored
	int chromatic; // used colors count
	float fitness;
	float likelihood;
	int errors; // count of edges (one direction) between vertices of same color.
} individual;

individual * generation;
individual * newGeneration;

individual* createIndividuals() { // returns empty generation
	individual * result = (individual *)malloc(sizeof(individual)*GENERATION_SIZE);
	for (int i = 0; i < GENERATION_SIZE; i++) {
		result[i].colors = (int *)malloc(sizeof(int)*verticesCount);
		for (int j = 0; j < verticesCount; j++) {
			result[i].colors[j] = 0; // uncolored
		}
		result[i].fitness = -1; // not computed
		result[i].chromatic = 0; // no colors yet
		result[i].errors = -1;
		result[i].likelihood = 0;
	}
	return result;
}

float calcFitness(individual* ind) {
	float COEF1 = verticesCount*edgesCount*1000; // value of miscolored edges fraction in answer rating
	float COEF2 = .1; // value of chromatic number in answer rating

	int miscoloredEdges = 0;
	for (int i = 0; i < verticesCount; i++) {
		for (int j = 0; j < vertices[i].adjacentsCount; j++) {
			if (ind->colors[i] == ind->colors[vertices[i].adjacents[j]])
				miscoloredEdges++;
		}
	}
	ind->errors = miscoloredEdges;
	miscoloredEdges /= 2; // halve it, because edges are bidirectional
	ind->fitness = ( ((float)miscoloredEdges)/edgesCount ) * COEF1 + ind->chromatic * COEF2;
	return ind->fitness;
}

//ubrat durki, perekrashivaet s poslednego, s4itaet kolvo cvetov
void shrinkAndCalcChromatic(individual* ind) {
	int* existingColors = (int *)malloc(sizeof(int)*(verticesCount+1));
	while (1) {
		for (int i = 1; i <= verticesCount; i++)
			existingColors[i] = 0;
		int maxColor = 0;
		for (int i = 0; i < verticesCount; i++) {
			if (ind->colors[i] > maxColor)
				maxColor = ind->colors[i];
			existingColors[ind->colors[i]] = 1;
		}

		int colorToRemove = 0;
		for(int i = 1; i < maxColor; i++) {
			if (existingColors[i] == 0) {
				colorToRemove = i;
				break;
			}
		}
		if (colorToRemove == 0) {
			ind->chromatic = maxColor;
			free(existingColors);
			return;
		}
		for (int i = 0; i < verticesCount; i++) {
			if (ind->colors[i] == maxColor)
				ind->colors[i] = colorToRemove;
		}
	}
}

void setRandomColoring(int num) {
	individual* ind = &generation[num];
	int colorsRange = rand() % verticesCount + 1;
	for (int i = 0; i < verticesCount; i++) {
		ind->colors[i] = rand() % colorsRange + 1;
	}
	shrinkAndCalcChromatic(ind);
}

int isAlmostZero(float f) {
	return (
		(f > 0 && ZERO_THRESHOLD > f) ||
		(f < 0 && -ZERO_THRESHOLD < f)
	) ? 1 : 0;
}

// will return index of solution, if any; -1 othervise
int createFirstGeneration () {
	generation = createIndividuals();
	for (int i = 0; i < GENERATION_SIZE; i++) {
		setRandomColoring(i);
		if (isAlmostZero(calcFitness(&generation[i]))) {
			return i;
		}
	}
	return -1;
}

// will return index of solution, if any; -1 othervise
int growNewGeneration() {
	for (int i = GENERATION_SIZE - 1; i >= 0; i--) {
		free(generation[i].colors);
	}
	free(generation);
	generation = newGeneration;
	for (int i = 0; i < GENERATION_SIZE; i++) {
		if (isAlmostZero(calcFitness(&generation[i]))) {
			return i;
		}
	}
	return -1;
}

float multiinv = 0; // sum of inverse fitness values

// select random individual for crossover (in accordance with likelhoods)
int selectMate() {
	float random = (rand()%100001) / 100000.0; // 0.000% - 100.000%
	float processedLikelihood = 0;
	for (int i = 0; i < GENERATION_SIZE; i++) {
		if (random > processedLikelihood && random < generation[i].likelihood + processedLikelihood)
			return i;
		processedLikelihood += generation[i].likelihood;
	}
	return GENERATION_SIZE - 1;
}

/*
 * Mark subgraph as visited (1s in visited array)
 * startedFrom - if its value is 1 - was passed to fillGenes as startVertex
 * (need to avoid going to already processed vertices)
 */
int fillGenes(int startVertex, int currentCount, int totalCount, int * visited, int* startedFrom) {
	for (int i = 0; i < vertices[startVertex].adjacentsCount; i++) {
		int newGene = vertices[startVertex].adjacents[i];
		if (visited[newGene] == 1)
			continue;
		visited[newGene] = 1;
		currentCount++;
		if (currentCount == totalCount) {
			return currentCount;
		}
	}
	for (int i = 0; i < vertices[startVertex].adjacentsCount; i++) {
		int newStartVertex = vertices[startVertex].adjacents[i];
		if (startedFrom[newStartVertex] == 1)
			continue;
		startedFrom[newStartVertex] = 1;
		currentCount = fillGenes(newStartVertex, currentCount, totalCount, visited, startedFrom); // recursive call
	}

	return currentCount;
}

void mutate(individual* ind) {
	for (int i = 0; i < verticesCount; i++) {
		if ((rand() % 100001)/100000.0 < MUTATION_RISK) {
			ind->colors[i] = rand() % ind->chromatic + 1;
		}
	}
}

int newGenerationSize;

// Father genes are selected as random vertex, and its adjacents, and their adjacents, and so on.
void breed(int father, int mother, int* visited, int* startedFrom) {
	int son = newGenerationSize;
	int daughter = newGenerationSize + 1;

	// fathers genes
	int balance = rand() % verticesCount; // domination in accordance with likelhood can be emulated by, for example, = generation[father].likelhood * verticesCount
	for (int i = 0; i < verticesCount; i++) {
		visited[i] = 0;
		startedFrom[i] = 0;
	}
	int currentCount = 0;
	do {
		int initialVertex;
		do {
			initialVertex = rand() % verticesCount;
		} while (visited[initialVertex] != 0); // it means that we hit already added detached subgraph
		visited[initialVertex] = 1;
		currentCount++;
		startedFrom[initialVertex] = 1;
		currentCount = fillGenes(initialVertex, currentCount, balance, visited, startedFrom);
	}while(currentCount < balance);

    // son gets visited genes from father, not visited from mother, daughter - vice versa
	for (int i = 0; i < verticesCount; i++) {
		if (visited[i]) {
			newGeneration[son].colors[i] = generation[father].colors[i];
			newGeneration[daughter].colors[i] = generation[mother].colors[i];
		}
		else {
			newGeneration[son].colors[i] = generation[mother].colors[i];
			newGeneration[daughter].colors[i] = generation[father].colors[i];
		}
	}
	shrinkAndCalcChromatic(&newGeneration[son]);
	shrinkAndCalcChromatic(&newGeneration[daughter]);
	mutate(&newGeneration[son]);
	mutate(&newGeneration[daughter]);
	shrinkAndCalcChromatic(&newGeneration[son]);
	shrinkAndCalcChromatic(&newGeneration[daughter]);

	newGenerationSize += 2;
}

void crossover() {
	for (int i = 0; i < GENERATION_SIZE; i++) {
		multiinv += 1.0/generation[i].fitness;
	}
	for (int i = 0; i < GENERATION_SIZE; i++) {
		generation[i].likelihood = (1.0 / generation[i].fitness) / multiinv;
	}

	newGeneration = createIndividuals();
	newGenerationSize = 0;
	int* visited = (int *)malloc(sizeof(int)*verticesCount); // avoid memory re-allocation in each loop
	int* startedFrom = (int *)malloc(sizeof(int)*verticesCount);
	for (int i = 0; i < GENERATION_SIZE / 2; i++) {
		int a = selectMate();
		int b;
		do { // avoid self-breed
			b = selectMate();
		} while (b == a);
		breed(a, b, visited, startedFrom);
	}
	free(visited);
	free(startedFrom);
};

void printIndividual(individual* ind) {
	printf("chromatic number = %i and %i errors is:\n", ind->chromatic, ind->errors);	
	// Comment out this loop to get more compact output	
	for (int i = 0; i < verticesCount; i++) {
		printf("%i\t=\t%i\n", i, ind->colors[i]);
	}
}

void findAndPrintBest(int generationNumber) {
	float minFitness = generation[0].fitness;
	int minFitnessIndex = 0;
	for (int i = 0; i < GENERATION_SIZE; i++) {
		if (generation[i].fitness < minFitness)
			minFitnessIndex = i;
	}
	printf("In generation %i the best coloring: ", generationNumber);
	printIndividual(&generation[minFitnessIndex]);
}

void freeVertices() {
for (int i = 0; i < verticesCount; i++) {
		free(vertices[i].adjacents);
	}
	free(vertices);
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("Usage: GeneticVertexColoring [input-file]\n");
        exit(0);
    }

	srand((unsigned) time(NULL));

	readFile(argv[1]);

	int result = createFirstGeneration();
	if (result != -1) {
		printf("Got solution: ");
		printIndividual(&generation[result]);
	}
	else {
		findAndPrintBest(0);
		for (int i = 0; i < GENERATIONS_COUNT; i++) {
			crossover();
			result = growNewGeneration();
			if (result != -1) {
				printf("Got solution: ");
				printIndividual(&generation[result]);
				break;
			}
			else {
				findAndPrintBest(i+1);
			}
		}
	}
	freeVertices();
	char c;
	scanf("%c", &c);
	return 0;
}
