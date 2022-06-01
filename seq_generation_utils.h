#include <iostream>
#include <string.h>
#include <limits.h>
#include <limits>
#include <cmath>
#include "defs.h"

using namespace std;

// debug binary
void bin(edge_id_t n)
{
    edge_id_t i;
    for (i = edge_id_t(1) << (sizeof(edge_id_t) * CHAR_BIT - 1); i > 0; i = i / 2)
    {

      if((n & i) != 0)
      {
        cout << "1";
      }
      else
      {
        cout << "0";
      }
    }
    cout << endl;
}

void generate_unique_weights(const edge_id_t m,  weight_t* weights,
                             const weight_t min_weight,
                             const weight_t max_weight,
                             const edge_id_t offset = 0){
  char significant_bits = 0;
  edge_id_t current = m;
  edge_id_t double_part_mask = 1;
  while(current){
    double_part_mask <<= 1;
    current = current / 2;
    significant_bits++;
  }
  edge_id_t shift = double_part_mask;

  for(char i=significant_bits; i < sizeof(edge_id_t) * CHAR_BIT; i++){
    shift <<= 1;
    double_part_mask = double_part_mask | shift;
  }
  double random;
  edge_id_t bits;
  for(edge_id_t i=0; i < m; i++){
    random = min_weight + (max_weight-min_weight)*drand48();
    memcpy(&bits, &random, sizeof(double));
    bits = bits & double_part_mask | (i + offset);
    memcpy(&weights[i], &bits, sizeof(double));

  }
}

/* write graph to file */
void writeBinaryGraph(graph_t *G, char *filename)
{
    FILE *F = fopen(filename, "wb");
    if (!F) error(EXIT_FAILURE, 0, "Error in opening file %s", filename);

    assert(fwrite(&G->n, sizeof(vertex_id_t), 1, F) == 1);

    assert(fwrite(&G->m, sizeof(edge_id_t), 1, F) == 1);
    assert(fwrite(&G->directed, sizeof(bool), 1, F) == 1);
    uint8_t align = 0;
    assert(fwrite(&align, sizeof(uint8_t), 1, F) == 1);

    assert(fwrite(G->rowsIndices, sizeof(edge_id_t), G->n+1, F) == G->n+1);
    assert(fwrite(G->endV, sizeof(vertex_id_t), G->rowsIndices[G->n], F) == G->rowsIndices[G->n]);
    assert(fwrite(G->weights, sizeof(weight_t), G->m, F) == G->m);
    fclose(F);
}


/* print graph */
void printGraph(graph_t *G)
{
	vertex_id_t i;
    edge_id_t j;
	for (i = 0; i < G->n; ++i) {
		printf("%d:", i);
		for (j=G->rowsIndices[i]; j < G->rowsIndices[i+1]; ++j)
			printf("%d (%f), ", G->endV[j], G->weights[j]);
		printf("\n");
	}
}

void freeGraph(graph_t *G) {
    free(G->rowsIndices);
    free(G->endV);
    free(G->weights);
}

void readGraph(graph_t *G, char *filename)
{
    uint8_t align;
    FILE *F = fopen(filename, "rb");
    if (!F) error(EXIT_FAILURE, 0, "Error in opening file %s", filename);

    assert(fread(&G->n, sizeof(vertex_id_t), 1, F) == 1);
    G->scale = log(G->n) / log (2);

    assert(fread(&G->m, sizeof(edge_id_t), 1, F) == 1);
    assert(fread(&G->directed, sizeof(bool), 1, F) == 1);
    assert(fread(&align, sizeof(uint8_t), 1, F) == 1);

    G->rowsIndices = (edge_id_t *)malloc((G->n+1) * sizeof(edge_id_t));
    assert(G->rowsIndices);
    assert(fread(G->rowsIndices, sizeof(edge_id_t), G->n+1, F) == (G->n+1));
    G->endV = (vertex_id_t *)malloc(G->rowsIndices[G->n] * sizeof(vertex_id_t));
    assert(G->endV);
    assert(fread(G->endV, sizeof(vertex_id_t), G->rowsIndices[G->n], F) == G->rowsIndices[G->n]);
    G->weights = (weight_t *)malloc(G->m * sizeof(weight_t));
    assert(G->weights);

    assert(fread(G->weights, sizeof(weight_t), G->m, F) == G->m);
    fclose(F);
}
