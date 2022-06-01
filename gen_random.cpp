#include "defs.h"
#include <stdio.h>
#include <cstdlib>
#include <assert.h>
#include <string.h>
#include <error.h>
#include "seq_generation_utils.h"


using namespace std;

char outFilename[FNAME_LEN];

/* helper */
void usage(int argc, char **argv)
{
    printf("Random graph generator\n");
    printf("Usage:\n");
    printf("%s -s <scale> [other options]\n", argv[0]);
    printf("Options:\n");
    printf("   -s <scale>, number of vertices is 2^<scale>\n");
    printf("   -k <half the average vertex degree>, default value is 16\n");
    printf("   -out <output filename>, file for the graph storage\n");
    exit(1);
}

/* initialization */
void init(int argc, char **argv, graph_t *G)
{
	bool no_out_filename = true;
    G->scale = -1;
    G->directed = false; // WARNING Graph is always undirected here. TODO: add random generation of directed graphs.
    G->permute_vertices = true;
    G->min_weight = 0;
    G->max_weight = 1;
    /* default value */
    G->avg_vertex_degree = DEFAULT_ARITY;
    if (argc == 1) {
        usage(argc, argv);
    }

    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-s")) {
            G->scale = (int)atoi(argv[++i]);
        }

        if (!strcmp(argv[i], "-k")) {
            G->avg_vertex_degree = (int)atoi(argv[++i]);
        }

		if (!strcmp(argv[i], "-out")) {
            no_out_filename = false;
			sprintf(outFilename, argv[++i]);
        }
    }

	if (no_out_filename) {
    	sprintf(outFilename, "random-%d", G->scale);
	}

    if (G->scale == -1) {
        usage(argc, argv);
    }

    G->n = (vertex_id_t)1 << G->scale;
    G->m = G->n * G->avg_vertex_degree;
}

/* random graph generator */
void gen_random_graph(graph_t *G)
{
    /* init */
    vertex_id_t n;
    edge_id_t m;
    edge_id_t offset;
    bool permute_vertices;
    vertex_id_t *permV, tmpVal;
    vertex_id_t u, v;
    vertex_id_t *src;
    vertex_id_t *dest;
    unsigned *degree;
    permute_vertices = G->permute_vertices;
    double *dbl_weight;
    double min_weight, max_weight;
    n = G->n;
    m = G->m;
    src = new vertex_id_t[m];
    assert(src != NULL);
    dest = new vertex_id_t[m];
    assert(dest != NULL);
    degree = new unsigned[n];
    assert(degree != NULL);
    memset(degree, 0, sizeof(unsigned) * n);

    dbl_weight = (double *) malloc(m * sizeof(double));
    assert(dbl_weight != NULL);

    srand48(2387);

    /* generate edges */
    for (edge_id_t i = 0; i < m; i++) {
        vertex_id_t u = rand() % n;
        vertex_id_t v = rand() % n;
        src[i] = u;
        dest[i] = v;
    }

    /* reshuffle */
    if (permute_vertices) {
        srand48(4791);
        permV = new vertex_id_t[n];
        assert(permV != NULL);

        for (vertex_id_t i = 0; i < n; i++) {
            permV[i] = i;
        }

        for (vertex_id_t i = 0; i < n; i++) {
            vertex_id_t j = n * drand48();
            tmpVal = permV[i];
            permV[i] = permV[j];
            permV[j] = tmpVal;
        }

        for (edge_id_t i = 0; i < m; i++) {
            src[i] = permV[src[i]];
            dest[i] = permV[dest[i]];
        }

        delete[] permV;
    }

    min_weight = G->min_weight;
    max_weight = G->max_weight;

    /* Generate edge weights */
    generate_unique_weights(G->m, dbl_weight, G->min_weight, G->max_weight,0);

    // for (edge_id_t i=0; i<m; i++) {
    //     dbl_weight[i]  = min_weight + (max_weight-min_weight)*drand48();
    // }

    /* update graph data structure */
    for (edge_id_t i = 0; i < m; i++) {
        degree[src[i]]++;
        degree[dest[i]]++;
    }

    G->endV = new vertex_id_t[2 * m];
    assert(G->endV != NULL);

    G->rowsIndices = new edge_id_t[n + 1];
    assert(G->rowsIndices != NULL);

    G->n = n;
    /* undirected graph, each edge is stored twice; if edge is (u, v), then it's
     * stored at the vertex u and at the vertex v */
    G->m = 2 * m;

    G->weights = (double *) malloc(G->m * sizeof(double));
    assert(G->weights != NULL);

    G->rowsIndices[0] = 0;
    for (vertex_id_t i = 1; i <= G->n; i++) {
        G->rowsIndices[i] = G->rowsIndices[i - 1] + degree[i - 1];
    }

    for (edge_id_t i = 0; i < m; i++) {
        u = src[i];
        v = dest[i];
        offset = degree[u]--;
        G->endV[G->rowsIndices[u] + offset - 1] = v;
        G->weights[G->rowsIndices[u]+offset-1] = dbl_weight[i];
        offset = degree[v]--;
        G->endV[G->rowsIndices[v] + offset - 1] = u;
        G->weights[G->rowsIndices[v]+offset-1] = dbl_weight[i];

    }

    delete[] src;
    delete[] dest;
    delete[] degree;
}


int main(int argc, char **argv) {
    graph_t g;

    init(argc, argv, &g);
    gen_random_graph(&g);
    // printGraph(&g);
    writeBinaryGraph(&g, outFilename);
    freeGraph(&g);
    return 0;
}
