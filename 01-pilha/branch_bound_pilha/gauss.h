#ifndef ANT_SYSTEM_H_INCLUDED
#define ANT_SYSTEM_H_INCLUDED

#include <stdbool.h>
#define TAMANHO 500 // Não sei porquê, mas não aceita valor maior que 80 !!!
#define INFINITO 900000

struct str_vertice{
    float valor_Z;
    float vetor_x[TAMANHO];
    int pai;
    int no;
};
typedef struct str_vertice vertice;

struct str_problema{
    float matriz_A[TAMANHO][TAMANHO];
    int vetor_relacional[TAMANHO];
    float vetor_b[TAMANHO];
    int minimiza;
    float vetor_custo[TAMANHO];
    int num_linhas;
    int no;
    struct str_problema *proximo; // Só vai saber o que "conteudo" é depois que ler a linha 27 a seguir.
    //struct str_problema *anterior;
};
typedef struct str_problema conteudo;

// COLOCAR CABEÇALHO DE FUNÇÕES AQUI.

// Funções de leitura e escrita de matriz.
void imprime_matriz(float matriz[TAMANHO][TAMANHO], int num_linhas, int num_colunas);
void imprime_vetor(float vetor[TAMANHO], int num_linhas);
void imprime_vetor_int(int vetor[TAMANHO], int num_linhas);
void le_matriz(float matriz[TAMANHO][TAMANHO], int num_linhas, int num_colunas);
void le_vetor(float vetor[TAMANHO], int num_linhas);

// Funções de escrita de matriz em arquivo.
void imprime_matriz_arq(float matriz[TAMANHO][TAMANHO], int num_linhas, int num_colunas, FILE *arqtxt);
void imprime_vetor_arq(float vetor[TAMANHO], int num_linhas, FILE *arqtxt);
void imprime_vetor_int_arq(int vetor[TAMANHO], int num_linhas, FILE *arqtxt);

// Inicializações.
void inicializa_vetor(float vetor[TAMANHO]);
void inicializa_vetor_int(int vetor[TAMANHO]);
void inicializa_matriz(float matriz[TAMANHO][TAMANHO]);

// Funções de manipulação de matrizes e vetores.
void transposta(float matriz[TAMANHO][TAMANHO], float matriz_transposta[TAMANHO][TAMANHO], int num_linhas, int num_colunas);
float multiplica_vetores(float vetor1[TAMANHO], float vetor2[TAMANHO], int num_linhas);
void coluna_matriz(float matriz[TAMANHO][TAMANHO], float vetor[TAMANHO], int num_linhas, int coluna);
int obtem_menor_custo(float vetor[TAMANHO], int num_linhas, float *custo_menor);
void diagonaliza_um(float matriz[TAMANHO][TAMANHO], int num_linhas);

// Eliminação de Gauss e Simplex.
void gauss(float matriz[TAMANHO][TAMANHO], float vetor_b[TAMANHO], float vetor_x[TAMANHO], int num_linhas, int num_colunas, FILE *arqtxt);
int simplex(float matriz_A[TAMANHO][TAMANHO], int vetor_relacional[TAMANHO], float vetor_b[TAMANHO], int minimiza, float vetor_custo[TAMANHO], float vetor_x[TAMANHO], float *valor_Z, int num_linhas, int num_colunas, FILE *arqtxt);

// Funções de cópia de vetores e matrizes.
void copia_vetor_int(int vetor_A[TAMANHO], int vetor_B[TAMANHO], int num_linhas);
void copia_vetor_float(float vetor_A[TAMANHO], float vetor_B[TAMANHO], int num_linhas);
void copia_matriz(float matriz_A[TAMANHO][TAMANHO], float matriz_B[TAMANHO][TAMANHO], int num_linhas, int num_colunas);

// Funções de manipulação de pilhas dinâmicas duplamente ligadas.
void cria_pilha_vazia(conteudo *pilha);
void empilha(conteudo *umVertice, conteudo *pilha);
conteudo *copia_topo_pilha(conteudo *pilha);
conteudo *desempilha(conteudo *pilha);
bool pilha_vazia(conteudo *pilha);

// Branch & Bound
int localiza_real(float vetor[TAMANHO], int num_colunas, int *real);
int branch_bound(float matriz_A[TAMANHO][TAMANHO], int vetor_relacional[TAMANHO], float vetor_b[TAMANHO], int minimiza, float vetor_custo[TAMANHO], vertice *resultados[TAMANHO], int num_linhas, int num_colunas, vertice *sol_otima, FILE *arqtxt);

#endif // ANT_SYSTEM_H_INCLUDED
