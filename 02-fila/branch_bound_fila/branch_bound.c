#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include "gauss.h"

// Copia vetores de inteiros
void copia_vetor_int(int vetor_A[TAMANHO], int vetor_B[TAMANHO], int num_linhas){
    int i;
    for(i = 0; i < num_linhas; i++){
        vetor_A[i] = vetor_B[i];
    }
}

void copia_vetor_float(float vetor_A[TAMANHO], float vetor_B[TAMANHO], int num_linhas){
    int i;
    for(i = 0; i < num_linhas; i++){
        vetor_A[i] = vetor_B[i];
    }
}

// Copia matrizes
void copia_matriz(float matriz_A[TAMANHO][TAMANHO], float matriz_B[TAMANHO][TAMANHO], int num_linhas, int num_colunas){
    int i, j;
    for(i = 0; i < num_linhas; i++){
        for(j = 0; j < num_colunas; j++){
            matriz_A[i][j] = matriz_B[i][j];
        }
    }
}

// FUNÇÕES DE MANIPULAÇÃO DE PILHA (lista duplamente ligada, pilha dinâmica):

void cria_pilha_vazia(conteudo *pilha){
    pilha = NULL;
}

void empilha(conteudo *umVertice, conteudo *pilha){
    if(pilha == NULL){
        pilha = (conteudo *) malloc(sizeof(conteudo));
        *pilha = *umVertice;
    }
    else{
        if(umVertice->proximo == NULL){
            umVertice->proximo = (conteudo *) malloc(sizeof(conteudo));
        }
        *(umVertice->proximo) = *pilha; // Atribuição de conteúdo de um ponteiro para outro é sempre assim !
        //pilha->anterior = umVertice;
        *pilha = *umVertice;
    }
}

conteudo *copia_topo_pilha(conteudo *pilha){
    int i = 0, valor;
    conteudo *umConteudo = (conteudo *) malloc(sizeof(conteudo));
    *umConteudo = *pilha;
    if(umConteudo == NULL){
        printf("\nPilha esta vazia.\n");
    }
    return umConteudo;
}

conteudo *desempilha(conteudo *pilha){
    conteudo *umConteudo = (conteudo *) malloc(sizeof(conteudo));
    *umConteudo = *pilha;
    if(umConteudo != NULL){
        *pilha = *(umConteudo->proximo); // Atribuição de conteúdo de um ponteiro para outro é sempre assim !
        //if(pilha != NULL){
            //pilha->anterior = NULL;
        //}
    }
    else{
        printf("\nPilha esta vazia.\n");
    }
    return umConteudo;
}

bool pilha_vazia(conteudo *pilha){
    return pilha == NULL;
}

// FUNÇÕES DE MANIPULAÇÃO DE FILA (fila dinâmica elegante):

void enfileira(fila *umaFila, conteudo *umProblema, int *tamanho_fila){
    conteudo *aux = (conteudo*) malloc(sizeof(conteudo));
    *aux = *umProblema;

    if(*tamanho_fila == 0){ // Caso seja o primeiro elemento, faz com que os ponteiros inicio e fim apontem para ele.
        umaFila->inicio = aux;
        umaFila->fim = aux;
    }

    umaFila->fim->proximo = aux; // faz com que o ponteiro proximo do fim agora aponte para aux.
    umaFila->fim = aux; // Daí, fim recebe aux.
    aux->proximo = NULL; // ponteiro proximo de aux recebe null (ou seja fim da fila).

    (*tamanho_fila)++;
    printf("\n\n\nElemento inserido com sucesso");
}

conteudo *desenfileira(fila *umaFila, int *tamanho_fila){
	conteudo *aux;
	conteudo *resultado = (conteudo*) malloc(sizeof(conteudo));

	if(umaFila->inicio == NULL) {
		printf("\tERRO: Fila vazia");
		resultado = NULL;
	}

	aux = umaFila->inicio;
	*resultado = *aux;
	umaFila->inicio = umaFila->inicio->proximo; // inicio recebe o proximo elemento da fila.
	free(aux);
	(*tamanho_fila)--;
	printf("\n\nElemento removido com sucesso\n");
	return resultado;
}

// Verifica se um vetor contém apenas variáveis inteiras, e devolve endereço da primeira variável que não é inteira.
int localiza_real(float vetor[TAMANHO], int num_colunas, int *real){
    int i;
    float resto;
    i = 0;
    while(i < num_colunas && *real == 0){
        resto = vetor[i] - floor(vetor[i]); // PAU AQUI !!! Devido à flutuações, há ocasiões em que 0.000... não é igual a 0 como deveria ser.
        if(resto > 0.0001){ // Aqui dá pra colocar grau de precisão, tipo maior que 0.0001, por exemplo.
            *real = 1;
        }
        else{
            i++;
        }
    }
    return i;
}

// *solucao_otima e resultados TÊM QUE SEREM INICIALIZADOS ANTES NO MAIN.
// Branch & Bound em profundidade.
int branch_bound(float matriz_A[TAMANHO][TAMANHO], int vetor_relacional[TAMANHO], float vetor_b[TAMANHO], int minimiza, float vetor_custo[TAMANHO], vertice *resultados[TAMANHO], int num_linhas, int num_colunas, vertice *sol_otima, FILE *arqtxt){
    int i, j, variavel, conta_nos, tem_real, num_no, novo_num_linhas, iteracao;
    int infactivel = 0;
    int conta_pilha = 0; // Aqui cumpre o papel de contar nós na fila.
    conteudo *umaFila[TAMANHO];
    conteudo *umProblema = (conteudo *) malloc(sizeof(conteudo));

    num_no = 0;
    conta_nos = 1; // Apenas um controle da quantidade de nós gerados no total, não influencia no funcionamento.
    iteracao = 0;
    (resultados[num_no])->no = conta_nos;
    (resultados[num_no])->pai = num_no;

    copia_matriz(umProblema->matriz_A, matriz_A, num_linhas, num_colunas);
    copia_vetor_int(umProblema->vetor_relacional, vetor_relacional, num_linhas);
    copia_vetor_float(umProblema->vetor_b, vetor_b, num_linhas);
    umProblema->minimiza = minimiza;
    copia_vetor_float(umProblema->vetor_custo, vetor_custo, num_colunas);
    umProblema->num_linhas = num_linhas;
    umProblema->no = conta_nos;
    //umProblema->anterior = NULL;
    umProblema->proximo = NULL;

    if(minimiza == 0){
        sol_otima->valor_Z = -INFINITO;
    }
    else{
        sol_otima->valor_Z = INFINITO;
    }

    enfileira(umaFila, umProblema, &conta_pilha);
    //empilha(umProblema, pilha);
    //conta_pilha++; // Apenas uma bandeira de controle do tamanho da pilha. Não influencia o funcionamento.

    while(conta_pilha > 0){
        umProblema = desenfileira(umaFila, &conta_pilha);
        //umProblema = desempilha(pilha);
        //conta_pilha--;
        num_no = umProblema->no;
        novo_num_linhas = umProblema->num_linhas;

        printf("\n\n**************************ITERACAO %d DO BRANCH & BOUND *******************************\n\n", iteracao);
        fprintf(arqtxt, "\n\n**************************ITERACAO %d, no %d, DO BRANCH & BOUND ********************************\n\n", iteracao, num_no);
        fprintf(arqtxt, "conta_pilha = %d\n\n", conta_pilha);
        imprime_matriz_arq(umProblema->matriz_A, novo_num_linhas, num_colunas, arqtxt);

        infactivel = simplex(umProblema->matriz_A, umProblema->vetor_relacional, umProblema->vetor_b, umProblema->minimiza, umProblema->vetor_custo, (resultados[num_no])->vetor_x, &((resultados[num_no])->valor_Z), novo_num_linhas, num_colunas, arqtxt);
        fprintf(arqtxt, "\n\ninfactivel = %d\n", infactivel);

        if(infactivel == 0){ // Ou seja, tem solução.

            i = 0;
            tem_real = 0;

            variavel = localiza_real((resultados[num_no])->vetor_x, num_colunas, &tem_real);

            fprintf(arqtxt, "\nvariavel = %d, tem_real = %d\n", variavel, tem_real);

            if(tem_real == 0){
                if(minimiza == 0){
                    if((resultados[num_no])->valor_Z > sol_otima->valor_Z){
                        *sol_otima = *(resultados[num_no]); // ISSO DEVE ACONTECER SEMPRE, E NÃO APENAS QUANDO ENCONTRA NÓ CANDIDATO !!!
                    }
                }
                else{
                    if((resultados[num_no])->valor_Z < sol_otima->valor_Z){
                        *sol_otima = *(resultados[num_no]);
                    }
                }
            }
            else{
                // Aqui decide se vale a pena continuar bifurcando.
                fprintf(arqtxt, "\nnum_no = %d, e num_no do pai = %d\n", num_no, (resultados[num_no])->pai);
                if((resultados[num_no])->valor_Z != (resultados[(resultados[num_no])->pai])->valor_Z && (resultados[num_no])->vetor_x[variavel] != (resultados[(resultados[num_no])->pai])->vetor_x[variavel]){
                    if((minimiza == 0 && (resultados[num_no])->valor_Z > sol_otima->valor_Z) || (minimiza == 1 && (resultados[num_no])->valor_Z < sol_otima->valor_Z)){
                        float nova_matriz_A[TAMANHO][TAMANHO];
                        inicializa_matriz(nova_matriz_A);
                        copia_matriz(nova_matriz_A, umProblema->matriz_A, novo_num_linhas, num_colunas);
                        nova_matriz_A[novo_num_linhas][variavel] = 1;

                        // Bifurcação 1
                        int novo_vetor_relacional[TAMANHO];
                        inicializa_vetor_int(novo_vetor_relacional);
                        copia_vetor_int(novo_vetor_relacional, umProblema->vetor_relacional, novo_num_linhas);
                        novo_vetor_relacional[novo_num_linhas] = -1; // menor ou igual.

                        float novo_vetor_b[TAMANHO];
                        inicializa_vetor(novo_vetor_b);
                        copia_vetor_float(novo_vetor_b, umProblema->vetor_b, novo_num_linhas);
                        novo_vetor_b[novo_num_linhas] = floor((resultados[num_no])->vetor_x[variavel]); // AQUI ENTRA O PISO DE X.

                        (resultados[conta_nos + 1])->pai = num_no;
                        (resultados[conta_nos + 1])->no = conta_nos + 1;

                        copia_matriz(umProblema->matriz_A, nova_matriz_A, novo_num_linhas + 1, num_colunas);
                        copia_vetor_int(umProblema->vetor_relacional, novo_vetor_relacional, novo_num_linhas + 1);
                        copia_vetor_float(umProblema->vetor_b, novo_vetor_b, novo_num_linhas + 1);
                        umProblema->minimiza = minimiza;
                        copia_vetor_float(umProblema->vetor_custo, vetor_custo, num_colunas);
                        umProblema->num_linhas = novo_num_linhas + 1;
                        umProblema->no = conta_nos + 1;
                        //umProblema->anterior = NULL;
                        umProblema->proximo = NULL;

                        enfileira(umaFila, umProblema, &conta_pilha);
                        conta_nos++;
                        //conta_pilha++;

                        // Bifurcação 2
                        novo_vetor_relacional[novo_num_linhas] = 1; // maior ou igual.

                        novo_vetor_b[novo_num_linhas] = ceil((resultados[num_no])->vetor_x[variavel]); // AQUI ENTRA O TETO DE X.

                        (resultados[conta_nos + 1])->pai = num_no;
                        (resultados[conta_nos + 1])->no = conta_nos + 1;

                        copia_matriz(umProblema->matriz_A, nova_matriz_A, novo_num_linhas + 1, num_colunas);
                        copia_vetor_int(umProblema->vetor_relacional, novo_vetor_relacional, novo_num_linhas + 1);
                        copia_vetor_float(umProblema->vetor_b, novo_vetor_b, novo_num_linhas + 1);
                        umProblema->minimiza = minimiza;
                        copia_vetor_float(umProblema->vetor_custo, vetor_custo, num_colunas);
                        umProblema->num_linhas = novo_num_linhas + 1;
                        umProblema->no = conta_nos + 1;
                        //umProblema->anterior = NULL;
                        umProblema->proximo = NULL;

                        enfileira(umaFila, umProblema, &conta_pilha);
                        conta_nos++;
                        //conta_pilha++;
                    }
                }
            }
        }
        else{
            printf("\n\nO sub-problema %d nao tem solucao.\n\n", num_no);
            fprintf(arqtxt, "\n\nO sub-problema %d nao tem solucao.\n\n", num_no);
        }
        iteracao++;
    }
    return conta_nos;
}

