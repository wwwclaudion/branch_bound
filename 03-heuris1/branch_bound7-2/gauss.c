#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include "gauss.h"

// Fun��es de leitura e escrita de matriz.

void imprime_matriz(float matriz[TAMANHO][TAMANHO], int num_linhas, int num_colunas){
    int i, j;
    for(i = 0; i < num_linhas; i++){
        for(j = 0; j < num_colunas; j++){
            printf(" %4.2f", matriz[i][j]);
        }
        printf("\n");
    }
}

void imprime_vetor(float vetor[TAMANHO], int num_linhas){
    int i;
    for(i = 0; i < num_linhas; i++){
        printf(" %4.2f \n", vetor[i]);
    }
}

void imprime_vetor_int(int vetor[TAMANHO], int num_linhas){
    int i;
    for(i = 0; i < num_linhas; i++){
        printf(" %4d \n", vetor[i]);
    }
}

void le_matriz(float matriz[TAMANHO][TAMANHO], int num_linhas, int num_colunas){
    int i,j;
    for(i = 0; i < num_linhas; i++){
        for(j = 0; j < num_colunas; j++){
            printf("Digite o valor da posicao A[%d][%d]=", i + 1, j + 1);
            scanf("%f",&matriz[i][j]);
        }
    }
}

void le_vetor(float vetor[TAMANHO], int num_linhas){
    int i;
    for(i = 0; i < num_linhas; i++){
        printf("Digite o valor da posicao vetor[%d]=", i + 1);
        scanf("%f",&vetor[i]);
    }
}

// Fun��es de escrita de matrizes em arquivos.

void imprime_matriz_arq(float matriz[TAMANHO][TAMANHO], int num_linhas, int num_colunas, FILE *arqtxt){
    int i, j;
    for(i = 0; i < num_linhas; i++){
        for(j = 0; j < num_colunas; j++){
            fprintf(arqtxt, " %4.2f", matriz[i][j]);
        }
        fprintf(arqtxt, "\n");
    }
}

void imprime_vetor_arq(float vetor[TAMANHO], int num_linhas, FILE *arqtxt){
    int i;
    for(i = 0; i < num_linhas; i++){
        fprintf(arqtxt, " %4.2f \n", vetor[i]);
    }
}

void imprime_vetor_int_arq(int vetor[TAMANHO], int num_linhas, FILE *arqtxt){
    int i;
    for(i = 0; i < num_linhas; i++){
        fprintf(arqtxt, " %4d \n", vetor[i]);
    }
}

// Inicializa��es.
void inicializa_vetor(float vetor[TAMANHO]){
	int i;
	for(i = 0; i < TAMANHO; i++){
		vetor[i] = 0;
	}
}

void inicializa_vetor_int(int vetor[TAMANHO]){
	int i;
	for(i = 0; i < TAMANHO; i++){
		vetor[i] = 0;
	}
}

void inicializa_matriz(float matriz[TAMANHO][TAMANHO]){
	int i, j;
	for(i = 0; i < TAMANHO; i++){
		for(j = 0; j < TAMANHO; j++){
			matriz[i][j] = 0;
		}
	}
}

// Transposi��o de matriz.
void transposta(float matriz[TAMANHO][TAMANHO], float matriz_transposta[TAMANHO][TAMANHO], int num_linhas, int num_colunas){
    int i, j;
    for(i = 0; i < num_linhas; i++){
        for(j = 0; j < num_colunas; j++){
            matriz_transposta[j][i] = matriz[i][j];
        }
    }
}

// Multiplica vetores.
float multiplica_vetores(float vetor1[TAMANHO], float vetor2[TAMANHO], int num_linhas){
    int i;
    float resultado = 0;
    for(i = 0; i < num_linhas; i++){
        resultado = resultado + vetor1[i] * vetor2[i];
    }
    return resultado;
}

// Obten��o da coluna j de uma matriz:
void coluna_matriz(float matriz[TAMANHO][TAMANHO], float vetor[TAMANHO], int num_linhas, int coluna){
    int i;
    for(i = 0; i < num_linhas; i++){
        vetor[i] = matriz[i][coluna];
    }
}

// Localiza��o da vari�vel de custo m�nimo.
int obtem_menor_custo(float vetor[TAMANHO], int num_linhas, float *custo_menor){
    int i;
    *custo_menor = INFINITO;
    int local;
    for(i = 0; i < num_linhas; i++){
        if(*custo_menor > vetor[i]){
            *custo_menor = vetor[i];
            local = i;
        }
    }
    return local;
}

// Diagonaliza uma matriz inicializada com 1.
void diagonaliza_um(float matriz[TAMANHO][TAMANHO], int num_linhas){
	int i;
	for(i = 0; i < num_linhas; i++){
		matriz[i][i] = 1;
	}
}

// Elimina��o de Gauss
void gauss(float matriz[TAMANHO][TAMANHO], float vetor_b[TAMANHO], float vetor_x[TAMANHO], int num_linhas, int num_colunas, FILE *arqtxt){ // num_colunas da matriz A !
	int i, j, k;
	//int solucao_unica;
	float subtrativo;
	float matriz_expandida[TAMANHO][TAMANHO]; // Esta j� vai se transformar na matriz U.
	float matriz_L[TAMANHO][TAMANHO];
	//float matriz_U[TAMANHO][TAMANHO];

	//solucao_unica = 0;

	// Constru��o da matriz expandida.
	inicializa_matriz(matriz_expandida);
	for(i = 0; i < num_linhas; i++){
		for(j = 0; j < num_colunas; j++){
			matriz_expandida[i][j] = matriz[i][j];
		}
		matriz_expandida[i][j] = vetor_b[i];
	}

	//printf("Matriz expandida: \n");
	//imprime_matriz(matriz_expandida, num_linhas, num_colunas + 1);

	// Inicializa matrizes L.
	inicializa_matriz(matriz_L);
	diagonaliza_um(matriz_L, num_linhas);

	// Elimina��o de Gauss (e preenchimento da matriz L).
	for(j = 0; j < num_colunas + 1; j++){
		for(i = j + 1; i < num_linhas; i++){
            // Consistindo quando pivo = 0: troca pela pr�xima linha at� encontrar (ii) linha que n�o possua 0 nessa posi��o.
			int ii = 1;
			while(matriz_expandida[j][j] == 0 && (j + ii) < num_linhas){
                float temp; // ATEN��O NISTO AQUI !!! Estava int antes, e isto estava lascando com tudo !!!
                for(k = j; k < num_colunas + 1; k++){
                    temp = matriz_expandida[j][k];
                    matriz_expandida[j][k] = matriz_expandida[j + ii][k];
                    matriz_expandida[j + ii][k] = temp;
                }
                ii++;
			}
			if(matriz_expandida[j][j] == 0){ // ESTAVA matriz aqui antes, outro errinho tosco.
                printf("\n\nSegundo Gauss, este sistema possui infinitas solucoes.\n\n");
                fprintf(arqtxt, "\n\nSegundo Gauss, este sistema possui infinitas solucoes.\n\n");
			}
			else{
                //solucao_unica = 1;
                // Agora usa piv� n�o nulo, consistido antes, para obter multiplicador.
                matriz_L[i][j] = matriz_expandida[i][j] / matriz_expandida[j][j];
                if(matriz_expandida[i][j] != 0){ // Se j� for zero, n�o precisa fazer nada.
                    for(k = j; k < num_colunas + 1; k++){
                        matriz_expandida[i][k] = matriz_expandida[i][k] - matriz_L[i][j] * matriz_expandida[j][k]; // j � coluna, mas neste caso corresponde a linha base.
                    }
                }
			}
		}
	}

	//printf("\n\nMatriz expandida apos eliminacaoo de Gauss: \n");
	//imprime_matriz(matriz_expandida, num_linhas, num_colunas + 1);

	// Isolar x de Ax = b.
    for(i = num_linhas - 1; i >= 0; i--){
        subtrativo = 0;
        for(j = i + 1; j < num_colunas; j++){
            subtrativo = subtrativo + vetor_x[j] * matriz_expandida[i][j];
            //printf("\nsubtrativo = %4.2f, vetor_x[%d] = %4.2f, matriz_expandida[%d][%d] = %4.2f\n", subtrativo, j, vetor_x[j], i, j, matriz_expandida[i][j]);
        }
        vetor_x[i] = (matriz_expandida[i][num_colunas] - subtrativo) / matriz_expandida[i][i];
        //printf("\nmatriz_expandida[%d][%d] = %4.2f, subtrativo = %4.2f, matriz_expandida[%d][%d] = %4.2f", i, num_colunas + 1, matriz_expandida[i][num_colunas + 1], subtrativo, i, i, matriz_expandida[i][i]);
        //printf("\nVetor[%d] = %4.2f\n\n", i, vetor_x[i]);
    }
	//return solucao_unica;
}

// Simplex
// num_linhas = n�mero de restri��es.
// num_colunas = n�mero de vari�veis.
int simplex(float matriz_A[TAMANHO][TAMANHO], int vetor_relacional[TAMANHO], float vetor_b[TAMANHO], int minimiza, float vetor_custo[TAMANHO], float vetor_x[TAMANHO], float *valor_Z, int num_linhas, int num_colunas, FILE *arqtxt){
    int i, j, k, n, it1, it2, temp;
    float custo_minimo = (float)INFINITO * (-1);
    int num_colunas_mp, num_variaveis; // n�mero de colunas da matriz padronizada.
    float epsilon;

    float matriz_padronizada[TAMANHO][TAMANHO];
    float matriz_base[TAMANHO][TAMANHO];
    float matriz_base_t[TAMANHO][TAMANHO];

    float vetor_x_temp[TAMANHO];
    float vetor_lambda[TAMANHO];
    float vetor_custo_base[TAMANHO];
    float vetor_custos_relativos[TAMANHO];
    float vetor_colunar[TAMANHO];
    float vetor_y[TAMANHO];

    int indica_fase1 = 0;
    int nao_entra_fase2 = 0;
    int infactibilidade = 0;

    inicializa_matriz(matriz_padronizada);
    inicializa_matriz(matriz_base);
    inicializa_matriz(matriz_base_t);

    inicializa_vetor(vetor_x_temp);
    inicializa_vetor(vetor_lambda);
    inicializa_vetor(vetor_custo_base);
    inicializa_vetor(vetor_custos_relativos);
    inicializa_vetor(vetor_colunar);
    inicializa_vetor(vetor_y);

    // Controladores posicionais.
    int vetor_variaveis[TAMANHO];
    int vetor_base[TAMANHO];

    // Ajuste dos b[i] negativos.
    for(i = 0; i < num_linhas; i++){
        if(vetor_b[i] < 0){
            for(j = 0; j < num_colunas; j++){
                matriz_A[i][j] = (-1) * matriz_A[i][j];
            }
            vetor_relacional[i] = (-1) * vetor_relacional[i];
            vetor_b[i] = (-1) * vetor_b[i];
        }
    }
    printf("\nVetor relacional: \n");
    imprime_vetor_int(vetor_relacional, num_linhas);
    printf("\nVetor b: \n");
    imprime_vetor(vetor_b, num_linhas);

    fprintf(arqtxt, "\nVetor relacional: \n");
    imprime_vetor_int_arq(vetor_relacional, num_linhas, arqtxt);
    fprintf(arqtxt, "\nVetor b: \n");
    imprime_vetor_arq(vetor_b, num_linhas, arqtxt);

    // Constru��o da matriz na forma padr�o.
    num_colunas_mp = num_linhas + num_colunas; // num_colunas_mp = n�mero de colunas na matriz padronizada.
    num_variaveis = num_colunas;
    for(i = 0; i < num_linhas; i++){
        for(j = 0; j < num_colunas; j++){
            matriz_padronizada[i][j] = matriz_A[i][j];
            vetor_variaveis[j] = j; // Vai ser �til depois pra saber qual coluna trocar por qual.
        }
        if(vetor_relacional[i] == -1){ // menor ou igual que.
            matriz_padronizada[i][j + i] = 1;
            vetor_base[i] = j + i;
        }
        else{
            if(vetor_relacional[i] == 1){ // maior ou igual que.
                matriz_padronizada[i][j + i] = -1;
                vetor_variaveis[num_variaveis] = j + i;
                num_variaveis++;
                matriz_padronizada[i][num_colunas_mp] = 1;
                vetor_base[i] = num_colunas_mp;
                num_colunas_mp++;
                indica_fase1 = 1;
            }
            else{
                if(vetor_relacional[i] == 0){
                    matriz_padronizada[i][num_colunas_mp] = 1;
                    vetor_base[i] = num_colunas_mp;
                    indica_fase1 = 1;
                }
            }
        }
    }
    printf("\nNumero de colunas aumentado = %d\n\n", num_colunas_mp);
    printf("\n\nMatriz expandida:\n\n");
    imprime_matriz(matriz_padronizada, num_linhas, num_colunas_mp);
    printf("\nVetor Base:\n");
    imprime_vetor_int(vetor_base, num_linhas);
    printf("\nNumero de variaveis = %d\n\n", num_variaveis);
    printf("\n\nVetor Variaveis:\n");
    imprime_vetor_int(vetor_variaveis, num_variaveis);

    // Fase 1, se necess�rio.
    if(indica_fase1 == 1){
        it1 = 0;

        // Constru��o do vetor de custo artificial.
        float vetor_custo_artificial[TAMANHO];
        inicializa_vetor(vetor_custo_artificial);

        for(i = num_linhas + num_colunas; i < num_colunas_mp; i++){
            vetor_custo_artificial[i] = 1;
        }
        printf("\n\nVetor Custo Artificial:\n");
        imprime_vetor(vetor_custo_artificial, num_colunas_mp);

        while(custo_minimo < 0 && num_colunas_mp > num_linhas + num_colunas){
            printf("\n********** ITERACAO %d da FASE 1 **********\n", it1 + 1);
            // Passo 1.
            // Obten��o da matriz base:
            for(j = 0; j < num_linhas; j++){
                for(i = 0; i < num_linhas; i++){
                    matriz_base[i][j] = matriz_padronizada[i][vetor_base[j]];
                }
            }
            printf("\n\nMatriz Base:\n");
            imprime_matriz(matriz_base, num_linhas, num_linhas);

            // Obten��o do vetor x da base.
            gauss(matriz_base, vetor_b, vetor_x_temp, num_linhas, num_linhas, arqtxt);
            printf("\n\nVetor x da base:\n");
            imprime_vetor(vetor_x_temp, num_linhas);

            // Passo 2.
            // Obten��o do vetor custo dos x da base:
            for(i = 0; i < num_linhas; i++){
                vetor_custo_base[i] = vetor_custo_artificial[vetor_base[i]];
            }
            printf("\n\nVetor custo da base:\n");
            imprime_vetor(vetor_custo_base, num_linhas);

            // Obten��o da transposta da matriz base.
            transposta(matriz_base, matriz_base_t, num_linhas, num_linhas);
            printf("\n\nMatriz Base Transposta:\n");
            imprime_matriz(matriz_base_t, num_linhas, num_linhas);

            // Passo 2.1: Obten��o do vetor lambda.
            gauss(matriz_base_t, vetor_custo_base, vetor_lambda, num_linhas, num_linhas, arqtxt);
            printf("\n\nVetor lambda:\n");
            imprime_vetor(vetor_lambda, num_linhas);

            // Passo 2.2: Obten��o dos custos relativos.
            for(i = 0; i < num_colunas_mp - num_linhas; i++){
                coluna_matriz(matriz_padronizada, vetor_colunar, num_linhas, vetor_variaveis[i]);
                vetor_custos_relativos[i] = vetor_custo_artificial[vetor_variaveis[i]] - multiplica_vetores(vetor_lambda, vetor_colunar, num_linhas);
            }
            printf("\n\nVetor custos relativos:\n");
            imprime_vetor(vetor_custos_relativos, num_colunas_mp - num_linhas);

            // Passo 2.3: Localiza��o do candidato a entrar na matriz_base.
            k = obtem_menor_custo(vetor_custos_relativos, num_colunas_mp - num_linhas, &custo_minimo);
            printf("\nk = %d\n", k);
            printf("custo minimo = %4.2f\n", custo_minimo);

            // Pssso 3.
            if(custo_minimo < 0){
            // Passo 4.
                // Passo 4.1: Dire��o Simplex.
                coluna_matriz(matriz_padronizada, vetor_colunar, num_linhas, vetor_variaveis[k]);
                gauss(matriz_base, vetor_colunar, vetor_y, num_linhas, num_linhas, arqtxt);
                printf("\n\nVetor y:\n");
                imprime_vetor(vetor_y, num_linhas);

                epsilon = INFINITO;
                n = -1;
                // Passo 4.2: C�lculo de �psilon.

                for(i = 0; i < num_linhas; i++){
                    if(vetor_y[i] > 0){
                        if(epsilon > (vetor_x_temp[i] / vetor_y[i])){
                            epsilon = vetor_x_temp[i] / vetor_y[i];
                            n = i; // endere�o no vetor x da base, do x que retorna menor �psilon.
                        }
                    }
                }
                printf("\nepsilon = %4.2f, n = %d\n", epsilon, n);

                if(n > -1){
                    // Pssso 4.3: Atualiza��o do vetor base e do vetor das vari�veis.
                    temp = vetor_base[n];
                    vetor_base[n] = vetor_variaveis[k];
                    if(temp < num_linhas + num_colunas){ // Ou seja, temp n�o recebeu vari�vel artificial.
                        vetor_variaveis[k] = temp;
                    }
                    else{ // temp recebeu vari�vel artificial e ser� descartado, com a matriz base sendo atualizada.
                        for(i = k; i < num_colunas_mp - num_linhas; i++){
                            vetor_variaveis[i] = vetor_variaveis[i + 1];
                        }
                        num_colunas_mp--;
                        num_variaveis--;
                    }
                    it1++;
                    printf("\nATENCAO\n");
                    printf("\nNumero de colunas aumentado = %d\n\n", num_colunas_mp);
                    printf("\nVetor Base:\n");
                    imprime_vetor_int(vetor_base, num_linhas);
                    printf("\nNumero de variaveis = %d\n\n", num_variaveis);
                    printf("\n\nVetor Variaveis:\n");
                    imprime_vetor_int(vetor_variaveis, num_variaveis);
                }
                else{
                    printf("\nPROBLEMA INFACT�VEL !\n");
                    custo_minimo = 1; // For�a sa�da do la�o.
                    num_colunas_mp = num_linhas + num_colunas; // For�a sa�da do la�o.
                    nao_entra_fase2 = 1; // N�o vai entrar na Fase 2.
                    infactibilidade = 1;
                }
            }
            else{
                if(num_colunas_mp > num_linhas + num_colunas){
                    printf("\nPROBLEMA INFACT�VEL !\n");
                    num_colunas_mp = num_linhas + num_colunas; // For�a sa�da do la�o.
                    nao_entra_fase2 = 1; // N�o vai entrar na Fase 2.
                    infactibilidade = 1;
                }
            }
        }
    }

    if(nao_entra_fase2 == 0){
        // Fase 2.
        // Atualiza��o do vetor de custos, que ser� utilizado para a posterior obten��o do vetor de custos relativos.
        float vetor_custo_padronizado[TAMANHO];
        inicializa_vetor(vetor_custo_padronizado);
        it2 = 0;

        // A flag minimiza � bin�ria: 0 se Z maximiza, 1 se Z minimiza.
        // Se for Max f(x), todo o vetor_custo � multiplicado por -1.
        if(minimiza == 0){
            for(j = 0; j < num_colunas; j++){
                vetor_custo_padronizado[j] = (-1) * vetor_custo[j];
            }
        }
        else{
            if(minimiza == 1){
                for(j = 0; j < num_colunas; j++){
                    vetor_custo_padronizado[j] = vetor_custo[j];
                }
            }
            else{
                printf("\nFlag minimiza preenchida de forma errada.\n");
            }
        }
        printf("\n\nVetor custos padronizado:\n");
        imprime_vetor(vetor_custo_padronizado, num_colunas_mp);

        while(custo_minimo < 0){
            printf("\n********** ITERACAO %d da FASE 2 **********\n", it2 + 1);
            // Passo 1.
            // Obten��o da matriz base:
            for(j = 0; j < num_linhas; j++){
                for(i = 0; i < num_linhas; i++){
                    matriz_base[i][j] = matriz_padronizada[i][vetor_base[j]];
                }
            }
            printf("\n\nMatriz Base:\n");
            imprime_matriz(matriz_base, num_linhas, num_linhas);

            // Obten��o do vetor x da base.
            gauss(matriz_base, vetor_b, vetor_x_temp, num_linhas, num_linhas, arqtxt);
            printf("\n\nVetor x da base:\n");
            imprime_vetor(vetor_x_temp, num_linhas);

            // Passo 2.
            // Obten��o do vetor custo dos x da base:
            for(i = 0; i < num_linhas; i++){
                vetor_custo_base[i] = vetor_custo_padronizado[vetor_base[i]];
            }
            printf("\n\nVetor custo da base:\n");
            imprime_vetor(vetor_custo_base, num_linhas);

            // Obten��o da transposta da matriz base.
            transposta(matriz_base, matriz_base_t, num_linhas, num_linhas);

            // Passo 2.1: Obten��o do vetor lambda.
            gauss(matriz_base_t, vetor_custo_base, vetor_lambda, num_linhas, num_linhas, arqtxt);
            printf("\n\nVetor lambda:\n");
            imprime_vetor(vetor_lambda, num_linhas);

            // Passo 2.2: Obten��o dos custos relativos.
            for(i = 0; i < num_colunas_mp - num_linhas; i++){
                coluna_matriz(matriz_padronizada, vetor_colunar, num_linhas, vetor_variaveis[i]);
                vetor_custos_relativos[i] = vetor_custo_padronizado[vetor_variaveis[i]] - multiplica_vetores(vetor_lambda, vetor_colunar, num_linhas);
            }
            printf("\n\nVetor custos relativos:\n");
            imprime_vetor(vetor_custos_relativos, num_colunas_mp - num_linhas);

            // Passo 2.3: Localiza��o do candidato a entrar na matriz_base.
            k = obtem_menor_custo(vetor_custos_relativos, num_colunas_mp - num_linhas, &custo_minimo);
            printf("\nk = %d\n", k);
            printf("custo minimo = %4.2f\n", custo_minimo);

            // Pssso 3.
            if(custo_minimo < 0){
            // Passo 4.
                // Passo 4.1: Dire��o Simplex.
                coluna_matriz(matriz_padronizada, vetor_colunar, num_linhas, vetor_variaveis[k]);
                gauss(matriz_base, vetor_colunar, vetor_y, num_linhas, num_linhas, arqtxt);
                printf("\n\nVetor y:\n");
                imprime_vetor(vetor_y, num_linhas);

                epsilon = INFINITO;
                n = -1;
                // Passo 4.2: C�lculo de �psilon.
                for(i = 0; i < num_linhas; i++){
                    if(vetor_y[i] > 0){
                        if(epsilon > (vetor_x_temp[i] / vetor_y[i])){
                            epsilon = vetor_x_temp[i] / vetor_y[i];
                            n = i; // endere�o no vetor x da base, do x que retorna menor �psilon.
                        }
                    }
                }
                if(n > -1){
                    printf("\nepsilon = %4.2f, n = %d\n", epsilon, n);

                    // Pssso 4.3: Atualiza��o do vetor base e do vetor das vari�veis.
                    temp = vetor_base[n];
                    vetor_base[n] = vetor_variaveis[k];
                    vetor_variaveis[k] = temp;

                    it2++;
                    printf("\nATENCAO\n");
                    printf("\nNumero de colunas aumentado = %d\n\n", num_colunas_mp);
                    printf("\nVetor Base:\n");
                    imprime_vetor_int(vetor_base, num_linhas);
                    printf("\nNumero de variaveis = %d\n\n", num_variaveis);
                    printf("\n\nVetor Variaveis:\n");
                    imprime_vetor_int(vetor_variaveis, num_variaveis);
                }
                else{
                    printf("\nSOLUCAO ILIMITADA !\n");
                    custo_minimo = 1;
                }
            }
        }

        // Preenchimento do vetor_x e de Z.
        for(i = 0; i < num_linhas; i++){
            vetor_x[vetor_base[i]] = vetor_x_temp[i];
        }

        *valor_Z = multiplica_vetores(vetor_x, vetor_custo_padronizado, num_linhas + num_colunas);
        if(minimiza == 0){
            *valor_Z = (-1) * (*valor_Z);
        }
    }

    fprintf(arqtxt, "\nO vetor x eh:\n\n");
    imprime_vetor_arq(vetor_x, num_linhas + num_colunas, arqtxt);

    if(minimiza == 0){
        fprintf(arqtxt, "\nMax f(x) = %4.2f\n\n", *valor_Z);
	}
	else{
        fprintf(arqtxt, "\nMin f(x) = %4.2f\n\n", *valor_Z);
	}
	return infactibilidade;
}
