// Implementação do algoritmo simplex.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "gauss.h"
#include "time.h"

int main(){
	system("color f9");
	clock_t t;

	float matriz_A[TAMANHO][TAMANHO];
	float vetor_x[TAMANHO];
	float vetor_b[TAMANHO];
	int vetor_relacional[TAMANHO];
	float vetor_custo[TAMANHO];
	//float valor_Z;
	int i, j;
	int total_nos;

	vertice *solucao_otima = (vertice *) malloc(sizeof(vertice));
	vertice *resultados[TAMANHO];
	for(i = 0; i < TAMANHO; i++){
        resultados[i] = (vertice *) malloc(sizeof(vertice));
        (resultados[i])->valor_Z = 0;
        for(j = 0; j < TAMANHO; j++){
            (resultados[i])->vetor_x[j] = 0;
        }
	}

	FILE *arqtxt;

	// Daqui em diante começa o preenchimento do problema.
	// Quarto exemplo http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/mknap1.txt
	// Versão binária // Versão unbounded resulta infinitas soluções.
	int num_linhas = 10; // 10 restrições
	int num_colunas = 80;
    int minimiza = 0; // 1 = Min, 0 = Max

	inicializa_matriz(matriz_A);

    matriz_A[0][0] = 40;   matriz_A[0][1] = 91;   matriz_A[0][2] = 10;   matriz_A[0][3] = 30;   matriz_A[0][4] = 160;  matriz_A[0][5] = 20;   matriz_A[0][6] = 3;    matriz_A[0][7] = 12;   matriz_A[0][8] = 3;    matriz_A[0][9] = 18;   matriz_A[0][10] = 9;   matriz_A[0][11] = 25;  matriz_A[0][12] = 1;   matriz_A[0][13] = 1;   matriz_A[0][14] = 10; matriz_A[0][15] = 280;  matriz_A[0][16] = 10;  matriz_A[0][17] = 8;  matriz_A[0][18] = 1;   matriz_A[0][19] = 1; 
	matriz_A[1][0] = 16;   matriz_A[1][1] = 92;   matriz_A[1][2] = 41;   matriz_A[1][3] = 16;   matriz_A[1][4] = 150;  matriz_A[1][5] = 23;   matriz_A[1][6] = 4;    matriz_A[1][7] = 18;   matriz_A[1][8] = 6;    matriz_A[1][9] = 0;    matriz_A[1][10] = 12;  matriz_A[1][11] = 8;   matriz_A[1][12] = 2;   matriz_A[1][13] = 1;   matriz_A[1][14] = 0;  matriz_A[1][15] = 200;  matriz_A[1][16] = 20;  matriz_A[1][17] = 6;  matriz_A[1][18] = 2;   matriz_A[1][19] = 1; 
	matriz_A[2][0] = 38;   matriz_A[2][1] = 39;   matriz_A[2][2] = 32;   matriz_A[2][3] = 71;   matriz_A[2][4] = 80;   matriz_A[2][5] = 26;   matriz_A[2][6] = 5;    matriz_A[2][7] = 40;   matriz_A[2][8] = 8;    matriz_A[2][9] = 12;   matriz_A[2][10] = 30;  matriz_A[2][11] = 15;  matriz_A[2][12] = 0;   matriz_A[2][13] = 1;   matriz_A[2][14] = 23; matriz_A[2][15] = 100;  matriz_A[2][16] = 0;   matriz_A[2][17] = 20; matriz_A[2][18] = 3;   matriz_A[2][19] = 0; 
	matriz_A[3][0] = 8;    matriz_A[3][1] = 71;   matriz_A[3][2] = 30;   matriz_A[3][3] = 60;   matriz_A[3][4] = 200;  matriz_A[3][5] = 18;   matriz_A[3][6] = 6;    matriz_A[3][7] = 30;   matriz_A[3][8] = 4;    matriz_A[3][9] = 8;    matriz_A[3][10] = 31;  matriz_A[3][11] = 6;   matriz_A[3][12] = 3;   matriz_A[3][13] = 0;   matriz_A[3][14] = 18; matriz_A[3][15] = 60;   matriz_A[3][16] = 21;  matriz_A[3][17] = 4;  matriz_A[3][18] = 0;   matriz_A[3][19] = 2; 
	matriz_A[4][0] = 38;   matriz_A[4][1] = 52;   matriz_A[4][2] = 30;   matriz_A[4][3] = 42;   matriz_A[4][4] = 170;  matriz_A[4][5] = 9;    matriz_A[4][6] = 7;    matriz_A[4][7] = 20;   matriz_A[4][8] = 1;    matriz_A[4][9] = 3;    matriz_A[4][10] = 21;  matriz_A[4][11] = 4;   matriz_A[4][12] = 1;   matriz_A[4][13] = 2;   matriz_A[4][14] = 14; matriz_A[4][15] = 310;  matriz_A[4][16] = 8;   matriz_A[4][17] = 4;  matriz_A[4][18] = 6;   matriz_A[4][19] = 1; 
	matriz_A[5][0] = 5;    matriz_A[5][1] = 10;   matriz_A[5][2] = 5;    matriz_A[5][3] = 15;   matriz_A[5][4] = 91;   matriz_A[5][5] = 24;   matriz_A[5][6] = 10;   matriz_A[5][7] = 15;   matriz_A[5][8] = 90;   matriz_A[5][9] = 15;   matriz_A[5][10] = 60;  matriz_A[5][11] = 5;   matriz_A[5][12] = 55;  matriz_A[5][13] = 60;  matriz_A[5][14] = 50; matriz_A[5][15] = 75;   matriz_A[5][16] = 100; matriz_A[5][17] = 65; matriz_A[5][18] = 15;  matriz_A[5][19] = 10; 
	matriz_A[6][0] = 10;   matriz_A[6][1] = 5;    matriz_A[6][2] = 6;    matriz_A[6][3] = 11;   matriz_A[6][4] = 41;   matriz_A[6][5] = 30;   matriz_A[6][6] = 5;    matriz_A[6][7] = 40;   matriz_A[6][8] = 2;    matriz_A[6][9] = 6;    matriz_A[6][10] = 100; matriz_A[6][11] = 10;  matriz_A[6][12] = 25;  matriz_A[6][13] = 39;  matriz_A[6][14] = 30; matriz_A[6][15] = 13;   matriz_A[6][16] = 30;  matriz_A[6][17] = 15; matriz_A[6][18] = 60;  matriz_A[6][19] = 5; 
	matriz_A[7][0] = 7;    matriz_A[7][1] = 0;    matriz_A[7][2] = 30;   matriz_A[7][3] = 22;   matriz_A[7][4] = 80;   matriz_A[7][5] = 94;   matriz_A[7][6] = 11;   matriz_A[7][7] = 81;   matriz_A[7][8] = 70;   matriz_A[7][9] = 64;   matriz_A[7][10] = 59;  matriz_A[7][11] = 18;  matriz_A[7][12] = 0;   matriz_A[7][13] = 36;  matriz_A[7][14] = 3;  matriz_A[7][15] = 8;    matriz_A[7][16] = 15;  matriz_A[7][17] = 42; matriz_A[7][18] = 9;   matriz_A[7][19] = 0; 
	matriz_A[8][0] = 19;   matriz_A[8][1] = 37;   matriz_A[8][2] = 27;   matriz_A[8][3] = 62;   matriz_A[8][4] = 39;   matriz_A[8][5] = 84;   matriz_A[8][6] = 16;   matriz_A[8][7] = 14;   matriz_A[8][8] = 21;   matriz_A[8][9] = 5;    matriz_A[8][10] = 21;  matriz_A[8][11] = 40;  matriz_A[8][12] = 1;   matriz_A[8][13] = 6;   matriz_A[8][14] = 82; matriz_A[8][15] = 91;   matriz_A[8][16] = 43;  matriz_A[8][17] = 30; matriz_A[8][18] = 62;  matriz_A[8][19] = 91; 
	matriz_A[9][0] = 15;   matriz_A[9][1] = 25;   matriz_A[9][2] = 0;    matriz_A[9][3] = 94;   matriz_A[9][4] = 53;   matriz_A[9][5] = 48;   matriz_A[9][6] = 27;   matriz_A[9][7] = 99;   matriz_A[9][8] = 6;    matriz_A[9][9] = 17;   matriz_A[9][10] = 69;  matriz_A[9][11] = 43;  matriz_A[9][12] = 0;   matriz_A[9][13] = 57;  matriz_A[9][14] = 7;  matriz_A[9][15] = 21;   matriz_A[9][16] = 78;  matriz_A[9][17] = 10; matriz_A[9][18] = 37;  matriz_A[9][19] = 26; 

	matriz_A[0][20] = 49;  matriz_A[0][21] = 8;   matriz_A[0][22] = 21;  matriz_A[0][23] = 6;   matriz_A[0][24] = 1;   matriz_A[0][25] = 5;   matriz_A[0][26] = 10;  matriz_A[0][27] = 8;   matriz_A[0][28] = 2;   matriz_A[0][29] = 1;   matriz_A[0][30] = 0;   matriz_A[0][31] = 10;  matriz_A[0][32] = 42;  matriz_A[0][33] = 6;   matriz_A[0][34] = 4;  matriz_A[0][35] = 8;   matriz_A[0][36] = 0;   matriz_A[0][37] = 10;  matriz_A[0][38] = 1;   matriz_A[0][39] = 40;
	matriz_A[1][20] = 70;  matriz_A[1][21] = 9;   matriz_A[1][22] = 22;  matriz_A[1][23] = 4;   matriz_A[1][24] = 1;   matriz_A[1][25] = 5;   matriz_A[1][26] = 10;  matriz_A[1][27] = 6;   matriz_A[1][28] = 4;   matriz_A[1][29] = 0;   matriz_A[1][30] = 4;   matriz_A[1][31] = 12;  matriz_A[1][32] = 8;   matriz_A[1][33] = 4;   matriz_A[1][34] = 3;  matriz_A[1][35] = 0;   matriz_A[1][36] = 10;  matriz_A[1][37] = 0;   matriz_A[1][38] = 6;   matriz_A[1][39] = 28;
	matriz_A[2][20] = 40;  matriz_A[2][21] = 6;   matriz_A[2][22] = 8;   matriz_A[2][23] = 0;   matriz_A[2][24] = 6;   matriz_A[2][25] = 4;   matriz_A[2][26] = 22;  matriz_A[2][27] = 4;   matriz_A[2][28] = 6;   matriz_A[2][29] = 1;   matriz_A[2][30] = 5;   matriz_A[2][31] = 14;  matriz_A[2][32] = 8;   matriz_A[2][33] = 2;   matriz_A[2][34] = 8;  matriz_A[2][35] = 0;   matriz_A[2][36] = 20;  matriz_A[2][37] = 0;   matriz_A[2][38] = 0;   matriz_A[2][39] = 6;
	matriz_A[3][20] = 32;  matriz_A[3][21] = 15;  matriz_A[3][22] = 31;  matriz_A[3][23] = 2;   matriz_A[3][24] = 2;   matriz_A[3][25] = 7;   matriz_A[3][26] = 8;   matriz_A[3][27] = 2;   matriz_A[3][28] = 8;   matriz_A[3][29] = 0;   matriz_A[3][30] = 2;   matriz_A[3][31] = 8;   matriz_A[3][32] = 6;   matriz_A[3][33] = 7;   matriz_A[3][34] = 1;  matriz_A[3][35] = 0;   matriz_A[3][36] = 0;   matriz_A[3][37] = 20;  matriz_A[3][38] = 8;   matriz_A[3][39] = 14;
	matriz_A[4][20] = 18;  matriz_A[4][21] = 15;  matriz_A[4][22] = 38;  matriz_A[4][23] = 10;  matriz_A[4][24] = 4;   matriz_A[4][25] = 8;   matriz_A[4][26] = 6;   matriz_A[4][27] = 0;   matriz_A[4][28] = 0;   matriz_A[4][29] = 3;   matriz_A[4][30] = 0;   matriz_A[4][31] = 10;  matriz_A[4][32] = 6;   matriz_A[4][33] = 1;   matriz_A[4][34] = 3;  matriz_A[4][35] = 0;   matriz_A[4][36] = 3;   matriz_A[4][37] = 5;   matriz_A[4][38] = 4;   matriz_A[4][39] = 0;
	matriz_A[5][20] = 30;  matriz_A[5][21] = 35;  matriz_A[5][22] = 50;  matriz_A[5][23] = 15;  matriz_A[5][24] = 45;  matriz_A[5][25] = 80;  matriz_A[5][26] = 40;  matriz_A[5][27] = 110; matriz_A[5][28] = 80;  matriz_A[5][29] = 80;  matriz_A[5][30] = 36;  matriz_A[5][31] = 20;  matriz_A[5][32] = 90;  matriz_A[5][33] = 50;  matriz_A[5][34] = 25; matriz_A[5][35] = 50;  matriz_A[5][36] = 35;  matriz_A[5][37] = 30;  matriz_A[5][38] = 60;  matriz_A[5][39] = 10;
	matriz_A[6][20] = 5;   matriz_A[6][21] = 10;  matriz_A[6][22] = 5;   matriz_A[6][23] = 15;  matriz_A[6][24] = 91;  matriz_A[6][25] = 24;  matriz_A[6][26] = 10;  matriz_A[6][27] = 15;  matriz_A[6][28] = 90;  matriz_A[6][29] = 15;  matriz_A[6][30] = 60;  matriz_A[6][31] = 5;   matriz_A[6][32] = 55;  matriz_A[6][33] = 60;  matriz_A[6][34] = 50; matriz_A[6][35] = 75;  matriz_A[6][36] = 100; matriz_A[6][37] = 65;  matriz_A[6][38] = 15;  matriz_A[6][39] = 10;
	matriz_A[7][20] = 42;  matriz_A[7][21] = 47;  matriz_A[7][22] = 52;  matriz_A[7][23] = 32;  matriz_A[7][24] = 26;  matriz_A[7][25] = 48;  matriz_A[7][26] = 55;  matriz_A[7][27] = 6;   matriz_A[7][28] = 29;  matriz_A[7][29] = 84;  matriz_A[7][30] = 8;   matriz_A[7][31] = 66;  matriz_A[7][32] = 98;  matriz_A[7][33] = 50;  matriz_A[7][34] = 0;  matriz_A[7][35] = 30;  matriz_A[7][36] = 0;   matriz_A[7][37] = 88;  matriz_A[7][38] = 15;  matriz_A[7][39] = 37;
	matriz_A[8][20] = 10;  matriz_A[8][21] = 41;  matriz_A[8][22] = 12;  matriz_A[8][23] = 4;   matriz_A[8][24] = 80;  matriz_A[8][25] = 77;  matriz_A[8][26] = 98;  matriz_A[8][27] = 50;  matriz_A[8][28] = 78;  matriz_A[8][29] = 35;  matriz_A[8][30] = 7;   matriz_A[8][31] = 1;   matriz_A[8][32] = 96;  matriz_A[8][33] = 67;  matriz_A[8][34] = 85; matriz_A[8][35] = 4;   matriz_A[8][36] = 23;  matriz_A[8][37] = 38;  matriz_A[8][38] = 2;   matriz_A[8][39] = 57;
	matriz_A[9][20] = 8;   matriz_A[9][21] = 66;  matriz_A[9][22] = 98;  matriz_A[9][23] = 50;  matriz_A[9][24] = 0;   matriz_A[9][25] = 30;  matriz_A[9][26] = 0;   matriz_A[9][27] = 88;  matriz_A[9][28] = 15;  matriz_A[9][29] = 37;  matriz_A[9][30] = 26;  matriz_A[9][31] = 72;  matriz_A[9][32] = 61;  matriz_A[9][33] = 57;  matriz_A[9][34] = 17; matriz_A[9][35] = 27;  matriz_A[9][36] = 83;  matriz_A[9][37] = 3;   matriz_A[9][38] = 9;   matriz_A[9][39] = 66;

    matriz_A[0][40] = 86;  matriz_A[0][41] = 11;  matriz_A[0][42] = 120; matriz_A[0][43] = 8;   matriz_A[0][44] = 3;   matriz_A[0][45] = 32;  matriz_A[0][46] = 28;  matriz_A[0][47] = 13;  matriz_A[0][48] = 2;   matriz_A[0][49] = 4;   matriz_A[0][50] = 45;  matriz_A[0][51] = 0;   matriz_A[0][52] = 85;  matriz_A[0][53] = 150;  matriz_A[0][54] = 65; matriz_A[0][55] = 95;  matriz_A[0][56] = 30;  matriz_A[0][57] = 0;   matriz_A[0][58] = 170; matriz_A[0][59] = 0; 
	matriz_A[1][40] = 93;  matriz_A[1][41] = 9;   matriz_A[1][42] = 30;  matriz_A[1][43] = 22;  matriz_A[1][44] = 0;   matriz_A[1][45] = 36;  matriz_A[1][46] = 45;  matriz_A[1][47] = 13;  matriz_A[1][48] = 2;   matriz_A[1][49] = 2;   matriz_A[1][50] = 75;  matriz_A[1][51] = 40;  matriz_A[1][52] = 365; matriz_A[1][53] = 95;   matriz_A[1][54] = 25; matriz_A[1][55] = 17;  matriz_A[1][56] = 125; matriz_A[1][57] = 20;  matriz_A[1][58] = 22;  matriz_A[1][59] = 84; 
	matriz_A[2][40] = 12;  matriz_A[2][41] = 6;   matriz_A[2][42] = 80;  matriz_A[2][43] = 13;  matriz_A[2][44] = 6;   matriz_A[2][45] = 22;  matriz_A[2][46] = 14;  matriz_A[2][47] = 0;   matriz_A[2][48] = 1;   matriz_A[2][49] = 2;   matriz_A[2][50] = 75;  matriz_A[2][51] = 0;   matriz_A[2][52] = 102; matriz_A[2][53] = 0;    matriz_A[2][54] = 0;  matriz_A[2][55] = 40;  matriz_A[2][56] = 60;  matriz_A[2][57] = 0;   matriz_A[2][58] = 165; matriz_A[2][59] = 0; 
	matriz_A[3][40] = 20;  matriz_A[3][41] = 2;   matriz_A[3][42] = 40;  matriz_A[3][43] = 6;   matriz_A[3][44] = 1;   matriz_A[3][45] = 14;  matriz_A[3][46] = 20;  matriz_A[3][47] = 12;  matriz_A[3][48] = 0;   matriz_A[3][49] = 1;   matriz_A[3][50] = 0;   matriz_A[3][51] = 0;   matriz_A[3][52] = 0;   matriz_A[3][53] = 0;    matriz_A[3][54] = 0;  matriz_A[3][55] = 80;  matriz_A[3][56] = 0;   matriz_A[3][57] = 110; matriz_A[3][58] = 0;   matriz_A[3][59] = 15; 
	matriz_A[4][40] = 30;  matriz_A[4][41] = 12;  matriz_A[4][42] = 16;  matriz_A[4][43] = 18;  matriz_A[4][44] = 3;   matriz_A[4][45] = 6;   matriz_A[4][46] = 22;  matriz_A[4][47] = 30;  matriz_A[4][48] = 4;   matriz_A[4][49] = 0;   matriz_A[4][50] = 0;   matriz_A[4][51] = 0;   matriz_A[4][52] = 0;   matriz_A[4][53] = 5;    matriz_A[4][54] = 10; matriz_A[4][55] = 10;  matriz_A[4][56] = 50;  matriz_A[4][57] = 2;   matriz_A[4][58] = 5;   matriz_A[4][59] = 5; 
	matriz_A[5][40] = 150; matriz_A[5][41] = 110; matriz_A[5][42] = 70;  matriz_A[5][43] = 10;  matriz_A[5][44] = 20;  matriz_A[5][45] = 30;  matriz_A[5][46] = 104; matriz_A[5][47] = 40;  matriz_A[5][48] = 40;  matriz_A[5][49] = 94;  matriz_A[5][50] = 150; matriz_A[5][51] = 50;  matriz_A[5][52] = 10;  matriz_A[5][53] = 50;   matriz_A[5][54] = 50; matriz_A[5][55] = 16;  matriz_A[5][56] = 10;  matriz_A[5][57] = 20;  matriz_A[5][58] = 50;  matriz_A[5][59] = 90; 
	matriz_A[6][40] = 30;  matriz_A[6][41] = 35;  matriz_A[6][42] = 50;  matriz_A[6][43] = 15;  matriz_A[6][44] = 45;  matriz_A[6][45] = 80;  matriz_A[6][46] = 40;  matriz_A[6][47] = 110; matriz_A[6][48] = 80;  matriz_A[6][49] = 80;  matriz_A[6][50] = 36;  matriz_A[6][51] = 20;  matriz_A[6][52] = 90;  matriz_A[6][53] = 50;   matriz_A[6][54] = 25; matriz_A[6][55] = 50;  matriz_A[6][56] = 35;  matriz_A[6][57] = 30;  matriz_A[6][58] = 60;  matriz_A[6][59] = 10; 
	matriz_A[7][40] = 26;  matriz_A[7][41] = 72;  matriz_A[7][42] = 61;  matriz_A[7][43] = 57;  matriz_A[7][44] = 17;  matriz_A[7][45] = 27;  matriz_A[7][46] = 83;  matriz_A[7][47] = 3;   matriz_A[7][48] = 9;   matriz_A[7][49] = 66;  matriz_A[7][50] = 97;  matriz_A[7][51] = 42;  matriz_A[7][52] = 2;   matriz_A[7][53] = 44;   matriz_A[7][54] = 71; matriz_A[7][55] = 11;  matriz_A[7][56] = 25;  matriz_A[7][57] = 74;  matriz_A[7][58] = 90;  matriz_A[7][59] = 20; 
	matriz_A[8][40] = 94;  matriz_A[8][41] = 86;  matriz_A[8][42] = 80;  matriz_A[8][43] = 92;  matriz_A[8][44] = 31;  matriz_A[8][45] = 17;  matriz_A[8][46] = 65;  matriz_A[8][47] = 51;  matriz_A[8][48] = 46;  matriz_A[8][49] = 66;  matriz_A[8][50] = 44;  matriz_A[8][51] = 3;   matriz_A[8][52] = 26;  matriz_A[8][53] = 0;    matriz_A[8][54] = 39; matriz_A[8][55] = 20;  matriz_A[8][56] = 11;  matriz_A[8][57] = 6;   matriz_A[8][58] = 55;  matriz_A[8][59] = 70; 
	matriz_A[9][40] = 97;  matriz_A[9][41] = 42;  matriz_A[9][42] = 2;   matriz_A[9][43] = 44;  matriz_A[9][44] = 71;  matriz_A[9][45] = 11;  matriz_A[9][46] = 25;  matriz_A[9][47] = 74;  matriz_A[9][48] = 90;  matriz_A[9][49] = 20;  matriz_A[9][50] = 0;   matriz_A[9][51] = 38;  matriz_A[9][52] = 33;  matriz_A[9][53] = 14;   matriz_A[9][54] = 9;  matriz_A[9][55] = 23;  matriz_A[9][56] = 12;  matriz_A[9][57] = 58;  matriz_A[9][58] = 6;   matriz_A[9][59] = 14; 

	matriz_A[0][60] = 40;  matriz_A[0][61] = 25;  matriz_A[0][62] = 20;  matriz_A[0][63] = 0;   matriz_A[0][64] = 0;   matriz_A[0][65] = 25;   matriz_A[0][66] = 0;   matriz_A[0][67] = 0;  matriz_A[0][68] = 25;  matriz_A[0][69] = 0;   matriz_A[0][70] = 30;  matriz_A[0][71] = 20;  matriz_A[0][72] = 125; matriz_A[0][73] = 5;   matriz_A[0][74] = 80; matriz_A[0][75] = 25; matriz_A[0][76] = 35;  matriz_A[0][77] = 73;   matriz_A[0][78] = 12;  matriz_A[0][79] = 15;
	matriz_A[1][60] = 75;  matriz_A[1][61] = 50;  matriz_A[1][62] = 15;  matriz_A[1][63] = 0;   matriz_A[1][64] = 0;   matriz_A[1][65] = 12;   matriz_A[1][66] = 0;   matriz_A[1][67] = 10; matriz_A[1][68] = 0;   matriz_A[1][69] = 50;  matriz_A[1][70] = 0;   matriz_A[1][71] = 0;   matriz_A[1][72] = 10;  matriz_A[1][73] = 0;   matriz_A[1][74] = 0;  matriz_A[1][75] = 50; matriz_A[1][76] = 60;  matriz_A[1][77] = 150;  matriz_A[1][78] = 0;   matriz_A[1][79] = 0;
	matriz_A[2][60] = 0;   matriz_A[2][61] = 0;   matriz_A[2][62] = 45;  matriz_A[2][63] = 0;   matriz_A[2][64] = 0;   matriz_A[2][65] = 0;    matriz_A[2][66] = 25;  matriz_A[2][67] = 0;  matriz_A[2][68] = 150; matriz_A[2][69] = 0;   matriz_A[2][70] = 0;   matriz_A[2][71] = 0;   matriz_A[2][72] = 158; matriz_A[2][73] = 0;   matriz_A[2][74] = 85; matriz_A[2][75] = 95; matriz_A[2][76] = 0;   matriz_A[2][77] = 89;   matriz_A[2][78] = 20;  matriz_A[2][79] = 1;
	matriz_A[3][60] = 0;   matriz_A[3][61] = 60;  matriz_A[3][62] = 5;   matriz_A[3][63] = 135; matriz_A[3][64] = 0;   matriz_A[3][65] = 0;    matriz_A[3][66] = 25;  matriz_A[3][67] = 0;  matriz_A[3][68] = 300; matriz_A[3][69] = 35;  matriz_A[3][70] = 100; matriz_A[3][71] = 0;   matriz_A[3][72] = 0;   matriz_A[3][73] = 25;  matriz_A[3][74] = 0;  matriz_A[3][75] = 0;  matriz_A[3][76] = 225; matriz_A[3][77] = 25;   matriz_A[3][78] = 0;   matriz_A[3][79] = 0;
	matriz_A[4][60] = 10;  matriz_A[4][61] = 5;   matriz_A[4][62] = 6;   matriz_A[4][63] = 11;  matriz_A[4][64] = 41;  matriz_A[4][65] = 30;   matriz_A[4][66] = 5;   matriz_A[4][67] = 40; matriz_A[4][68] = 2;   matriz_A[4][69] = 6;   matriz_A[4][70] = 100; matriz_A[4][71] = 10;  matriz_A[4][72] = 25;  matriz_A[4][73] = 39;  matriz_A[4][74] = 30; matriz_A[4][75] = 13; matriz_A[4][76] = 30;  matriz_A[4][77] = 15;   matriz_A[4][78] = 60;  matriz_A[4][79] = 5;
	matriz_A[5][60] = 75;  matriz_A[5][61] = 40;  matriz_A[5][62] = 365; matriz_A[5][63] = 95;  matriz_A[5][64] = 25;  matriz_A[5][65] = 17;   matriz_A[5][66] = 125; matriz_A[5][67] = 20; matriz_A[5][68] = 22;  matriz_A[5][69] = 84;  matriz_A[5][70] = 75;  matriz_A[5][71] = 50;  matriz_A[5][72] = 15;  matriz_A[5][73] = 0;   matriz_A[5][74] = 0;  matriz_A[5][75] = 12; matriz_A[5][76] = 1;   matriz_A[5][77] = 10;   matriz_A[5][78] = 0;   matriz_A[5][79] = 50;
	matriz_A[6][60] = 150; matriz_A[6][61] = 110; matriz_A[6][62] = 70;  matriz_A[6][63] = 10;  matriz_A[6][64] = 20;  matriz_A[6][65] = 30;   matriz_A[6][66] = 104; matriz_A[6][67] = 40; matriz_A[6][68] = 40;  matriz_A[6][69] = 94;  matriz_A[6][70] = 150; matriz_A[6][71] = 50;  matriz_A[6][72] = 10;  matriz_A[6][73] = 50;  matriz_A[6][74] = 50; matriz_A[6][75] = 16; matriz_A[6][76] = 10;  matriz_A[6][77] = 20;   matriz_A[6][78] = 50;  matriz_A[6][79] = 90;
	matriz_A[7][60] = 3;   matriz_A[7][61] = 74;  matriz_A[7][62] = 88;  matriz_A[7][63] = 50;  matriz_A[7][64] = 55;  matriz_A[7][65] = 19;   matriz_A[7][66] = 0;   matriz_A[7][67] = 6;  matriz_A[7][68] = 30;  matriz_A[7][69] = 62;  matriz_A[7][70] = 17;  matriz_A[7][71] = 81;  matriz_A[7][72] = 25;  matriz_A[7][73] = 46;  matriz_A[7][74] = 67; matriz_A[7][75] = 28; matriz_A[7][76] = 36;  matriz_A[7][77] = 8;    matriz_A[7][78] = 1;   matriz_A[7][79] = 52;
	matriz_A[8][60] = 11;  matriz_A[8][61] = 75;  matriz_A[8][62] = 82;  matriz_A[8][63] = 35;  matriz_A[8][64] = 47;  matriz_A[8][65] = 99;   matriz_A[8][66] = 5;   matriz_A[8][67] = 14; matriz_A[8][68] = 23;  matriz_A[8][69] = 38;  matriz_A[8][70] = 48;  matriz_A[8][71] = 14;  matriz_A[8][72] = 5;   matriz_A[8][73] = 72;  matriz_A[8][74] = 14; matriz_A[8][75] = 39; matriz_A[8][76] = 46;  matriz_A[8][77] = 27;   matriz_A[8][78] = 11;  matriz_A[8][79] = 91;
	matriz_A[9][60] = 78;  matriz_A[9][61] = 0;   matriz_A[9][62] = 12;  matriz_A[9][63] = 99;  matriz_A[9][64] = 84;  matriz_A[9][65] = 31;   matriz_A[9][66] = 16;  matriz_A[9][67] = 7;  matriz_A[9][68] = 33;  matriz_A[9][69] = 20;  matriz_A[9][70] = 5;   matriz_A[9][71] = 18;  matriz_A[9][72] = 96;  matriz_A[9][73] = 63;  matriz_A[9][74] = 31; matriz_A[9][75] = 0;  matriz_A[9][76] = 70;  matriz_A[9][77] = 4;    matriz_A[9][78] = 66;  matriz_A[9][79] = 9;
		
	vetor_b[0] = 5000;
	vetor_b[1] = 9000;
	vetor_b[2] = 8000;
	vetor_b[3] = 5500;
	vetor_b[4] = 6000;
	vetor_b[5] = 12000;
	vetor_b[6] = 10000;
	vetor_b[7] = 8000;
	vetor_b[8] = 7000;
	vetor_b[9] = 6500;
	
	// -1 -> menor ou igual. +1 -> maior ou igual. 0 -> igual.
	for(i = 0; i < num_linhas; i++){
        vetor_relacional[i] = -1;
	}

	vetor_custo[0] = 49;
	vetor_custo[1] = 420;
	vetor_custo[2] = 316;
	vetor_custo[3] = 72;
	vetor_custo[4] = 71;
	vetor_custo[5] = 49;
	vetor_custo[6] = 108;
	vetor_custo[7] = 116;
	vetor_custo[8] = 90;
	vetor_custo[9] = 738;
	vetor_custo[10] = 1811;
	vetor_custo[11] = 430;
	vetor_custo[12] = 3060;
	vetor_custo[13] = 215;
	vetor_custo[14] = 58;
	vetor_custo[15] = 296;
	vetor_custo[16] = 620;
	vetor_custo[17] = 418;
	vetor_custo[18] = 47;
	vetor_custo[19] = 81;
	vetor_custo[20] = 460;
	vetor_custo[21] = 200;
	vetor_custo[22] = 1900;
	vetor_custo[23] = 48;
	vetor_custo[24] = 27;
	vetor_custo[25] = 122;
	vetor_custo[26] = 21;
	vetor_custo[27] = 225;
	vetor_custo[28] = 216;
	vetor_custo[29] = 62;
	vetor_custo[30] = 431;
	vetor_custo[31] = 220;
	vetor_custo[32] = 22;
	vetor_custo[33] = 115;
	vetor_custo[34] = 71;
	vetor_custo[35] = 29;
	vetor_custo[36] = 116;
	vetor_custo[37] = 51;
	vetor_custo[38] = 88;
	vetor_custo[39] = 70;
	vetor_custo[40] = 1611;
	vetor_custo[41] = 2860;
	vetor_custo[42] = 38;
	vetor_custo[43] = 420;
	vetor_custo[44] = 27;
	vetor_custo[45] = 1325;
	vetor_custo[46] = 820;
	vetor_custo[47] = 631;
	vetor_custo[48] = 528;
	vetor_custo[49] = 322;
	vetor_custo[50] = 631;
	vetor_custo[51] = 132;
	vetor_custo[52] = 420;
	vetor_custo[53] = 86;
	vetor_custo[54] = 42;
	vetor_custo[55] = 103;
	vetor_custo[56] = 215;
	vetor_custo[57] = 81;
	vetor_custo[58] = 91;
	vetor_custo[59] = 26;
	vetor_custo[60] = 560;
	vetor_custo[61] = 1125;
	vetor_custo[62] = 300;
	vetor_custo[63] = 620;
	vetor_custo[64] = 100;
	vetor_custo[65] = 431;
	vetor_custo[66] = 68;
	vetor_custo[67] = 328;
	vetor_custo[68] = 47;
	vetor_custo[69] = 122;
	vetor_custo[70] = 322;
	vetor_custo[71] = 196;
	vetor_custo[72] = 41;
	vetor_custo[73] = 25;
	vetor_custo[74] = 425;
	vetor_custo[75] = 260;
	vetor_custo[76] = 416;
	vetor_custo[77] = 115;
	vetor_custo[78] = 82;
	vetor_custo[79] = 22;
	
	inicializa_vetor(vetor_x);

	if((arqtxt = fopen("arquivo_branch_bound.txt", "w")) == NULL){
        printf("\nErro abertura arquivo texto - criacao\n\n");
    }
    else{
        t = clock();
        //gauss(matriz_A, vetor_b, vetor_x, num_linhas, num_colunas);
        //simplex(matriz_A, vetor_relacional, vetor_b, minimiza, vetor_custo, vetor_x, &valor_Z, num_linhas, num_colunas);
        total_nos = branch_bound(matriz_A, vetor_relacional, vetor_b, minimiza, vetor_custo, resultados, num_linhas, num_colunas, solucao_otima, arqtxt);
        t = clock() - t;
    }

	//printf("\nO vetor x eh:\n\n");
	//imprime_vetor(vetor_x, num_linhas + num_colunas);


	printf("\nNos gerados: %d\n\n", total_nos);
	fprintf(arqtxt, "\nNos gerados %d:\n\n", total_nos);
	for(i = 0; i < total_nos; i++){
        if((resultados[i])->valor_Z != 0){
            printf("\nNo %d, filho do no %d\n", (resultados[i])->no, (resultados[i])->pai);
            fprintf(arqtxt, "\nNo %d, filho do no %d\n", (resultados[i])->no, (resultados[i])->pai);
            printf("\nO vetor x[%d] eh:\n\n", i);
            imprime_vetor((resultados[i])->vetor_x, num_colunas);
            fprintf(arqtxt, "\nO vetor x[%d] eh:\n\n", i);
            imprime_vetor_arq((resultados[i])->vetor_x, num_colunas, arqtxt);
            if(minimiza == 0){
                printf("\nMax f(x)[%d] = %4.2f\n\n", i, (resultados[i])->valor_Z);
                fprintf(arqtxt, "\nMax f(x)[%d] = %4.2f\n\n", i, (resultados[i])->valor_Z);
            }
            else{
                printf("\nMin f(x)[%d] = %4.2f\n\n", i, (resultados[i])->valor_Z);
                fprintf(arqtxt, "\nMin f(x)[%d] = %4.2f\n\n", i, (resultados[i])->valor_Z);
            }
        }
	}

	printf("\n\n\nO vetor x da solucao otima eh o do no %d:\n\n", solucao_otima->no);
	fprintf(arqtxt, "\n\n\nO vetor x da solucao otima eh o do no %d:\n\n", solucao_otima->no);
	imprime_vetor(solucao_otima->vetor_x, num_colunas);
	imprime_vetor_arq(solucao_otima->vetor_x, num_colunas, arqtxt);

    printf("\nO valor Z da solucao otima, no %d:\n", solucao_otima->no);
	fprintf(arqtxt, "\nO valor Z da solucao otima, no %d:\n", solucao_otima->no);
	if(minimiza == 0){
        printf("\nMax f(x) = %4.2f\n\n", solucao_otima->valor_Z);
        fprintf(arqtxt, "\nMax f(x) = %4.2f\n\n", solucao_otima->valor_Z);
	}
	else{
	    printf("\nMin f(x) = %4.2f\n\n", solucao_otima->valor_Z);
        fprintf(arqtxt, "\nMin f(x) = %4.2f\n\n", solucao_otima->valor_Z);
	}

	printf("\nTempo de execucao: %lf ms\n\n", ((double)t)/(CLOCKS_PER_SEC/1000));
	fprintf(arqtxt, "\nTempo de execucao: %lf ms\n\n", ((double)t)/(CLOCKS_PER_SEC/1000));

    fclose(arqtxt);
	system("pause");
    return 0;
}
