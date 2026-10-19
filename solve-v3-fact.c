#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// sigma permitation 4 paires (0,1) (2,3) (4,5) (6,7) soit {0,1,3,2,5,4,7,6} 
int8_t sigma_t3[56] ; // correspondacne triplet <-> triplet (sigma)
int8_t sigma_p2[28] ;

int8_t sigma_1_t3[56] ; 
int8_t sigma_1_p2[28] ;
int8_t tp[56][3]; // 3 paires possibles par triplet 
int8_t t3[56][3]; // valeur des triplets
int8_t p2[28][2]; // valeur des paires 

// triplets et leurs 3 paires
int8_t coverage[28];
long long count = 0;
long long count_sigma = 0 ;

int8_t histot2p[56] ; // correspondance triplet to paire retenu

int8_t bestUp = 0 ; 
// Pour chaque paire, combien de triplets restants (>= ti) la contiennent
// On va maintenir un compteur "remaining_uses[pi]" = nb de triplets >= ti qui ont pi comme option
// Si coverage[pi]==1 et remaining_uses[pi]==0 -> impossible d'atteindre 2 -> coupe
// Si coverage[pi]==0 et remaining_uses[pi]<2 -> impossible -> coupe

// Pré-calculer pour chaque paire, la liste des triplets qui l'utilisent (triée)
// On utilise une approche: pour chaque triplet ti, on décrémente les remaining au début et incrémente au retour

int8_t remaining[28]; // nb de triplets >= current ti pouvant encore choisir cette paire

clock_t begin ;

/* here, do your time-consuming job */



void init() {
    int8_t pi = 0;
    int8_t pair_idx[8][8];
    memset(pair_idx, -1, sizeof(pair_idx));
    for (int8_t i = 0; i < 8; i++)
        for (int8_t j = i+1; j < 8; j++) {
			p2[pi][0] = i ; p2[pi][1] = j ;
            pair_idx[i][j] = pi++;
		}

    int8_t ti = 0;
    for (int8_t a = 0; a < 8; a++)
      for (int8_t b = a+1; b < 8; b++)
        for (int8_t c = b+1; c < 8; c++) {
			t3[ti][0] = a ; t3[ti][1] = b ; t3[ti][2] = c ;
            tp[ti][0] = pair_idx[a][b];
            tp[ti][1] = pair_idx[a][c];
            tp[ti][2] = pair_idx[b][c];
            ti++;
        }
    // remaining[pi] = nombre total de triplets qui contiennent la paire pi
    // = 6 pour chaque paire dans C(8,3)
    for (int8_t p = 0; p < 28; p++) remaining[p] = 6;
}

int best_it3 = 0 ;
void bt(int8_t ti) {
static long long last = 0 ;
	// printf("%d ",ti);
	if(ti < bestUp) {
		bestUp = ti ;
		clock_t end = clock();
		printf("BestUp=%d count=%lld %lld in %.1fs\n",bestUp,count,count_sigma,(float)(end - begin) / CLOCKS_PER_SEC);
		fflush(stdout);
	}
    if (ti == 56) {
		if (count == 0) bestUp = 56 ;
		int is_sigma = 1 ;		
		if (is_sigma) 
		{
			/*printf("---------------------- sol n° %d ---------------------\n",count_sigma);
			for(int8_t it3=0; it3 < 56 ; it3 ++ ) {
				printf("{%d,%d,%d}->(%d,%d) = sigma => {%d,%d,%d}->(%d,%d) sigma(%d,%d) = (%d,%d) \n"
				,t3[it3][0],t3[it3][1],t3[it3][2],p2[histot2p[it3]][0],p2[histot2p[it3]][1]
				,t3[sigma_t3[it3]][0],t3[sigma_t3[it3]][1],t3[sigma_t3[it3]][2]
				,p2[sigma_p2[histot2p[it3]]][0],p2[sigma_p2[histot2p[it3]]][1]
				,p2[histot2p[it3]][0],p2[histot2p[it3]][1]
				,p2[sigma_p2[histot2p[it3]]][0],p2[sigma_p2[histot2p[it3]]][1]
				);
			}
			printf("--------------------End sol n° %d ---------------------\n",count_sigma); */
			count_sigma++ ;
		}
		count++; 
		if (count <= 4) {
			printf("---------------------- sol n° %d  is_sigma= %d------------------\n",count,is_sigma);
			for(int8_t it3=0; it3 < 56 ; it3 ++ ) {
				int8_t it3s = it3 ;
				int8_t ip2 = histot2p[it3] ;
				int8_t ip2s = ip2 ; 
	
				printf("S{%d,%d,%d}=(%d,%d) => S(Sgm(T)={%d,%d,%d})=(%d,%d) Sgm((%d,%d))=(%d,%d) ?? S{%d,%d,%d}=(%d,%d) %c \n"
				,t3[it3][0],t3[it3][1],t3[it3][2],p2[ip2][0],p2[ip2][1]
				,t3[sigma_t3[it3]][0],t3[sigma_t3[it3]][1],t3[sigma_t3[it3]][2]
				,p2[ip2s][0],p2[ip2s][1]
				,p2[ip2][0],p2[ip2][1],p2[sigma_p2[ip2]][0],p2[sigma_p2[ip2]][1]
				,t3[it3s][0],t3[it3s][1],t3[it3s][2],p2[ip2s][0],p2[ip2s][1]
				,(sigma_p2[ip2] ==ip2s) ? '*' : ' '
				);
			}
			printf("--------------------End sol n° %d ---------------------\n",count);

		}	
		last++ ;
		if (last >= 10000000000LL) {
			last = 0 ;
			printf("count=%lld count_sigma=%lld\n",count,count_sigma) ;
			fflush(stdout);
		}
		return; 
	}

    // Retirer ce triplet du "remaining" de ses 3 paires
    int8_t p0 = tp[ti][0], p1 = tp[ti][1], p2 = tp[ti][2];
    remaining[p0]--;
    remaining[p1]--;
    remaining[p2]--;
	if ((coverage[p0] < 2) &&  (remaining[p1]+coverage[p1] >= 2) && (remaining[p2]+coverage[p2] >= 2)) {
		coverage[p0]++ ;
		histot2p[ti] = p0 ;
		bt(ti + 1);
		coverage[p0]-- ;
	}
	/*
	remaining[p0]++;
	remaining[p1]++;
    remaining[p2]++;

    remaining[p0]--;
	remaining[p1]--;
    remaining[p2]--;
	*/
	if ((coverage[p1] < 2) &&  (remaining[p0]+coverage[p0] >= 2) && (remaining[p2]+coverage[p2] >= 2)) {
		coverage[p1]++ ;
//		histot2p[ti] = p1 ;
//		bt(ti + 1);
		histot2p[ti] = p1 ;
		bt(ti + 1);
		coverage[p1]-- ;
	}
	/*
	remaining[p0]++;
	remaining[p1]++;
    remaining[p2]++;

    remaining[p0]--;
    remaining[p1]--;
	remaining[p2]--;
	*/
	if ((coverage[p2] < 2) &&  (remaining[p0]+coverage[p0] >= 2) && (remaining[p1]+coverage[p1] >= 2)) {
		coverage[p2]++ ;
		histot2p[ti] = p2 ;
		bt(ti + 1);
		coverage[p2]-- ;
	}
	remaining[p0]++;
    remaining[p1]++;
	remaining[p2]++;
}

int main() {
    init();
	begin = clock() ;
    memset(coverage, 0, sizeof(coverage));
    bt(0);
    printf("Solutions: %lld %lld\n", count,count_sigma);
    return 0;
}


