#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
// triplets et leurs 3 paires
typedef struct {
	int8_t coverage ;
	int8_t remaining ;
} pair ;


int8_t tp[56][3];
pair pairs[28] ;

// int8_t coverage[28];
// int8_t remaining[28]; // nb de triplets >= current ti pouvant encore choisir cette paire
long long count = 0;

int8_t bestUp = 0 ;
// Pour chaque paire, combien de triplets restants (>= ti) la contiennent
// On va maintenir un compteur "remaining_uses[pi]" = nb de triplets >= ti qui ont pi comme option
// Si coverage[pi]==1 et remaining_uses[pi]==0 -> impossible d'atteindre 2 -> coupe
// Si coverage[pi]==0 et remaining_uses[pi]<2 -> impossible -> coupe

// Pré-calculer pour chaque paire, la liste des triplets qui l'utilisent (triée)
// On utilise une approche: pour chaque triplet ti, on décrémente les remaining au début et incrémente au retour


clock_t begin ;

/* here, do your time-consuming job */



void init() {
    int8_t pi = 0;
    int8_t pair_idx[8][8];
    memset(pair_idx, -1, sizeof(pair_idx));
    for (int8_t i = 0; i < 8; i++)
        for (int8_t j = i+1; j < 8; j++)
            pair_idx[i][j] = pi++;

    int8_t ti = 0;
    for (int8_t a = 0; a < 8; a++)
      for (int8_t b = a+1; b < 8; b++)
        for (int8_t c = b+1; c < 8; c++) {
            tp[ti][0] = pair_idx[a][b];
            tp[ti][1] = pair_idx[a][c];
            tp[ti][2] = pair_idx[b][c];
            ti++;
        }
    // remaining[pi] = nombre total de triplets qui contiennent la paire pi
    // = 6 pour chaque paire dans C(8,3)
    for (int8_t p = 0; p < 28; p++) { pairs[p].remaining = 6; pairs[p].coverage = 0 ; }
}


void bt2(int8_t ti) {
static long long last = 0 ;
	if(ti < bestUp) {
		bestUp = ti ;
		clock_t end = clock();
		printf("BestUp=%d count=%lld in %.1fs\n",bestUp,count,(float)(end - begin) / CLOCKS_PER_SEC);
		fflush(stdout);
	}
    if (ti == 56) {
		if (count == 0) bestUp = 56 ;
		count++;
		last++ ;
		if (last >= 10000000000LL) {
			last = 0 ;
			printf("count=%lld \n",count) ;
			fflush(stdout);
		}
		return;
	}

    // Retirer ce triplet du "remaining" de ses 3 paires
    int8_t p0 = tp[ti][0], p1 = tp[ti][1], p2 = tp[ti][2];

	pairs[p0].remaining -- ;
    pairs[p1].remaining -- ;
    pairs[p2].remaining -- ;
	int8_t p0_1 = tp[ti+1][0], p1_1 = tp[ti+1][1], p2_1 = tp[ti+1][2];
	// test p0_1 p1_1 et p2_1 sont tous differents de p0,p1,p2
	if ( ((1<<p0) | (1<<p1) | (1<<p2)) & ((1<<p0_1) | (1<<p1_1) | (1<<p2_1)) ==0  ) {
		pairs[p0_1].remaining -- ;
		pairs[p1_1].remaining -- ;
		pairs[p2_1].remaining -- ;
		int Test = 0 ;
		if ((pairs[p0].coverage < 2) &&  (pairs[p1].remaining+pairs[p1].coverage >= 2) && (pairs[p2].remaining+pairs[p2].coverage >= 2)) {
			pairs[p0].coverage++ ;
				if ((pairs[p0_1].coverage < 2) && (pairs[p1_1].remaining + pairs[p1_1].coverage >=2) && (pairs[p2_1].remaining + pairs[p2_1].coverage >=2)) {
					pairs[p0_1].coverage++ ; bt2(ti + 2); pairs[p0_1].coverage-- ;
					Test |= 1 ;
				}
				if ((pairs[p1_1].coverage < 2) && (pairs[p0_1].remaining + pairs[p0_1].coverage >=2) && (pairs[p2_1].remaining + pairs[p2_1].coverage >=2)) {
					pairs[p1_1].coverage++ ;bt2(ti + 2); pairs[p1_1].coverage-- ;
					Test |= 2 ;
				}
				if ((pairs[p2_1].coverage < 2) && (pairs[p0_1].remaining + pairs[p0_1].coverage >=2) && (pairs[p1_1].remaining + pairs[p1_1].coverage >=2)) {
					pairs[p2_1].coverage++ ; bt2(ti + 2); 	pairs[p2_1].coverage-- ;
					Test |= 4 ;
				}
			if( ti==0) {
				bestUp = ti ;
				clock_t end = clock();
				printf("BestUp=%d* count=%lld in %.1fs\n",bestUp,count,(float)(end - begin) / CLOCKS_PER_SEC);
				fflush(stdout);
			}
			pairs[p0].coverage-- ; ;
		}
		if ((pairs[p1].coverage < 2) &&  (pairs[p0].remaining+pairs[p0].coverage >= 2) && (pairs[p2].remaining+pairs[p2].coverage >= 2)) {
			pairs[p1].coverage++ ;
				if(Test & 1) { pairs[p0_1].coverage++ ; bt2(ti + 2);  pairs[p0_1].coverage-- ;	}
				if(Test & 2) { pairs[p1_1].coverage++ ; bt2(ti + 2);  pairs[p1_1].coverage-- ;	}
				if(Test & 4) { pairs[p2_1].coverage++ ;	bt2(ti + 2);  pairs[p2_1].coverage-- ;	}
			pairs[p1].coverage-- ;
		}
		if ((pairs[p2].coverage < 2) &&  (pairs[p0].remaining+pairs[p0].coverage >= 2) && (pairs[p1].remaining+pairs[p1].coverage >= 2)) {
			pairs[p2].coverage++ ;
				if(Test & 1) { pairs[p0_1].coverage++ ; bt2(ti + 2); pairs[p0_1].coverage-- ; }
				if(Test & 2) { pairs[p1_1].coverage++ ;	bt2(ti + 2); pairs[p1_1].coverage-- ; }
				if(Test & 4) { pairs[p2_1].coverage++ ;	bt2(ti + 2); pairs[p2_1].coverage-- ; }
			pairs[p2].coverage-- ;
		}
		pairs[p0_1].remaining++;
		pairs[p1_1].remaining++;
		pairs[p2_1].remaining++;
	} else {
		if ((pairs[p0].coverage < 2) &&  (pairs[p1].remaining+pairs[p1].coverage >= 2) && (pairs[p2].remaining+pairs[p2].coverage >= 2)) {
			pairs[p0].coverage++ ;
				pairs[p0_1].remaining -- ;
				pairs[p1_1].remaining -- ;
				pairs[p2_1].remaining -- ;
				if ((pairs[p0_1].coverage < 2) && (pairs[p1_1].remaining + pairs[p1_1].coverage >=2) && (pairs[p2_1].remaining + pairs[p2_1].coverage >=2)) {
					pairs[p0_1].coverage++ ; bt2(ti + 2); pairs[p0_1].coverage-- ;
				}
				if ((pairs[p1_1].coverage < 2) && (pairs[p0_1].remaining + pairs[p0_1].coverage >=2) && (pairs[p2_1].remaining + pairs[p2_1].coverage >=2)) {
					pairs[p1_1].coverage++ ; bt2(ti + 2); pairs[p1_1].coverage-- ;
				}
				if ((pairs[p2_1].coverage < 2) && (pairs[p0_1].remaining + pairs[p0_1].coverage >=2) && (pairs[p1_1].remaining + pairs[p1_1].coverage >=2)) {
					pairs[p2_1].coverage++ ; bt2(ti + 2); pairs[p2_1].coverage-- ;
				}
				pairs[p0_1].remaining ++ ;
				pairs[p1_1].remaining ++ ;
				pairs[p2_1].remaining ++ ;
			if( ti==0) {
				bestUp = ti ;
				clock_t end = clock();
				printf("BestUp=%d* count=%lld in %.1fs\n",bestUp,count,(float)(end - begin) / CLOCKS_PER_SEC);
				fflush(stdout);
			}

			pairs[p0].coverage-- ; ;
		}
		if ((pairs[p1].coverage < 2) &&  (pairs[p0].remaining+pairs[p0].coverage >= 2) && (pairs[p2].remaining+pairs[p2].coverage >= 2)) {
			pairs[p1].coverage++ ;
				pairs[p0_1].remaining--;
				pairs[p1_1].remaining--;
				pairs[p2_1].remaining--;
				if ((pairs[p0_1].coverage < 2) && (pairs[p1_1].remaining + pairs[p1_1].coverage >=2) && (pairs[p2_1].remaining + pairs[p2_1].coverage >=2)) {
					pairs[p0_1].coverage++ ; bt2(ti + 2); pairs[p0_1].coverage-- ;
				}
				if ((pairs[p1_1].coverage < 2) && (pairs[p0_1].remaining + pairs[p0_1].coverage >=2) && (pairs[p2_1].remaining + pairs[p2_1].coverage >=2)) {
					pairs[p1_1].coverage++ ; bt2(ti + 2); pairs[p1_1].coverage-- ;
				}
				if ((pairs[p2_1].coverage < 2) && (pairs[p0_1].remaining + pairs[p0_1].coverage >=2) && (pairs[p1_1].remaining + pairs[p1_1].coverage >=2)) {
					pairs[p2_1].coverage++ ; bt2(ti + 2); pairs[p2_1].coverage-- ;
				}
				pairs[p0_1].remaining++;
				pairs[p1_1].remaining++;
				pairs[p2_1].remaining ++;
			pairs[p1].coverage-- ;
		}
		if ((pairs[p2].coverage < 2) &&  (pairs[p0].remaining+pairs[p0].coverage >= 2) && (pairs[p1].remaining+pairs[p1].coverage >= 2)) {
			pairs[p2].coverage++ ;
				pairs[p0_1].remaining--;
				pairs[p1_1].remaining--;
				pairs[p2_1].remaining--;
				if ((pairs[p0_1].coverage < 2) && (pairs[p1_1].remaining + pairs[p1_1].coverage >=2) && (pairs[p2_1].remaining + pairs[p2_1].coverage >=2)) {
					pairs[p0_1].coverage++ ; bt2(ti + 2); pairs[p0_1].coverage-- ;
				}
				if ((pairs[p1_1].coverage < 2) && (pairs[p0_1].remaining + pairs[p0_1].coverage >=2) && (pairs[p2_1].remaining + pairs[p2_1].coverage >=2)) {
					pairs[p1_1].coverage++ ; bt2(ti + 2); pairs[p1_1].coverage-- ;
				}
				if ((pairs[p2_1].coverage < 2) && (pairs[p0_1].remaining + pairs[p0_1].coverage >=2) && (pairs[p1_1].remaining + pairs[p1_1].coverage >=2)) {
					pairs[p2_1].coverage++ ; bt2(ti + 2); pairs[p2_1].coverage-- ;
				}
				pairs[p0_1].remaining++;
				pairs[p1_1].remaining++;
				pairs[p2_1].remaining++;
			pairs[p2].coverage-- ;
		}
	}
	pairs[p0].remaining++;
    pairs[p1].remaining++;
	pairs[p2].remaining++;
}

int main() {
    init();
	begin = clock() ;
    bt2(0);
    printf("Solutions: %lld\n", count);
    return 0;
}


