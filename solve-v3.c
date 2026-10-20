#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// triplets et leurs 3 paires
int8_t tp[56][3];
int8_t coverage[28];
long long count = 0;

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
    for (int8_t p = 0; p < 28; p++) remaining[p] = 6;
}

// Vérifie si les paires non encore satisfaites peuvent encore l'être
// (élagage anticipé)
static inline int feasible() {
    for (int8_t p = 0; p < 28; p++) {
		if (remaining[p]+coverage[p] < 2) return 0 ;
    }
    return 1;
}

void bt(int8_t ti) {
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
    remaining[p0]--;
    remaining[p1]--;
    remaining[p2]--;
	if ((coverage[p0] < 2) &&  (remaining[p1]+coverage[p1] >= 2) && (remaining[p2]+coverage[p2] >= 2)) {
		coverage[p0]++ ;
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
    printf("Solutions: %lld\n", count);
    return 0;
}


