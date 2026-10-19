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
		/*
Non, le type (2,2,2,2) n'est pas le seul avec Fix > 0. Il y en a 6 en tout :
TypeOrdre 			σ		Fix(σ)	Contribution
(1,1,1,1,1,1,1,1)	1		N		N
(2,2,2,1,1)			2		187120	78590400 = 420 x 187120 ; 420 = 28 x 5 x 3 
(2,2,2,2)			2		 	166750080 105 x 1588096 ; 105 = 7x5x3 
(4,4)				4		2960	3729600 = 1260 x 2970 ; 1260 = 35 x 3! x 3!
(7,1)				7		432		2488320 = 5760 x 432 ; 5760 = 8 x 6! 
(8,)				8		112		564480 = 5040 x 112 ; 5040 = 7! 
		*/
				
//	int8_t compl[8] = {1,2,3,4,5,6,7,0} ; // (8) = 0  150 150			
	int8_t compl[8] = {1,2,3,4,5,6,0,7} ; // (7,1) =  432 432			
//	int8_t compl[8] = {1,2,3,4,5,0,7,6} ; // (6,2) = 0				
//	int8_t compl[8] = {1,2,3,4,5,0,6,7} ; // (6,1,1) = 0				
//	int8_t compl[8] = {1,2,3,4,0,6,7,5} ; // (5,3) = 0				
//	int8_t compl[8] = {1,2,3,4,0,6,5,7} ; // (5,2,1) = 0		
//	int8_t compl[8] = {1,2,3,4,0,5,6,7} ; // (5,1,1,1) = 0		
//	int8_t compl[8] = {1,2,3,0,5,6,7,4} ; // (4,4) = 6444 6444
//	int8_t compl[8] = {1,2,3,0,5,6,4,7} ; // (4,3,1) = 0
//	int8_t compl[8] = {1,0,3,2,5,6,7,4} ; // (4,2,2) = 0
//	int8_t compl[8] = {1,2,3,0,5,4,6,7} ; // (4,2,1,1) = 0
//	int8_t compl[8] = {1,2,3,0,4,5,6,7} ; // (4,1,1,1,1) = 0
//	int8_t compl[8] = {1,2,0,4,5,3,7,6} ; // (3,3,2) = 0	
//	int8_t compl[8] = {1,2,0,4,5,3,6,7} ; // (3,3,1,1) = 0	
//	int8_t compl[8] = {1,2,0,4,3,6,5,7} ; // (3,2,2,1) = 0	
//	int8_t compl[8] = {1,2,0,4,3,5,6,7} ; // (3,2,1,1,1) = 0
//	int8_t compl[8] = {1,2,0,3,4,5,6,7} ; // (3,1,1,1,1,1) = 0
//	int8_t compl[8] = {1,0,3,2,5,4,7,6} ; // (2,2,2,2) = 8468664 8468664	
//	int8_t compl[8] = {1,0,3,2,5,4,6,7} ; // (2,2,2,1,1) = 280680 280680
//	int8_t compl[8] = {2,3,0,1,5,4,7,6} ; // (2,2,2,2) = 8468664
//	int8_t compl[8] = {1,0,3,2,4,5,6,7} ; // (2,2,1,1,1,1) = 0
//	int8_t compl[8] = {1,0,2,3,4,5,6,7} ; // (2,1,1,1,1,1,1) = 0
	for(int k =0; k<8;k++) printf("%d ",compl[k]) ; printf("\n");
	for(int8_t nt3 = 0 ; nt3 < 56 ;nt3++ ){
		int8_t v0 = compl[t3[nt3][0]] ;
		int8_t v1 = compl[t3[nt3][1]] ;
		int8_t v2 = compl[t3[nt3][2]] ;
		if(v0 > v1) { int8_t vtmp = v0 ; v0 = v1 ; v1= vtmp ; }
		if(v1 > v2) { 
			int8_t vtmp = v1 ; v1 = v2 ; v2= vtmp ;
			if(v0 > v1) {int8_t vtmp = v0 ; v0 = v1 ; v1= vtmp ;  }
		}
		for (int8_t nt3_c = 0 ; nt3_c < 56; nt3_c++) {
			if(v0 == t3[nt3_c][0] && v1 == t3[nt3_c][1] && v2 == t3[nt3_c][2]) {
				sigma_t3[nt3] = nt3_c ; sigma_1_t3[nt3_c] = nt3 ; break ;
			}
		}
	}
	for(int8_t np2 = 0 ; np2 < 28 ;np2++ ){
		int8_t v0 = compl[p2[np2][0]] ;
		int8_t v1 = compl[p2[np2][1]] ;
		if(v0 > v1) { int8_t vtmp = v0 ; v0 = v1 ; v1= vtmp ; }
		for (int8_t np2_c = 0 ; np2_c < 28; np2_c++) {
			if(v0 == p2[np2_c][0] && v1 == p2[np2_c][1]) {
				sigma_p2[np2] = np2_c ; sigma_1_p2[np2_c] = np2 ; break ;
			}
		}
	}
	
	for(int8_t nt3 = 0 ; nt3 < 56 ;nt3++) {
		int8_t nt3_c = sigma_t3[nt3] ;
		int8_t nt3_c_1 = sigma_1_t3[nt3_c] ;
		printf("{%d,%d,%d}->{%d,%d,%d}->{%d,%d,%d}\n",t3[nt3][0],t3[nt3][1],t3[nt3][2]
		,t3[nt3_c][0],t3[nt3_c][1],t3[nt3_c][2]
		,t3[nt3_c_1][0],t3[nt3_c_1][1],t3[nt3_c_1][2]);
	}
	for(int8_t np2 = 0 ; np2 < 28 ;np2++) {
		int8_t np2_c = sigma_p2[np2] ;
		int8_t np2_c_1 = sigma_1_p2[np2_c] ;
		printf("{%d,%d}->{%d,%d}->{%d,%d}\n",p2[np2][0],p2[np2][1]
		,p2[np2_c][0],p2[np2_c][1]  ,p2[np2_c_1][0],p2[np2_c_1][1]);
	}
	
	// exit(0) ; 
	
	
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
		for(int8_t it3=0; it3 < 56 ; it3 ++ ) {
			int8_t it3s = sigma_t3[it3] ;
			int8_t ip2 = histot2p[it3] ;
			int8_t ip2s = histot2p[it3s] ; 
		//	if(sigma_p2[histot2p[it3s]] != ip2) {
			if(ip2s != sigma_p2[ip2] ) {
//			if(sigma_p2[histot2p[it3s]] != ip2) {			
				is_sigma= 0 ; 
				 if(it3 > best_it3) {
					printf("Best it3 %d %lld %lld \n",it3,count,count_sigma);
					best_it3= it3 ;
				}
				
			if(count <=3) printf("T(it3=%d)={%d,%d,%d} sigma(it3=%d)=%d=(%d,%d,%d) S(T)=%d=(%d,%d) S(sigma(T))=%d=(%d,%d)\n"
				,it3,t3[it3][0],t3[it3][1],t3[it3][2]
				,it3,it3s,t3[it3s][0],t3[it3s][1],t3[it3s][2]
				,ip2,p2[ip2][0],p2[ip2][1]
				,ip2s,p2[ip2s][0],p2[ip2s][1]);
				
				break ;
			}
		}
		
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
				int8_t it3s = sigma_t3[it3] ;
				int8_t ip2 = histot2p[it3] ;
				int8_t ip2s = histot2p[it3s] ; 
	
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
		if ( sigma_1_t3[ti] < ti ) {	
			int8_t ti_s = sigma_1_t3[ti] ;
			int8_t ip2s = histot2p[ti_s] ; 
			if(sigma_p2[ip2s] == p0) {
				//histot2p[ti] = sigma_1_p2[p0] ;
				histot2p[ti] = p0 ;
				bt(ti + 1);
			}
		} else if (sigma_t3[ti] < ti ) {
			int8_t ti_s = sigma_t3[ti] ;
			int8_t ip2s = histot2p[ti_s] ; 
			if (sigma_1_p2[ip2s] == p0) {
				histot2p[ti] = p0 ;
				bt(ti + 1);
			}
		} else if (sigma_t3[ti] == ti ) {
			//int8_t ip2 = histot2p[ti] ; 
			if (sigma_p2[p0] == p0) 
			{
				histot2p[ti] = p0 ;
				bt(ti + 1);
			}
		} else {
			histot2p[ti] = p0 ;
			bt(ti + 1);
		}
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
		if ( sigma_1_t3[ti] < ti ) {	
			int8_t ti_s = sigma_1_t3[ti] ;
			int8_t ip2s = histot2p[ti_s] ; 
			if(sigma_p2[ip2s] == p1) {
				//histot2p[ti] = sigma_1_p2[p1] ;
				histot2p[ti] = p1 ;
				bt(ti + 1);
			}
		} else if (sigma_t3[ti] < ti ) {
			int8_t ti_s = sigma_t3[ti] ;
			int8_t ip2s = histot2p[ti_s] ; 
			if (sigma_1_p2[ip2s] == p1) {
				histot2p[ti] = p1 ;
				bt(ti + 1);
			}
		} else if (sigma_t3[ti] == ti ) {
			// int8_t ip2 = histot2p[ti] ; 
			// if (ip2 == p1) 
			if (sigma_p2[p1] == p1) 
			{
				histot2p[ti] = p1 ;
				bt(ti + 1);
			}

		} else {
			histot2p[ti] = p1 ;
			bt(ti + 1);
		}
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
//		histot2p[ti] = p2 ;
//		 bt(ti + 1);
		if ( sigma_1_t3[ti] < ti ) {	
			int8_t ti_s = sigma_1_t3[ti] ;
			int8_t ip2s = histot2p[ti_s] ; 
			if(sigma_p2[ip2s] == p2) {
				//histot2p[ti] = sigma_1_p2[p2] ;
				histot2p[ti] = p2 ;
				bt(ti + 1);
			}
		} else if (sigma_t3[ti] < ti ) {
			int8_t ti_s = sigma_t3[ti] ;
			int8_t ip2s = histot2p[ti_s] ; 
			if (sigma_1_p2[ip2s] == p2) {
				histot2p[ti] = p2 ;
				bt(ti + 1);
			}
		} else if (sigma_t3[ti] == ti ) {
			//int8_t ip2 = histot2p[ti] ; 
			//if (ip2 == p2) 
			if (sigma_p2[p2] == p2) 
			{
				histot2p[ti] = p2 ;
				bt(ti + 1);
			}

		} else {
			histot2p[ti] = p2 ;
			bt(ti + 1);
		}

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


