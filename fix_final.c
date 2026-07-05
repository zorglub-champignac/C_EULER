//
// Created by Jeannot on 05/07/2026.
//

// /tmp/fix_final.c
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
/* ============================================================
   Calcul de Fix(sigma) pour tous les types de S8.

   Principe : on itère sur les orbites de triplets sous sigma.
   Pour chaque orbite {t0, sigma(t0), ..., sigma^{k-1}(t0)},
   on choisit f(t0) = p, ce qui force f(sigma^j(t0)) = sigma^j(p).
   Chaque paire sigma^j(p) reçoit +1 dans cov[].
   Contrainte finale : cov[p] == 2 pour toute paire p.
   ============================================================ */

int sg[8];

/* Paires et triplets */
int P[28][2]; // liste des paires
int pidx[8][8]; // numero ds la liste d'une paire
int Tr[56][3]; // liste des triplets
int tidx[8][8][8]; // numero ds la liste d'un triplet

// #define Ident
// nombre de decomposition en cycles
#if defined Ident
	#define NB_PART	22
#else
	#define NB_PART 21
#endif
/* Orbites de triplets */
#define MAXORB 56
#define MAXK   15
int norbs;
int oksz[MAXORB];           /* taille de l'orbite o              */
int onc[MAXORB];            /* nb de choix valides pour l'orbite */
int ochoice[MAXORB][3][MAXK]; /* ochoice[o][m][j] = paire à l'étape j pour le choix m */

int cov[28];
long long cnt;

/* --- Fonctions de base --- */
// remet dans l'ordre alpha un triplet issu de l'application de la permutation sigma (sg)
void appt(int in[3], int out[3]){
    int a=sg[in[0]], b=sg[in[1]], c=sg[in[2]];
    if(a>b){int t=a;a=b;b=t;} if(b>c){int t=b;b=c;c=t;} if(a>b){int t=a;a=b;b=t;}
    out[0]=a; out[1]=b; out[2]=c;
}
// remet dans l'ordre alpha une paire issue de l'application de la permutation sigma (sg)
// et retourne le numero de la paire (parmi 28)
int appp(int pi){
    int a=sg[P[pi][0]], b=sg[P[pi][1]];
    if(a>b){int t=a;a=b;b=t;} return pidx[a][b];
}
// retourne True si pi, la paire d'indice i, fait partie du triplet t[3]
int pin(int pi, int t[3]){
    int a=P[pi][0], b=P[pi][1];
    return (a==t[0]&&b==t[1])||(a==t[0]&&b==t[2])||(a==t[1]&&b==t[2]);
}
// #define DEBUG

/* --- Construction des orbites et de leurs choix pour chaque sigma--- */
void build(){
    int np=0;
    memset(pidx, 0, sizeof(pidx));
    for(int i=0;i<8;i++) {
		for(int j=i+1;j<8;j++){
			P[np][0]=i; P[np][1]=j; pidx[i][j]=np++;
		}
	}
    int nt=0;
    for(int a=0;a<8;a++) {
		for(int b=a+1;b<8;b++) {
			for(int c=b+1;c<8;c++){
				Tr[nt][0]=a; Tr[nt][1]=b; Tr[nt][2]=c; tidx[a][b][c]=nt++;
			}
		}
	}
    int vis[56]={0};
    norbs=0;
	// int orb[56][MAXK][3];
	for(int it = 0; it < 56; it ++) {
		if(vis[it]) continue ;
        /* Construire l'orbite du triplet c-a-d la liste ses triplets sigma^k (Tr=cur)*/
        int orb[MAXK][3];
		int k=0 ;
		int cur[3]={Tr[it][0],Tr[it][1],Tr[it][2]};
		// calcul du cycle par zima du tripet n° i
        while(!vis[tidx[cur[0]][cur[1]][cur[2]]]){
            int idx=tidx[cur[0]][cur[1]][cur[2]];
            vis[idx]=1;
            memcpy(orb[k], cur,sizeof(cur));
            k++;
            int nxt[3]; appt(cur,nxt); memcpy(cur,nxt,sizeof(cur));
        }
#if defined DEBUG
		printf("it=%d ",it);
		for(int ik=0;ik<k;ik++) {
			printf("(%d,%d,%d)",orb[ik][0],orb[ik][1],orb[ik][2]);
		}
		printf("\n");
#endif
		// on stocke la longueur du cycle du triplet n° it
        oksz[norbs]=k;
        /* Choix valides : paires p de orb[0] tq sigma^j(p) dans orb[j] */
        int t0[3]; memcpy(t0, orb[0], 12);
        int ps[3][2]={{t0[0],t0[1]},{t0[0],t0[2]},{t0[1],t0[2]}};
        int nc=0;
        for(int m=0;m<3;m++){ // boucle sur les 3 paires associees au triplet n° it
            int pi=pidx[ps[m][0]][ps[m][1]]; // n° de la paire
            int cp=pi, valid=1, tmp[MAXK];
			tmp[0] = cp ; cp = appp(cp); // cas j=0, on sait que la paire est ds le triplet.
            for(int j=1;j<k;j++){ // on teste que sigma^j[ pi ] est dans orb[j]
                if(!pin(cp,orb[j])){
					valid=0; break;
				}
                tmp[j]=cp;
                cp=appp(cp);
            }
            if(valid && cp==pi){ // paire valide qui boucle
                for(int j=0;j<k;j++) {
					ochoice[norbs][nc][j]=tmp[j];
#if defined DEBUG
#endif
				}
                nc++;
            }
        }
#if defined DEBUG
		printf("it=%d => onc[%d]=%d \n",it,norbs,nc); fflush(stdout);
#endif
        onc[norbs]=nc; // nombre de paires parmi 3 qui bouclent (si zéro triplet non valide pour sigma
        norbs++;
    }
#if defined DEBUG
	printf("Dump ochoice norbs=%d\n",norbs);
	for (int inorbs=0;inorbs<norbs;inorbs++) {
		printf("inorbs=%d ",inorbs);
		for (int inc=0;inc<onc[inorbs];inc++) {
			printf("nc=%d {",inc);
			for (int i=0;i<oksz[inorbs];i++) {
				printf("%c%d",i==0 ?' ':',',ochoice[inorbs][inc][i]);
			}
			printf("}");
		}
		printf("\n");
	}
	printf("\n");
#endif
}

/* --- Backtracking SANS élagage (version sûre) --- */
// developpement recursif
#if defined DEBUG
int t2p[56] ;
#endif

void bt(int idx){
#if defined DEBUG
	printf("\nbt_%d norbs=%d",idx,norbs);
#endif
    if(idx==norbs){ // on a developpe toutes les orbites correspondant a sigma
        for(int p=0;p<28;p++) if(cov[p]!=2) return;
#if defined DEBUG
		{
			for(int i=0;i<56;i++) {
				printf("{%d,%d,%d}-[%d](%d,%d) ",Tr[i][0],Tr[i][1],Tr[i][2],t2p[i],P[t2p[i]][0],P[t2p[i]][1]);
			}
			printf("\n");
		}
		exit(0) ;
#endif
        cnt++;
        return;
    }
#if defined DEBUG
	printf("oksz[%d]=%d ",idx,oksz[idx]);
#endif
    int k=oksz[idx]; // longueur de l'orbite
    for(int m=0;m<onc[idx];m++){
        /* Calcul du delta pour ce choix */
        int d[28]={0};
        for(int j=0;j<k;j++) {
			d[ochoice[idx][m][j]]++ ;
		}
        /* Vérification : ne pas dépasser 2 */
        int ok=1;
        for(int p=0;p<28 && ok ;p++) {
			if(cov[p]+d[p]>2) ok=0;
		}
        if(!ok) continue;
        /* Appliquer */

        for(int p=0;p<28;p++) {
			cov[p]+=d[p];
#if defined(DEBUG)
			if(d[p]) {
				t2p[idx] = p;
				if(idx == 0) {
					printf("t2p[%d]=%d ",idx,p); fflush(stdout);
				}
			}
#endif
		}
        bt(idx+1);
#if defined(DEBUG)
		printf("rewind");
#endif
        for(int p=0;p<28;p++) cov[p]-=d[p];
    }
}

long long fix(int s[8]){
    memcpy(sg, s, 32);
    build(); // calcul des orbites des triplets stockées ds
	// ochoice[MAXORB][3][MAXK]
    memset(cov, 0, sizeof(cov));
    cnt=0;
    for(int o=0;o<norbs;o++) { // on verifie que tous les triplet sont valides pour ce sigma
		if(onc[o]==0) return 0;
	}
	// on va developper et compter les fixes
    bt(0);
    return cnt;
}
// a partir de la succession des cycles (part) et du nombre de cyles (n)
// on construit la permutation correspondante
void mksig(int *part, int n, int s[8]){
    for(int i=0;i<8;i++) s[i]=i; // identite
    int pos=0;
    for(int i=0;i<n;i++){ // on balaye les cycles
        int l=part[i];
        if(l>1){ // si cycle de longueur > 1 , on implante le cycle
            for(int j=0;j<l-1;j++) s[pos+j]=pos+j+1;
            s[pos+l-1]=pos;
        }
        pos+=l;
    }
}

int main(){
    int types[NB_PART][8]={
#if defined(Ident)
	{1,1,1,1,1,1,1,1},
#endif
        {2,1,1,1,1,1,1},{2,2,1,1,1,1},{3,1,1,1,1,1},
        {2,2,2,1,1},{3,2,1,1,1},{4,1,1,1,1},{2,2,2,2},
        {3,2,2,1},{3,3,1,1},{4,2,1,1},{5,1,1,1},
        {3,3,2},{4,2,2},{4,3,1},{4,4},
        {5,2,1},{6,1,1},{5,3},{6,2},{7,1},{8}
    };
    int npart[NB_PART]={
#if defined(Ident)
		8,
#endif
		7,6,6,5,5,5,4,4,4,4,4,3,3,3,3,3,3,2,2,2,1};
    const char *L[NB_PART]={
#if defined(Ident)
		"(1,1,1,1,1,1,1,1)",
#endif
        "(2,1,1,1,1,1,1)","(2,2,1,1,1,1)","(3,1,1,1,1,1)",
        "(2,2,2,1,1)","(3,2,1,1,1)","(4,1,1,1,1)","(2,2,2,2)",
        "(3,2,2,1)","(3,3,1,1)","(4,2,1,1)","(5,1,1,1)",
        "(3,3,2)","(4,2,2)","(4,3,1)","(4,4)",
        "(5,2,1)","(6,1,1)","(5,3)","(6,2)","(7,1)","(8,)"
    };
    for(int i=0;i<NB_PART;i++){
        int s[8];
        mksig(types[i], npart[i], s); // on construit la permutation
#if defined DEBUG
		for(int i =0;i<8;i++) printf("%d->%d ",i,s[i]);
		printf("\n");
#endif
        long long f=fix(s); //
        printf("Fix%-18s = %10lld\n", L[i], f);
        fflush(stdout);
    }
    return 0;
}

