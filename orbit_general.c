/*
 * A) Initialisation de certains tableaux par init_basic();
 *   - tableau des 56 paires P2[28][2]={i,j} par ordre alphabetique i < j
 *     fonction reciproque pidx[j][i] = pidx[i][j] = np (rang dans P2)
 *   - tableau des triplets tri[n_tri] ={a,b,c} a < b < c  (ordre alphabetique)
 *     fonction reciproque tidx[a][b][c] = n_tri fonction riche pas de contrainte sur l'ordre de a,b,c
 * B) Initialisation de ST1 et ST2 (Stage) decrivant le level 1 et le level 2
 *   - par build_ST1 et build_ST2
 *   - pour le level 1 :
 *     - 720 permutations de {0,1,2,3,4,5}
 *     - les 20=6x5x4/6 triplets de {0,1,2,3,4,5} sont pris (ds W1)
 *   - pour le level 2 :
 *     - 2 permutations de {6,7}
 *     - les triplets sont donnés par la macro LISTE_W2
 *        - doit respecter chaque triplet doit contenir un 6 ou un 7 ou les deux (hors de W1
 *        - doit etre stable par la permutation 6<->7
 *   - pour les 2 levels on precalcule pour chaque permutation :
 *        - recip_tri[num Triplet] : la fonction reciproque sur les triplets
*            recip[tidx[perm(a)][perm(b)][perm(c)] = tidx[a][b][c]
*            Cela marche pour ST2 car les triplets de la liste W2 doivent etre stable par le groupe des permutations
*         - direct_pair[num Pair] :
*           direct_pair[pidx[a][b]] = pidx[perm(a)][perm(b)]
*         Exemple pour ST1 et la permutation (0->1->2->4->0).
*             Le triplet (0,1,2) devient (1,2,4)
*             - (0,1,2) (n° 0) devient (1,2,4) (n° 22) : donc recip[22] = 0
*             -
 *
 *
 * Meme methode que orbit_count.c, avec des groupes plus gros :
 *
 * Etage 1 : W1 = les 20 triplets de {0,...,5,6}, groupe G1 = S6 (720 perms,
 *           6 et 7 fixes).
 * Etage 2 : W2 = les 15 triplets {a,b,7}, a,b=0..6 groupe G2 = S2 = {id1,id2, 7)}.
 *           G2 fixe point par point toute affectation de W1.
 *
 *   N = somme sur les representants r1 (orbites de G1 sur les affectations
 *       faisables de W1) de  |O1| *
 *       somme sur les representants r2 (orbites de G2 sur les affectations
 *       faisables de W2) de  |O2| * count(r1, r2)
 *
 * Facteur nominal 720 * 15 = 10800. Le comptage est exact (orbites courtes
 * ponderees par leur vraie taille).
 *
 * Canonicalisation : le representant d'une orbite est l'affectation
 * lexicographiquement minimale (chiffres 0,1,2 dans l'ordre des triplets de W).
 * Le test est paresseux : comparaison chiffre par chiffre avec sortie
 * anticipee, et la derniere permutation "tueuse" est essayee en premier.
 *
 * UTILISATION
 *   gcc -O3 -march=native -o orbit_count6 orbit_count6.c
 *   ./orbit_count6 --info                  structure (construit les orbites)
 *   ./orbit_count6 --sample K              estime le temps total et N avec K jobs
 *   ./orbit_count6                         tous les jobs (fichier res6v2_0.txt)
 *   ./orbit_count6 PART NPARTS             jobs d'id % NPARTS == PART
 *                                          (fichier res6v2_PART.txt), reprise auto
 *   ./orbit_count6_v2 --sum res6v2_0.txt ...    total et jobs manquants
 *
 * Un job = un representant r1 de l'etage 1 ; il boucle sur les representants
 * de l'etage 2 et ecrit une ligne "id poids1 T" avec T = somme poids2*count.
 * Resultat attendu : N = 125 047 355 835 960.
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef unsigned long long u64;

typedef int8_t typeRemCov ;
#define nbPairs 28
#define nbTriplets 56

#define CIBLE 19
#include "cible.h"
/* ---------------------------------------------------------------- base */
static typeRemCov pidx[8][8]; // tableau des index[0,28[ des paires ((i,j) et (j,i) => index )
static typeRemCov P2[nbPairs][2]; // liste des paires
static typeRemCov tri[nbTriplets][3]; // liste des triplets
static typeRemCov tidx[8][8][8];//  tableau des index [0,56[ des triplets ordonnés (a<b<c => index )
static typeRemCov tp[nbTriplets][3]; // liste des index de paires associées au triplet pid[a][b], pid[a][c], pid[b][c]
/* ------------------------------------------------ backtracking */


typedef typeRemCov remT_Type[nbRemT][3];




/************************* structures **********************************/
typedef struct {
    typeRemCov pairs[MAXT] ;
    int16_t weight;
} Leaf;

typedef struct {
    typeRemCov recip_tri[MAXT] ;
    typeRemCov direct_pair[nbPairs]  ;
}Perm;
typedef struct {
    typeRemCov nbT;
    int nbPerm;                          // nbT : nombre de triplets, nbPerm : nombre de permutations.
    const typeRemCov *W;                       // liste des triplets (invariant)
    typeRemCov cov[nbPairs];                        // couverture des paires
    typeRemCov curPair[MAXT];
    Leaf *leafs; long nleaf, cap;
    long wAssign;                       // poids total du developpement
    int killer;
    Perm Pe[NBPERM1] ;
} Stage;

static void init_basic(void) ;
static Stage * build_ST1(void) ;
static Stage * build_ST2(void) ;
static Stage *stage_build(Stage * S,int perms[][8]) ;
void stage_phase1(Stage *S) ;
//void build_remT(const Stage *ST1,const Stage *ST2) ;
remT_Type * build_remT(const Stage *ST1,const Stage *ST2) ;
const typeRemCov *   build_rem0(const remT_Type * remT) ;

void check_sum(long nb_leaf,char *name_fp ) ;
static void print_info(FILE * fout,const Stage *ST1,const Stage *ST2) ;
void compute_sample( int K,const typeRemCov *rem0,const Stage *ST1,const Stage *ST2) ;

static int gen_group(int lo, int m, int perms[][8]) ;

void compute_jobs( const char *fileName,const typeRemCov *rem0,const Stage *ST1,const Stage *ST2) ;
remT_Type * get_remT(remT_Type * remT) {
    static remT_Type * initRemT = NULL ;
    if ( remT == NULL ) {
        return initRemT;
    } else {
        initRemT = remT ;
        return remT ;
    }
}
static void bt(int k) ;

int main(int argc, char **argv){
    const clock_t c0 = clock();
    init_basic();
    printf("Execution %s\n",CIBLE_NAME);
    // ************* construction deterministe   ************
    Stage * ST1 = build_ST1(); printf("buid_st1 OK\n");
    Stage * ST2 = build_ST2();  printf("buid_st2 OK\n");
    const typeRemCov * rem0 ;
    {
        remT_Type * remT = build_remT(ST1,ST2); // remT temporaire
        rem0 = build_rem0(remT); // calcul de rem0 a partir de remT
        get_remT(remT); // initialise get_remT avec pour remT pour bt
        bt(-1); // premier appel de bt qui utilise get_remT pour recuperer remT en lovcal
    }
    // ************* phase 1 pour ST2  puis ST1 ************
    // developpement des triplets par recursion de stage_dfs () => leafs
    // stage_leaf :
    //choix de killer la permutation representative  d'indice inferieur (stage_leaf)
    // calcul de weight= nbPerm/stab et choix a chaque triplet ds choice[MAXT]
    stage_phase1(ST1);  printf("ST1 phase 1 OK\n");
    stage_phase1(ST2);  printf("ST2 phase 1 OK\n");

    fprintf(stderr, "[orbites construites en %.1f s]\n", (double)(clock() - c0) / CLOCKS_PER_SEC);

    if (argc >= 2 && strcmp(argv[1], "--info") == 0) { print_info(stdout,ST1,ST2); return 0; }
    if (argc >= 3 && strcmp(argv[1], "--sum") == 0) { check_sum(ST1->nleaf,argv[2] ) ; return 0; }
    if (argc >= 3 && strcmp(argv[1], "--sample") == 0) { int K = atoi(argv[2]); compute_sample(K,rem0,ST1,ST2); return 0; }
    char fname[64];
    if ( argc > 1) { strcpy(fname, argv[1]); } else {
        sprintf(fname, "res_%s_%d.txt",CIBLE_NAME, 0);
    }
    // *********  phase 2 **************************
    // developpement des feuilles par developpement recursif sur les paires de bt()
    // suivant ls
    compute_jobs(fname,rem0,ST1,ST2) ;
    return 0;
}
/************ BEGIN ***************** phase 1 pour ST2 puis ST1 *****************************/
static void stage_dfs(Stage *S, int j ) ;
static void stage_leaf(Stage *S) ;
void stage_phase1(Stage *S) {
    stage_dfs(S, 0);
    {
        long sw = 0;
        for (long q = 0; q < S->nleaf; q++) sw += S->leafs[q].weight;
        if (sw != S->wAssign) {
            fprintf(stderr, "ERREUR: somme des poids %ld != %ld affectations\n", sw, S->wAssign);
            exit(1);
        }
    }
}

// algo recursif de developpement des triplets choisis par permutations
// i niveau courant du developpement (de 0 -> 20 ou 0 ->6
// appelle pour chaque terminaison stage_leaf pour calculer le poids et le stabilisateur
// utilise uniquement dans las phase 1
// developpe du niveau 0 au niveau nbT (nombre de triplets 20)
static void stage_dfs(Stage *S, int j )
{
    if (j == S->nbT) {
        stage_leaf(S);
        return;
    }
    const typeRemCov *pair_choice = tp[S->W[j]] ;
    // for (typeRemCov c = 0; c < 3; c++) {
    typeRemCov p = pair_choice[0]; // paire
    if (S->cov[p] < 2) {
        S->cov[p]++;  S->curPair[j] = p ;
        stage_dfs(S, j + 1);
        S->cov[p]--;
    }
    p = pair_choice[1]; // paire
    if (S->cov[p] < 2) {
        S->cov[p]++;  S->curPair[j] = p ;
        stage_dfs(S, j + 1);
        S->cov[p]--;
    }
    p = pair_choice[2]; // paire
    if (S->cov[p] < 2) {
        S->cov[p]++;  S->curPair[j] = p ;
        stage_dfs(S, j + 1);
        S->cov[p]--;
    }

    // }
}
// calcul du developpement des feuilles (terminaux des developpements de permutation)
// calcule dans wAssign le poids total (nombre de leaf)
static void stage_leaf(Stage *S)
{
    int  stab = 0;
    typeRemCov nT = S->nbT;
   S->wAssign++;
    { // tueuse d'abord
        Perm *peKiller = &S->Pe[S->killer];
        for (int j = 0; j < nT; j++) {
            typeRemCov dk = peKiller->direct_pair[S->curPair[peKiller->recip_tri[j]]];
            if (dk < S->curPair[j]) { return ; }
            if (dk > S->curPair[j]) { break ; }
        }
    }
    for (int nPerm = 0; nPerm < S->nbPerm; nPerm++) {  // boucle sur toutes les permutations du stage
        int less = 0, greater = 0;
        Perm *pe = &S->Pe[nPerm];
        for (int j = 0; j < nT; j++) {
            typeRemCov d = pe->direct_pair[S->curPair[pe->recip_tri[j]]];
            if (d < S->curPair[j]) { less = 1 ;break ; }
            if (d > S->curPair[j]) { greater = 1 ; break ; }
        }
        if (less) {
            S->killer = nPerm ;
             return;
        } //perm lexicalement inferieur => devient killer
        if (!greater) stab++;           /* choice identique : stabilisateur */
    }
    //printf("newRes=%d/%d",S->nbPerm,stab);
    /* canonique */
    if (S->nleaf == S->cap) {
        // allocation supplémentaire ?
        S->cap = S->cap ? S->cap * 2 : 1024;
        S->leafs = (Leaf *)realloc(S->leafs, (size_t)S->cap * sizeof(Leaf));
        if (!S->leafs) { fprintf(stderr, "memoire insuffisante\n"); exit(1); }
    }
    {
        Leaf *r = &S->leafs[S->nleaf++];
        memcpy(r->pairs, S->curPair, nT*sizeof(S->curPair[0]));
        r->weight = (short)(S->nbPerm / stab);
    }
}
/************ END ***************** phase 1 pour ST2 puis ST1 *****************************/

/**** BEGIN ************************** phase 2 (jobs) *****************************/
static u64 run_job(long id,const typeRemCov rem0[nbPairs],const Stage *ST1,const Stage *ST2) ;
static clock_t jobs_total = 0;
static clock_t bt_total = 0;
void compute_jobs( const char *fileName,const typeRemCov *rem0,const Stage *ST1,const Stage *ST2)   {
    long nb_leaf = ST1->nleaf;
    char fname[64];
    time_t t0 = time(NULL);
    long ndo = 0, nskip = 0;
    if (fileName != NULL) {
        strcpy(fname, fileName);
    } else {
        sprintf(fname, "res_%s_%d.txt",CIBLE_NAME,0);
    }
    char *done = (char *)calloc((size_t)nb_leaf, 1);
    FILE *fp = fopen(fname, "r");
    if (fp) {
        long jid; int w; u64 T;
        char line[200] ;
        while (fgets(line,200,fp)) {
            if (line[0] != '#') {
                if (sscanf(line, "%ld %d %llu", &jid, &w, &T) == 3) {
                    if (jid >= 0 && jid < nb_leaf) {
                        done[jid] = 1; nskip++;
                    }
                }
            }
        }
        fclose(fp);
    }
    print_info(stdout,ST1,ST2);
    fp = fopen(fname, "a");
    if (!fp) { perror("res"); return ; }
    print_info(fp,ST1,ST2);
    printf(" fichier %s, deja faits : %ld\n", fname, nskip);
    int freq_ndo = (nb_leaf / 100) +1  ;
    int mod_ndo = 0;
    for (long id = 0; id < nb_leaf; id ++) {
        // for (long id = nb_leaf-1; id >= 0; id --) {
        if (done[id]) continue;
        u64 T = run_job(id,rem0,ST1,ST2);
        fprintf(fp, "%ld %d %llu\n", id, (int)ST1->leafs[id].weight, T);
        mod_ndo++;
        if (mod_ndo == freq_ndo) {
            fflush(fp);
            printf("[%lds](%.3f,%.3f)T=%lld %ld jobs faits dans cette session (dernier id %ld)\n",
                   (long)(time(NULL) - t0),jobs_total /(double) CLK_TCK,bt_total/(double) CLK_TCK,T, ndo+mod_ndo, id);
            ndo += mod_ndo;
            mod_ndo = 0 ;
        }
    }
    fprintf(fp,"#[%lds] %ld jobs faits dans cette session %s (dernier id %ld)\n",
           (long)(time(NULL) - t0), ndo,CIBLE_NAME, nb_leaf);
    fclose(fp);
    check_sum(nb_leaf,fname ) ;
}
// donnees partagées par run_job() et bt()
static typeRemCov  cov[nbPairs] ;
static typeRemCov rem_cov[nbPairs] ;
static u64 cnt ;

// un job = un representant de l'etage 1 ; retourne T = somme w2*count
// rem0 contient les paires disponibles pour les triplets non developpés
// rem0 est fixe
// les triplets developpés sont fixe : ceux dans ST1 (20 si permut [01..6[) et ceux ds la liste de ST2
static u64 run_job(long id,const typeRemCov rem0[nbPairs],const Stage *ST1,const Stage *ST2)
{
    clock_t deb = clock() ;
    typeRemCov cov1[nbPairs];
    u64 T = 0;
    const Leaf *r1 = &ST1->leafs[id];
    memset(cov1, 0, sizeof(cov1));
    for (int i = 0; i < MAXT; i++) {  // contrib stage 1
        cov1[r1->pairs[i]]++ ;
    }
    for (int q = 0; q < ST2->nleaf; q++) {
        const Leaf *r2 = &ST2->leafs[q];
        int ok = 1;
        memcpy(cov, cov1, sizeof(cov1)); // init stage1
        for (typeRemCov i = 0; i < SIZE_W2; i++) {   // add stage2
            typeRemCov p = r2->pairs[i] ;
            if (cov[p]>=2 ) { ok = 0; break ;}
            cov[p]++;
        }
        if (ok) {
            for (typeRemCov p =0; p< nbPairs;  p++) {
                if (cov[p] + rem0[p] < 2) { ok = 0; break; }
                rem_cov[p] = cov[p] + rem0[p] ;
            }
        }
        if (!ok) continue;
        clock_t deb_bt = clock() ;
        cnt = 0;
        bt(0);
        T += (u64)r2->weight * cnt;
        clock_t end_bt = clock() ;
        bt_total += end_bt - deb_bt;
    }
    clock_t end = clock();
    jobs_total += end- deb ;
    return T;
}
// recursif, developpe les triplets complementaires (nbRemT) en phase 2 apres ceux developpes par les permutations
// va de niveau 0 a nbRemT (56-(20+6)
// calcule cnt le nombre de developpements possibles.
// Utilise dans la phase 2 (jobs)
static int bt_init = 1 ;
static remT_Type remT_bt = {0};
static void bt(const int k)
{
    if (bt_init) {
        bt_init = 0 ;
        memcpy(&remT_bt,get_remT(NULL),sizeof(remT_Type)) ;
        return ;
    }
    if (k == nbRemT) { cnt++; return ; }
    const typeRemCov pa = remT_bt[k][0]; const typeRemCov  pb = remT_bt[k][1]; const typeRemCov  pc = remT_bt[k][2];
    if (rem_cov[pa] >= 2 && rem_cov[pb] >= 2 &&  rem_cov[pc] >= 2) {
        rem_cov[pa]-- ; rem_cov[pb]-- ; rem_cov[pc]-- ;
        if (cov[pa] < 2 && rem_cov[pb] >= 2 && rem_cov[pc] >= 2) { // pas encore de double attribution
                cov[pa]++;rem_cov[pa]++; bt(k + 1); cov[pa]--;rem_cov[pa]--;
        }
         if (cov[pb] < 2 && rem_cov[pa] >= 2 && rem_cov[pc] >= 2) { // pas encore de double attribution
            cov[pb]++;rem_cov[pb]++; bt(k + 1); cov[pb]--;rem_cov[pb]--;
        }
        if (cov[pc] < 2 && rem_cov[pa] >= 2 && rem_cov[pb] >= 2) { // pas encore de double attribution
            cov[pc]++;rem_cov[pc]++; bt(k + 1); cov[pc]--;rem_cov[pc]--;
        }
        rem_cov[pa]++ ; rem_cov[pb]++ ; rem_cov[pc]++ ;
    }
}
/***** END ************************* phase 2 (jobs) *****************************/

/*************************************** routines tritement options de main() **************/
/* -------------------------------------------------------- main */
static void print_info(FILE * fout,const Stage *ST1,const Stage *ST2) {
    long  tw1 = 0, tw2 = 0, full = 0;
    for (long i = 0; i < ST1->nleaf; i++) { tw1 += ST1->leafs[i].weight; if (ST1->leafs[i].weight == ST1->nbPerm) full++; }
    for (long i = 0; i < ST2->nleaf; i++) tw2 += ST2->leafs[i].weight;
    fprintf(fout,"# Etage 1 : %ld affectations faisables de W1 (%d triplets), %ld orbites sous S5 "
           "(somme poids = %ld, orbites pleines de %d : %ld)\n", ST1->wAssign,ST1->nbT, ST1->nleaf, tw1,ST1->nbPerm, full);
    fprintf(fout,"# Etage 2 : %ld affectations faisables de W2 (%d triplets), %ld orbites sous S2 "
           "(somme poids = %ld)\n", ST2->wAssign,ST2->nbT, ST2->nleaf, tw2);
    fprintf(fout,"# Nombre de jobs (representants de l'etage 1) = %ld ; sous-comptages par job = %ld\n",
           ST1->nleaf, ST2->nleaf);
    fprintf(fout,"# Gain nominal sur les configurations de depart : %.1f\n",
           (double)ST1->wAssign * ST2->wAssign / ((double)ST1->nleaf * ST2->nleaf));
    fprintf(fout,"# Triplets restants a parcourir : %d\n", nbRemT);
    fflush(fout);
}
void check_sum(long nb_leaf,char *name_fp ) {
    char *done = calloc((size_t)nb_leaf, 1);
    unsigned long long total = 0;
    long ndone = 0, missing = 0;

    FILE *fp = fopen(name_fp, "r");
    long jid; int w; u64 T;
    if (!fp) {
        fprintf(stderr, "ouverture impossible : %s\n", name_fp);
    }else {
        char line[200] ;
        while (fgets(line,200,fp)) {
            if (line[0] != '#') {
                if (sscanf(line, "%ld %d %llu", &jid, &w, &T) == 3) {
                    if (jid >= 0 && jid < nb_leaf && !done[jid]) {
                        done[jid] = 1; ndone++;
                        total += (u64)w * T;
                    }
                }
            }
        }
    }
    fclose(fp);

    for (long id = 0; id < nb_leaf; id++) if (!done[id]) missing++;
    printf("jobs termines : %ld / %ld   manquants : %ld\n", ndone, nb_leaf, missing);
    printf("TOTAL (partiel si manquants) = %llu\n", total);
    if (missing == 0)
        printf("Comparaison avec N attendu 125047355835960 : %s\n",
               total == 125047355835960ULL ? "IDENTIQUE" : "DIFFERENT");
    free(done);

}
void compute_sample( int K,const typeRemCov *rem0,const Stage *ST1,const Stage *ST2) {
    long nb_leaf = ST1->nleaf ;
    long step = nb_leaf / K;
    double t0 = (double)clock() / CLOCKS_PER_SEC;
    double sumT = 0.0;
    for (int i = 0; i < K; i++) {
        long id = (long)i * step + step / 2;
        sumT += (double)ST1->leafs[id].weight * (double)run_job(id,rem0,ST1,ST2);
    }
    double dt = (double)clock() / CLOCKS_PER_SEC - t0;
    printf("%d jobs echantillonnes sur %ld : %.2f s  (%.4f s/job)\n", K, nb_leaf, dt, dt / K);
    printf("Temps total estime (1 coeur) : %.0f s = %.1f h\n", dt / K * nb_leaf, dt / K * nb_leaf / 3600.0);
    printf("Estimation de N (tres approximative) : %.3e\n", sumT / K * nb_leaf);

}

/*********************** routines de construction *****************************************/
static Stage * build_ST1(void) {
    static int perms1[NBPERM1][8] ; // liste des NBPERM permutation de [0..PERM1[
    static typeRemCov  W1[MAXT] ; // liste des triplets de [0,PERM1[
    // 012 013 014 023 024 034 123 124 134 234


    // int k = 0;
    for (int k=0, a = 0; a < PERM1; a++)
        for (int b = a + 1; b < PERM1; b++)
            for (int c = b + 1; c < PERM1; c++) W1[k++] = tidx[a][b][c];
    Stage *ST1 = (Stage *)calloc(1, sizeof(Stage));
    int nbPerm1 = gen_group(0, PERM1, perms1); // permutations de {0,1,..,PERM1-1}
    ST1->nbT = MAXT; ST1->nbPerm = nbPerm1; ST1->W = W1;
    ST1 = stage_build(ST1, perms1);
    return ST1;
}
static Stage * build_ST2(void) {
    const  int8_t  liste_w2[SIZE_W2][3]=LIST_W2 ;
    static typeRemCov  W2[SIZE_W2] ; // liste des triplets de W2
    for (int k=0;k<SIZE_W2;k++) {
        W2[k]= tidx[liste_w2[k][0]][liste_w2[k][1]][liste_w2[k][2]];
    }
    static int perms2[NBPERM2][8];  // liste des 2 permutations 01234567 01234576
    Stage *ST2 = (Stage *)calloc(1, sizeof(Stage));
    int nbPerm2 = gen_group(PERM2,G_PERM2,perms2);
    if (nbPerm2 != NBPERM2) {
        fprintf(stderr,"ERROR nbPerm2=%d expexted %d\n",nbPerm2,NBPERM2);
        exit(1);
    }
    ST2->nbT =SIZE_W2; ST2->nbPerm = nbPerm2; ST2->W = W2;
    ST2 = stage_build(ST2, perms2);
    return ST2 ;
}
// *W : Liste de triplets
// n : nombre de triplets
// perms : liste de permutations [0,8[
// ng : nombre de permutations
// retourne la structure Stage correspondante
// applique au niveau 1 (SIZE_W1 20 triplets) et au niveau 2 (SIZE_W2 triplets)
static Stage *stage_build(Stage * S,int perms[][8]) {
    // const typeRemCov *W = S->W ;
    for (int nP = 0; nP < S->nbPerm; nP++) { // boucle sur les permutations
        Perm *pe = &S->Pe[nP];
        for (typeRemCov i = 0; i < S->nbT; i++) { // boucle sur les triplets
            typeRemCov *t = tri[S->W[i]]; // on recupere le triplet
            typeRemCov gi = tidx[perms[nP][t[0]]][perms[nP][t[1]]][perms[nP][t[2]]]; // on applique la permutation
            // et l'on recupere l'index correspondant du triplet associe
            typeRemCov j ;
            for (j = 0; j < S->nbT; j++) if (S->W[j] == gi) break; // on recupere l'index du triplet associe
            if (j == S->nbT) { // ne doit jamais se produire
                fprintf(stderr, "W non invariant sous le groupe pour {%d,%d,%d}\n"
                    ,perms[nP][t[0]],perms[nP][t[1]],perms[nP][t[2]]);
                exit(1);
            }
            pe->recip_tri[j] = i; // i => j , on le note ds la correpondance inversée
        }
        for (typeRemCov i= 0; i<nbPairs;i++) {
            typeRemCov j =  pidx[perms[nP][P2[i][0]]] [perms[nP][P2[i][1]]] ;;
            pe->direct_pair[i] = j ;
        }
    }
    return S;
}

// const remT_Type * build_remT(const Stage *ST1,const Stage *ST2) {
remT_Type * build_remT(const Stage *ST1,const Stage *ST2) {
    remT_Type * remT = calloc(1, sizeof(remT_Type));
    int used[nbTriplets] = {0};
    for (int i = 0; i < MAXT; i++) used[ST1->W[i]] = 1;
    for (int i = 0; i < SIZE_W2; i++) used[ST2->W[i]] = 1;
    for (int j = 0; j < nbRemT; j++) {
        for (int i = 0; i < nbTriplets; i++){
            if (!used[i]) {
                used[i] = 1;
                for (uint8_t k = 0; k < 3; k++) {
                    (*remT)[j][k]  = tp[i][k];
                }
                break ;
            }
        }
    }
    return remT ;
}
static typeRemCov l_rem0[nbRemT];
const typeRemCov *   build_rem0(const remT_Type * remT0) {
    for (int i = 0; i < nbRemT; i++) { l_rem0[(*remT0)[i][0]]++; l_rem0[(*remT0)[i][1]]++; l_rem0[(*remT0)[i][2]]++; }
    return l_rem0 ;
}
/****************************** outils divers ********************************************/


// retourne 0 ou 1 suivant que la permutation de n elements dans a[] est la derniere ou non
// mais dans a[] la permutation suivante si ce n'est pas la derniere
// "ordre alphabeitque"
// la premiere correspond a a[0] < a[1] < a[2], ...,< a[n-1]
// la derniere correspond à a[0] > a[1] > a[2], ...,>a[n-1]
// exemple si a[]={'a','b','c','d'}
// abcd abdc acbd acdb adbc adcb bacd badc bcad bcda bdac bdca cabd cadb cbad cbda cdab cdba dabc dacb dbac dbca dcab dcba
static int next_perm(int *a, int n)
{
    int i = n - 2;
    while (i >= 0 && a[i] > a[i + 1]) i--;
    if (i < 0) return 0;
    int j = n - 1;
    while (a[j] < a[i]) j--;
    { int t = a[i]; a[i] = a[j]; a[j] = t; }
    for (j = i + 1, i = n - 1; j < i; j++, i--) { int t = a[j]; a[j] = a[i]; a[i] = t; }
    return 1;
}

static int gen_group(int lo, int m, int perms[][8])
{
    int a[8], ng = 0;
    for (int i = 0; i < m; i++) a[i] = lo + i;
    do {
        for (int i = 0; i < 8; i++) perms[ng][i] = i;
        for (int i = 0; i < m; i++) perms[ng][lo + i] = a[i];
        ng++;
    } while (next_perm(a, m));
    return ng;
}
static void init_basic(void)
{
    typeRemCov np = 0, nt = 0;
    for (typeRemCov i = 0; i < 8; i++)
        for (typeRemCov j = i + 1; j < 8; j++) {
            pidx[i][j] = pidx[j][i] = np;
            P2[np][0] = i; P2[np][1] = j; np++;
        }
    for (int a = 0; a < 8; a++)
        for (int b = a + 1; b < 8; b++)
            for (int c = b + 1; c < 8; c++) {
                tri[nt][0] = a; tri[nt][1] = b; tri[nt][2] = c;
                tidx[a][b][c] = nt;
                tidx[a][c][b] = nt;
                tidx[b][a][c] = nt;
                tidx[b][c][a] = nt;
                tidx[c][a][b] = nt;
                tidx[c][b][a] = nt;

                tp[nt][0] = pidx[a][b];
                tp[nt][1] = pidx[a][c];
                tp[nt][2] = pidx[b][c];
                nt++;
            }
}
