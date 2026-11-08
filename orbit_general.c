/*
 * orbit_count7_v1.c
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

typedef int32_t typeRemCov ;
#define nbPairs 28
#define nbTriplets 56
/* ---------------------------------------------------------------- base */
static int pid_[8][8]; // tableau des index[0,28[ des paires ((i,j) et (j,i) => index )
static int P2[nbPairs][2]; // liste des paires
static int tri[nbTriplets][3]; // liste des triplets
static int tidx[8][8][8];//  tableau des index [0,56[ des triplets ordonnés (a<b<c => index )
static int tp[nbTriplets][3]; // liste des index de paires associées au triplet pid[a][b], pid[a][c], pid[b][c]

// construit les tableaux precedents
static void init_basic(void)
{
    int np = 0, nt = 0;
    for (int i = 0; i < 8; i++)
        for (int j = i + 1; j < 8; j++) {
            pid_[i][j] = pid_[j][i] = np;
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

                tp[nt][0] = pid_[a][b];
                tp[nt][1] = pid_[a][c];
                tp[nt][2] = pid_[b][c];
                nt++;
            }
}
// retrouve l'index d'un triplet (a,b,c) non ordonné.
static int tindex3(int a, int b, int c)
{
    int t;
    if (a > b) { t = a; a = b; b = t; }
    if (b > c) { t = b; b = c; c = t; }
    if (a > b) { t = a; a = b; b = t; }
    return tidx[a][b][c];
}

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
// perm5_W4 CIBLE = 1
// perm6_W6 CIBLE = 2

#define CIBLE 20

#if CIBLE == 1
#   define CIBLE_NAME  "perm5_W4"
#   define PERM1 5  // permutation de [012..4[
#   define NBPERM1 120 //5!
#   define MAXT   10  // 5x4x3
#   define SIZE_W2 4
#   define LIST_W2 {{5,6,7},{0,5,6},{0,5,7},{0,6,7}}
#   define PERM2   5
#   define G_PERM2 3
#   define NBPERM2 6
#   define nbRemT 42  // NbTriplets - MAXT - SIZE_W2  56-(10+4)
#elif CIBLE == 2
#   define CIBLE_NAME  "perm6_W6"
#   define PERM1   6  // permutation de [012..4[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 6
#   define LIST_W2 {{0,6,7},{1,6,7},{2,6,7},{3,6,7},{4,6,7},{5,6,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 30  // NbTriplets - MAXT - SIZE_W2  56-(20+6)
#elif CIBLE == 3
#   define CIBLE_NAME  "perm6_W12"
#   define PERM1   6  // permutation de [012..4[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 12
#   define LIST_W2 {{0,1,6},{0,2,6},{0,3,6},{0,1,7},{0,2,7},{0,3,7},{1,2,6},{1,3,6},{1,2,7},{1,3,7},{2,3,6},{2,3,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 24  // NbTriplets - MAXT - SIZE_W2  56-(20+12)
#elif CIBLE == 4
#   define CIBLE_NAME  "perm6_W10"
#   define PERM1   6  // permutation de [012..4[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 10
#   define LIST_W2 {{0,1,6},{0,2,6},{0,3,6},{0,1,7},{0,2,7},{0,3,7},{1,2,6},{1,3,6},{1,2,7},{1,3,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 26  // NbTriplets - MAXT - SIZE_W2  56-(20+10)
#elif CIBLE == 5
#   define CIBLE_NAME  "perm6_W8"
#   define PERM1   6  // permutation de [012..4[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 8
#   define LIST_W2 {{0,1,6},{0,2,6},{0,3,6},{0,1,7},{0,2,7},{0,3,7},{1,2,6},{1,2,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 28  // NbTriplets - MAXT - SIZE_W2  56-(20+8)
#elif CIBLE == 6
#   define CIBLE_NAME  "perm6_W10a"
#   define PERM1   6  // permutation de [012..4[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 10
#   define LIST_W2 {{5,6,7},{0,5,6},{0,5,7},{0,6,7},{1,5,6},{1,5,7},{1,6,7},{2,5,6},{2,5,7},{2,6,7}}
#   define PERM2   5
#   define G_PERM2 3
#   define NBPERM2 6
#   define nbRemT 26  // NbTriplets - MAXT - SIZE_W2  56-(20+10)
#elif CIBLE == 7
#   define CIBLE_NAME  "perm6_W10b"
#   define PERM1   6  // permutation de [012..6[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 10
#   define LIST_W2 {{5,6,7},{0,5,6},{0,5,7},{0,6,7},{2,5,6},{2,5,7},{2,6,7},{4,5,6},{4,5,7},{4,6,7}}
#   define PERM2   5
#   define G_PERM2 3
#   define NBPERM2 6
#   define nbRemT 26  // NbTriplets - MAXT - SIZE_W2  56-(20+10)
#elif CIBLE == 8
#   define CIBLE_NAME  "perm6_W10c"
#   define PERM1   6  // permutation de [012..4[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 10
#   define LIST_W2 {{5,6,7},{0,5,6},{0,5,7},{0,6,7},{1,5,6},{1,5,7},{1,6,7},{2,5,6},{2,5,7},{2,6,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 26  // NbTriplets - MAXT - SIZE_W2  56-(20+10)
#elif CIBLE == 9
#   define CIBLE_NAME  "perm5_W10d"
#   define PERM1   5  // permutation de [012..4[
#   define NBPERM1 120 // 5!
#   define MAXT    10 // 5x4x3/6
#   define SIZE_W2 10
#   define LIST_W2 {{5,6,7},{0,5,6},{0,5,7},{0,6,7},{1,5,6},{1,5,7},{1,6,7},{2,5,6},{2,5,7},{2,6,7}}
#   define PERM2   5
#   define G_PERM2 3
#   define NBPERM2 6
#   define nbRemT 36  // NbTriplets - MAXT - SIZE_W2  56-(10+10)
#elif CIBLE == 10
#   define CIBLE_NAME  "perm6_W4"
#   define PERM1   6  // permutation de [012..4[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 4
#   define LIST_W2 {{5,6,7},{0,5,6},{0,5,7},{0,6,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 32  // NbTriplets - MAXT - SIZE_W2  56-(20+4)
#elif CIBLE == 11
#   define CIBLE_NAME  "perm6_W4a"
#   define PERM1   6  // permutation de [012..5[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 4
#   define LIST_W2 {{0,1,6},{0,1,7},{2,3,6},{2,3,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 32  // NbTriplets - MAXT - SIZE_W2  56-(20+4)
#elif CIBLE == 12
#   define CIBLE_NAME  "perm6_W14"
#   define PERM1   6  // permutation de [012..5[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 14
#   define LIST_W2 {{0,1,6},{0,2,6},{0,3,6},{0,4,6},{0,1,7},{0,2,7},{0,3,7},{0,4,7},{1,2,6},{1,3,6},{1,2,7},{1,3,7},{2,3,6},{2,3,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 22  // NbTriplets - MAXT - SIZE_W2  56-(20+14)
#elif CIBLE == 13
#   define CIBLE_NAME  "perm6_W4e"
#   define PERM1   6  // permutation de [012..5[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 4
#   define LIST_W2 {{0,4,6},{0,4,7},{1,5,6},{1,5,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 32  // NbTriplets - MAXT - SIZE_W2  56-(20+4)
#elif CIBLE == 14
#   define CIBLE_NAME  "perm6_W2"
#   define PERM1   6  // permutation de [012..5[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 2
#   define LIST_W2 {{0,4,6},{0,4,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 34  // NbTriplets - MAXT - SIZE_W2  56-(20+2)
#elif CIBLE == 15
#   define CIBLE_NAME  "perm6_W8b"
#   define PERM1   6  // permutation de [012..4[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 8
#   define LIST_W2 {{0,6,7},{0,1,6},{0,1,7},{0,2,6},{0,2,7},{1,2,6},{1,2,7},{1,6,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 28  // NbTriplets - MAXT - SIZE_W2  56-(20+8)

#elif CIBLE == 16
#   define CIBLE_NAME  "perm6_W12a"
#   define PERM1   6  // permutation de [012..4[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 12
#   define LIST_W2 {{0,1,6},{0,2,6},{0,3,6},{0,1,7},{0,2,7},{0,3,7},{1,2,6},{1,2,7},{0,6,7},{1,6,7},{2,6,7},{3,6,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 24  // NbTriplets - MAXT - SIZE_W2  56-(20+12)
#elif CIBLE == 17
#   define CIBLE_NAME  "perm6_W12b"
#   define PERM1   6  // permutation de [012..4[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 12
#   define LIST_W2 {{0,1,6},{0,1,7},{0,2,6},{0,2,7},{0,3,6},{0,3,7},{1,2,6},{1,2,7},{1,3,6},{1,3,7},{2,3,6},{2,3,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 24  // NbTriplets - MAXT - SIZE_W2  56-(20+12)
#elif CIBLE == 18
#   define CIBLE_NAME  "perm6_W4f"
#   define PERM1 6  // permutation de [012..4[
#   define NBPERM1 720 //5!
#   define MAXT   20  // 6x5x4/3
#   define SIZE_W2 4
#   define LIST_W2 {{0,4,6},{0,4,7},{2,5,6},{2,5,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 32  // NbTriplets - MAXT - SIZE_W2  56-(20+4)
#elif CIBLE == 19
#   define CIBLE_NAME  "perm6_W8c"
#   define PERM1   6  // permutation de [012..4[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 8
#   define LIST_W2 {{0,1,6},{0,1,7},{0,2,6},{0,2,7},{1,2,6},{1,2,7},{0,3,6},{0,3,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 28  // NbTriplets - MAXT - SIZE_W2  56-(20+8)
#elif CIBLE == 20
#   define CIBLE_NAME  "perm6_W9"
#   define PERM1   6  // permutation de [012..4[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 9
#   define LIST_W2 {{0,6,7},{1,6,7},{2,6,7},{0,1,6},{0,1,7},{0,2,6},{0,2,7},{1,2,6},{1,2,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 27  // NbTriplets - MAXT - SIZE_W2  56-(20+9)
#elif CIBLE == 21
#   define CIBLE_NAME  "perm6_W12c"
#   define PERM1   6  // permutation de [012..4[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 12
#   define LIST_W2 {{0,6,7},{1,6,7},{2,6,7},{3,6,7},{0,1,6},{0,1,7},{0,2,6},{0,2,7},{0,3,6},{0,3,7},{1,2,6},{1,2,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 24  // NbTriplets - MAXT - SIZE_W2  56-(20+12)
#elif CIBLE == 22
#   define CIBLE_NAME  "perm6_W1"
#   define PERM1   6  // permutation de [012..4[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 8
#   define LIST_W2 {{0,1,6},{0,2,6},{0,3,6},{0,4,6},{1,2,6},{1,3,6},{1,4,6},{1,5,6}}
#   define PERM2   6
#   define G_PERM2 1
#   define NBPERM2 1
#   define nbRemT 28  // NbTriplets - MAXT - SIZE_W2  56-(20+8)

#endif

/* ------------------------------------------------------- etages */
// #define MAXT 10
typedef struct {
    uint8_t choice[MAXT];  // liste des MAXT(35) choix pour la feuille
    int16_t weight;
} Leaf;

typedef struct {
    uint8_t recip[MAXT] ;
    uint8_t choiceByP[MAXT][3] ;
}Perm;
typedef struct {
    int nbT, nbPerm;                          // nbT : nombre de triplets, nbPerm : nombre de permutations.
    const int *W;                       // liste des triplets (invariant)
    int8_t cov[nbPairs];                        // couverture des paires
    uint8_t curChoice[MAXT];             // choix pour chaque triplet (a,b,c) de la paire associée (0 = (a,b) , 1=(a,c) 2=(b,c)
    Leaf *leafs; long nleaf, cap;
    long wAssign;                       // poids total du developpement
    int killer;
    Perm Pe[NBPERM1] ;
} Stage;

// calcul du developpement des feuilles (terminaux des developpements de permutation)
// calcule dans wAssign le poids total (nombre de leaf)
static void stage_leaf(Stage *S)
{
    int  stab = 0, nT = S->nbT;
    const uint8_t *curChoice = S->curChoice; // correspond au developpement en cours par  stage_dfs
    S->wAssign++;
    { // tueuse d'abord
        Perm *peKiller = &S->Pe[S->killer];
        for (int j = 0; j < nT; j++) {
            int dk = peKiller->choiceByP[j][curChoice[peKiller->recip[j]]];
            if (dk < curChoice[j]) return;    // killer est inferieure a curChoice, pas de changement
            if (dk > curChoice[j]) break;     // killer est superieur a curchoice => recalcul de killer
        }
    }
    for (int nPerm = 0; nPerm < S->nbPerm; nPerm++) {  // boucle sur toutes les permutations du stage
        int less = 0, greater = 0;
        Perm *pe = &S->Pe[nPerm];
        for (int j = 0; j < nT; j++) {
            int d = pe->choiceByP[j][curChoice[pe->recip[j]]];
            if (d < curChoice[j]) { less = 1; break; }
            if (d > curChoice[j]) { greater = 1; break; }
        }
        if (less) { S->killer = nPerm; return; } //perm lexicalement inferieur => devient killer
        if (!greater) stab++;           /* choice identique : stabilisateur */
    }
    /* canonique */
    if (S->nleaf == S->cap) {
        // allocation supplémentaire ?
        S->cap = S->cap ? S->cap * 2 : 1024;
        S->leafs = (Leaf *)realloc(S->leafs, (size_t)S->cap * sizeof(Leaf));
        if (!S->leafs) { fprintf(stderr, "memoire insuffisante\n"); exit(1); }
     }
    {
        Leaf *r = &S->leafs[S->nleaf++];
        memcpy(r->choice, curChoice, nT);
        r->weight = (short)(S->nbPerm / stab);
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
    for (uint8_t c = 0; c < 3; c++) {
        int p = tp[S->W[j]][c]; // paire
        if (S->cov[p] < 2) {
            S->cov[p]++; S->curChoice[j] = c;
            stage_dfs(S, j + 1);
            S->cov[p]--;
        }
    }
}
// *W : Liste de triplets
// n : nombre de triplets
// perms : liste de permutations [0,8[
// ng : nombre de permutations
// retourne la structure Stage correspondante
// applique au niveau 1 (35 triplets) et au niveau 2 (15 triplets)
static Stage *stage_build(Stage * S,int perms[][8]) {
    const int *W = S->W ;
    for (int nP = 0; nP < S->nbPerm; nP++) { // boucle sur les permutations
        for (int i = 0; i < S->nbT; i++) { // boucle sur les triplets
            int *t = tri[S->W[i]]; // on recupere le triplet
            int gi = tidx[perms[nP][t[0]]][perms[nP][t[1]]][perms[nP][t[2]]]; // on applique la permutation
            // et l'on recupere l'index correspondant du triplet associe
            int j ;
            for (j = 0; j < S->nbT; j++) if (S->W[j] == gi) break; // on recupere l'index du triplet associe
            if (j == S->nbT) { // ne doit jamais se produire
                fprintf(stderr, "W non invariant sous le groupe pour {%d,%d,%d}\n"
                    ,perms[nP][t[0]],perms[nP][t[1]],perms[nP][t[2]]);
                exit(1);
            }
            Perm *pe = &S->Pe[nP];
            pe->recip[j] = i; // i => j , on le note ds la correpondance inversée
            for (int8_t c = 0; c < 3; c++) { // boucle sur la paire associée
                int p = tp[W[i]][c]; // paire associee pour W[i]
                int q = pid_[perms[nP][P2[p][0]]][perms[nP][P2[p][1]]]; // q index de la paire en appliquant la permut g
                uint8_t k ;
                for (k = 0; k < 3; k++) if (tp[W[j]][k] == q) break;
                if (k == 3) {
                    // ne doit jamais se produire
                    fprintf(stderr, "paire image absente\n");
                    exit(1);
                }
                pe->choiceByP[j][c] = k; // on stocke le nouvel index (0, ou 1, ou 2)
            }
            // printf("Pe[%d](%d,%d,%d) \n",nP,pe->choiceByP[j][0],pe->choiceByP[j][1],pe->choiceByP[j][2]);
        }
    }
    return S;
}
void stage_phase1(Stage *S) {
    /* choiceByP[g][j][c] est indexe par le chiffre c du triplet ANTECEDENT : voir stage_leaf.
       Il faut aussi que choiceByP soit defini par l'antecedent i=inv[g][j] : c'est le cas ci-dessus. */
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

/* ------------------------------------------------- donnees globales */





static Stage * build_ST1(void) {
    static int perms1[NBPERM1][8] ; // liste des NBPERM permutation de [0..PERM1[
    static int  W1[MAXT] ; // liste des triplets de [0,PERM1[
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
    static int  W2[SIZE_W2] ; // liste des triplets de W2
    for (int k=0;k<SIZE_W2;k++) {
        W2[k]= tidx[liste_w2[k][0]][liste_w2[k][1]][liste_w2[k][2]];
    }
    static int perms2[NBPERM2][8];  // liste des 2 permutations 01234567 01234576
    Stage *ST2 = (Stage *)calloc(1, sizeof(Stage));
 //  int nbPerm2 = gen_group(5, 3, perms2); // permutations de {5,6,7}
 //   int nbPerm2 = gen_group(6, 2, perms2); // permutations de {5,6,7}
    int nbPerm2 = gen_group(PERM2,G_PERM2,perms2);
    if (nbPerm2 != NBPERM2) {
        fprintf(stderr,"ERROR nbPerm2=%d expexted %d\n",nbPerm2,NBPERM2);
        exit(1);
    }
    ST2->nbT =SIZE_W2; ST2->nbPerm = nbPerm2; ST2->W = W2;
    ST2 = stage_build(ST2, perms2);
    return ST2 ;
}




static void  build_remT(typeRemCov remT[][3],const Stage *ST1,const Stage *ST2) {
    int used[nbTriplets] = {0};
    for (int i = 0; i < MAXT; i++) used[ST1->W[i]] = 1;
    for (int i = 0; i < SIZE_W2; i++) used[ST2->W[i]] = 1;
    for (int j = 0; j < nbRemT; j++) {
        for (int i = 0; i < nbTriplets; i++){
            if (!used[i]) {
                used[i] = 1;
                for (uint8_t k = 0; k < 3; k++) {
                    remT[j][k] = tp[i][k];
                }
                break ;
            }
        }
    }
}

/* ------------------------------------------------ backtracking */

static typeRemCov  remT[nbRemT][3] ;
static typeRemCov  cov[nbPairs] ;
static typeRemCov  rem[nbPairs] ;

static u64 cnt ;


//static bt_dev BT = {{{0}},{0},{0},0};
// recursif, developpe les triplets complementaires (nbRemT) apres ceux developpes par les permutations
// va de niveau 0 a nbRemT (56-(20+6)
// calcule cnt le nombre de developpements possibles.
// Utilise dans la phase 2 (jobs)
static void bt(const int k)
{
    if (k == nbRemT) { cnt++; return ; }
    int const pa = remT[k][0]; int const pb = remT[k][1]; int const pc = remT[k][2];
    if (cov[pa] + rem[pa] >= 2 && cov[pb] + rem[pb] >= 2 && cov[pc] + rem[pc] >= 2) {
        rem[pa]--; rem[pb]--; rem[pc]--;

        if (cov[pa] < 2 && cov[pb] + rem[pb] >= 2 && cov[pc] + rem[pc] >= 2) { // pas encore de double attribution
            cov[pa]++; bt(k + 1); cov[pa]--;
        }
        if (cov[pb] < 2 && cov[pa] + rem[pa] >= 2 && cov[pc] + rem[pc] >= 2) { // pas encore de double attribution
            cov[pb]++; bt(k + 1); cov[pb]--;
        }
        if (cov[pc] < 2 && cov[pa] + rem[pa] >= 2 && cov[pb] + rem[pb] >= 2) { // pas encore de double attribution
            cov[pc]++; bt(k + 1); cov[pc]--;
        }
        rem[pa]++; rem[pb]++; rem[pc]++;
    }
}
/* un job = un representant de l'etage 1 ; retourne T = somme w2*count */
static u64 run_job(long id,const typeRemCov rem0[nbPairs],const Stage *ST1,const Stage *ST2)
{
    typeRemCov base[nbPairs];
    u64 T = 0;
    const Leaf *r1 = &ST1->leafs[id];
    const int * W1 = ST1->W;
    const int * W2 = ST2->W;
    memset(base, 0, sizeof(base));
    for (int i = 0; i < MAXT; i++) {
        base[tp[W1[i]][r1->choice[i]]]++;
    }
    for (int q = 0; q < ST2->nleaf; q++) {
        const Leaf *r2 = &ST2->leafs[q];
        int ok = 1;
        memcpy(cov, base, sizeof(base));
        for (int i = 0; i < SIZE_W2; i++) {
            cov[tp[W2[i]][r2->choice[i]]]++;
        }
        memcpy(rem, rem0, nbPairs *sizeof(typeRemCov));
        for (int p = 0; p < nbPairs; p++)
            if (cov[p] > 2 || cov[p] + rem[p] < 2) { ok = 0; break; }
        if (!ok) continue;
        cnt = 0;
        bt(0);
        T += (u64)r2->weight * cnt;
    }
    return T;
}
/* -------------------------------------------------------- main */
static void print_info(FILE * fout,const Stage *ST1,const Stage *ST2)
{
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
    char *done = (char *)calloc((size_t)nb_leaf, 1);
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
    // if (freq_ndo==0) freq_ndo= 1;
    int mod_ndo = 0;
    //for (long id = 0; id < nb_leaf; id ++) {
    for (long id = nb_leaf-1; id >= 0; id --) {
        if (done[id]) continue;
        u64 T = run_job(id,rem0,ST1,ST2);
        fprintf(fp, "%ld %d %llu\n", id, (int)ST1->leafs[id].weight, T);
 //       ndo++;
        mod_ndo++;
//        if (ndo % freq_ndo == 0) {
        if (mod_ndo == freq_ndo) {
            fflush(fp);
            printf("[%lds] %ld jobs faits dans cette session (dernier id %ld)\n",
                   (long)(time(NULL) - t0), ndo+mod_ndo, id);
            ndo += mod_ndo;
            mod_ndo = 0 ;
        }
    }
    fprintf(fp,"#[%lds] %ld jobs faits dans cette session %s (dernier id %ld)\n",
           (long)(time(NULL) - t0), ndo,CIBLE_NAME, nb_leaf);
    fclose(fp);
    check_sum(nb_leaf,fname ) ;
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

int main(int argc, char **argv){
    const clock_t c0 = clock();
    init_basic();
    printf("Execution %s\n",CIBLE_NAME);
    Stage * ST1 = build_ST1(); printf("buid_st1 OK\n");
    Stage * ST2 = build_ST2();  printf("buid_st2 OK\n");
    stage_phase1(ST2);  printf("ST2 phase 1 OK\n");
    stage_phase1(ST1);  printf("ST1 phase 1 OK\n");
    build_remT(remT,ST1,ST2);
    fprintf(stderr, "[orbites construites en %.1f s]\n", (double)(clock() - c0) / CLOCKS_PER_SEC);

    if (argc >= 2 && strcmp(argv[1], "--info") == 0) { print_info(stdout,ST1,ST2); return 0; }
    if (argc >= 3 && strcmp(argv[1], "--sum") == 0) { check_sum(ST1->nleaf,argv[2] ) ; return 0; }

    static typeRemCov rem0[nbPairs] = {0};
    for (int i = 0; i < nbRemT; i++) { rem0[remT[i][0]]++; rem0[remT[i][1]]++; rem0[remT[i][2]]++; }

    if (argc >= 3 && strcmp(argv[1], "--sample") == 0) { int K = atoi(argv[2]); compute_sample(K,rem0,ST1,ST2); return 0; }

    char fname[64];
    if ( argc > 1) { strcpy(fname, argv[1]); } else {
        sprintf(fname, "res_%s_%d.txt",CIBLE_NAME, 0);
    }
    compute_jobs(fname,rem0,ST1,ST2) ;
    return 0;
}