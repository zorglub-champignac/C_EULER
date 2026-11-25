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
 *        - doit etre stable par la permutation 6<->7)
 *   - pour les 2 levels on precalcule pour chaque permutation :
 *        - recip_tri[num Triplet] : la fonction reciproque sur les triplets
 *           recip[tidx[perm(a)][perm(b)][perm(c)] = tidx[a][b][c]
 *            Cela est possible pour ST2 car les triplets de la liste W2 doivent etre stable par le groupe des permutations
 *            (en fait groupe de taille 2 : Identite+6<->7
 *         - direct_pair[num Pair] :
 *           direct_pair[pidx[a][b]] = pidx[perm(a)][perm(b)]
 *         Exemple pour ST1 et la permutation (0->1->2->4->0).
 *             Le triplet (0,1,2) devient (1,2,4)
 *             - (0,1,2) (n° 0) devient (1,2,4) (n° 22) : donc recip[22] = 0
 *          Meme methode que orbit_count.c, avec des groupes plus gros :
 *            Etage 1 : W1 = les 20 triplets de {0,...,5}, groupe G1 = S6 (720 perms,
 *           6 et 7 fixes).
 *            Etage 2 : W2 = les 15 triplets contenant un 6, un 7 où les deux.
 *                G2 fixe point par point toute affectation de W1.
 *            N = somme sur les representants r1 (orbites de G1 sur les affectations
 *                  faisables de W1) de  |O1| *
 *                somme sur les representants r2 (orbites de G2 sur les affectations
 *                  faisables de W2) de  |O2| * count(r1, r2)
 *            Facteur nominal 720 * 2 = 1440. Le comptage est exact (orbites courtes
 *                  ponderees par leur vraie taille).
 *             Canonicalisation : le representant d'une orbite est l'affectation
 *                lexicographiquement minimale (chiffres 0,1,2 dans l'ordre des triplets de W).
 *                Le test est paresseux : comparaison chiffre par chiffre avec sortie
 *                  anticipee, et la derniere permutation "tueuse" est essayee en premier.
 *    - regroupement des développements.
*          L'astuce reside dans le fait que le calcul par lie au developpement d'une orbite juesqu'à la fin
*          depend seulement des paires consommées cov[i], i dans[0 28[ durant la phase de développement par triplets.
*          On peut donc économiser des developpements conduisant au meme valeurs de cov[i] en simplement ajoutant les poids.
*          Cette phase de regroupement est effectue par compute_merge_dev().
*          A noter que l'on peut suivant le meme procédé regrouper auparavant les developpements de ST1 et les developpements de ST2.
*          Pour effectuer ces regroupements on utilise une clé calculée sur cov[i], i[0,28[
*          en remarquant que cov[i] peut prendre les valeurs 0,1,2 . Base 3 sur 28 valeurs 3**28 < 2**64 donc la clé tient sur u64.
 *
 * UTILISATION
 *   ./orbit_general --info           structure (construit les orbites)
 *   ./orbit_general                         tous les jobs (fichier res6v2_0.txt)
 *   ./orbit_general --sum res6v2_0.txt ... total et jobs manquants
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


#define CIBLE 20
#include "cible.h"

#include "orbit.h"



/* ---------------------------------------------------------------- base */
static typeRemCov pidx[8][8]; // tableau des index[0,28[ des paires ((i,j) et (j,i) => index )
static typeRemCov P2[nbPairs][2]; // liste des paires
static typeRemCov tri[nbTriplets][3]; // liste des triplets
static typeRemCov tidx[8][8][8];//  tableau des index [0,56[ des triplets ordonnés (a<b<c => index)
static typeRemCov tp[nbTriplets][3]; // liste des index de paires associées au triplet pid[a][b], pid[a][c], pid[b][c]
/* ------------------------------------------------ backtracking */

#define MEMO_BITS 20
#define MEMO_LEVEL 12
static u64 hashKey(hashCov key) {
    __uint128_t hash = ( (((__uint128_t)key) * 0x9E3779B97F4A7C15)& 0xFFFFFFFFFFFFFFFF ) >> (64 - MEMO_BITS)  ;
    return (hashCov) hash ;
}

typedef struct {
    hashCov key ;
    u64 cnt ;

} hash_element;

remT_Type * get_remT(remT_Type * remT) {
    static remT_Type * initRemT = NULL ;
    if ( remT == NULL ) {
        return initRemT;
    } else {
        initRemT = remT ;
        return remT ;
    }
}

hash_element *tbHashElement = NULL ;
static void bt(int k) ;

static void bt_memo(int k,int32_t hash_key );

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
 //   if (argc >= 3 && strcmp(argv[1], "--sample") == 0) { int K = atoi(argv[2]); compute_sample(K,rem0,ST1,ST2); return 0; }
    all_key_tb all_tb = compute_merge_dev(rem0,ST1,ST2);
    if (argc >= 2 && strcmp(argv[1], "--info") == 0) {
        print_info(stdout,all_tb); return 0;
    }
    if (argc >= 3 && strcmp(argv[1], "--sum") == 0) {
    //    key_tb ktb = compute_merge_dev(rem0,ST1,ST2);
        check_sum(all_tb.ketb.nbKey,argv[2] ) ;
        return 0;
    }
    char fname[64];
    if ( argc > 1) { strcpy(fname, argv[1]); } else {
        sprintf(fname, "C%d_S%dxS%d_tri%dX%d.txt",CIBLE,PERM1,G_PERM2,MAXT,SIZE_W2);
    }
    // *********  phase 2 **************************
    // developpement des feuilles par developpement recursif sur les paires de bt()
    // suivant ls
    compute_jobs(fname,rem0,all_tb) ;
    return 0;
}
/************ BEGIN ***************** phase 1 pour ST2 et ST1 *****************************/
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
// donnees partagées par bt()
static typeRemCov  cov[nbPairs] ;
static typeRemCov rem_cov[nbPairs] ;

int cmp_k_element (const void * k1, const void * k2) {
    hashCov diff = ((const k_element *)k1)->key - ((const k_element *)k2)->key ;
    if (diff > 0) return 1 ;
    else if (diff < 0) return -1 ;
    else return 0 ;
}
void merge_k_element(k_element *k1,const k_element *k2) {
    k1->weight += k2->weight ;
}

int compute_key_tb_diff (key_tb * tb_elem) {
    qsort(tb_elem->tb,tb_elem->nbKey,sizeof(k_element),cmp_k_element) ;
    // hashCov ant_key = tb_elem->tb[0].key ;
    int nbZero =0 ;
    for (nbZero=0;nbZero<tb_elem->nbKey;nbZero++) {
        if (tb_elem->tb[nbZero].key) break ;
    }
    int nbElemDif = 1 ;
    hashCov ant_key = tb_elem->tb[nbZero].key ;
    tb_elem->tb[0] = tb_elem->tb[nbZero] ;
    for (int i=nbZero+1;i< tb_elem->nbKey ;i++) {
        if ( tb_elem->tb[i].key != ant_key)  {
            if (nbElemDif != i) {
               tb_elem->tb[nbElemDif++] = tb_elem->tb[i] ;
            } else {
                nbElemDif++ ;
            }
            ant_key = tb_elem->tb[i].key ;
        } else {
            merge_k_element(tb_elem->tb+nbElemDif-1,tb_elem->tb+i) ;
        }
    }
    // printf("nbKey=%d Diff=%d\n",tb_elem->nbKey,nbElemDif) ;
    tb_elem->nbKey = nbElemDif ;
    return nbElemDif ;
}

static u64 cnt ;
// static u64 cnt_cache ;
// static uint32_t has_key = 0 ;
key_tb compact_key(const Stage *ST);

/*
hashCov covToKey(const typeRemCov *cov_) {
    hashCov hC = 0 ;
    for (int i=0;i<nbPairs;i++) {hC += pow3[i] * cov_[i] ;}
    return hC ;
}
int keyToCov(hashCov key,typeRemCov *cov_) {
    memset(cov_,0,nbPairs);
    for (int i=nbPairs-1; i>=0; i--) {
        if (key >= pow3[i]) {
            cov_[i]++ ; key -= pow3[i] ;
            if (key >= pow3[i]) {
                cov_[i]++ ; key -= pow3[i] ;
            }
        }
    }
    return (key ? 1:0) ; // doit etre null
}
*/
hashCov covToKey(const typeRemCov *cov_) {
    hashCov hC = 0 ;
    for (int i=0;i<nbPairs;i++) {hC +=  ((u64)cov_[i] & 3) << (2*i) ;}
    return hC ;
}
void keyToCov(hashCov key,typeRemCov *cov_) {
    memset(cov_,0,nbPairs);
    for (int i=0;i<nbPairs;i++) {
        cov_[i] = (((u64)key) >> (2*i)) & 3 ;
    }
//    return (key ? 1:0) ; // doit etre null
}

time_t t0 ;
static all_key_tb compute_merge_dev(const typeRemCov *rem0,const Stage *ST1,const Stage *ST2) {
    t0 = time(NULL);
    const key_tb ketb1 = compact_key(ST1) ;
    const key_tb ketb2 = compact_key(ST2) ;
    all_key_tb all_key ;
    all_key.ST1 = ST1 ;
    all_key.ST2 = ST2 ;
    all_key.ktb1 = ketb1 ;
    all_key.ktb2 = ketb2 ;
    int64_t nbKeyAlloc = (int64_t)ketb1.nbKey * ketb2.nbKey ;
    int overbook = 0 ;
    int indexConflicts = 0 ;
    if ( nbKeyAlloc > MAX_HASH_ALLOC ) {
        nbKeyAlloc = MAX_HASH_ALLOC ;
        indexConflicts = nbKeyAlloc ;
        overbook = 1 ;
        nbKeyAlloc +=  MAX_CONFLITS ; // rajout de la zone pour gerer les conflits
    }
    all_key.ketb.tb = calloc(nbKeyAlloc,sizeof(all_key.ketb.tb[0]));
    int nb_conflits = 0 ;
    for (int i1=0; i1<ketb1.nbKey; i1++) {
        for (int i2=0; i2<ketb2.nbKey; i2++) {
            typeRemCov cov_[nbPairs] ;
            hashCov hC = ketb1.tb[i1].key + ketb2.tb[i2].key ;
            keyToCov(hC,cov_) ;
            int ok = 1;
            for (int p=0;p<nbPairs;p++) {
                if (cov_[p] > 2 || cov_[p] + rem0[p] < 2 ) { ok= 0; break ; }
                rem_cov[p] = (typeRemCov)(cov_[p] + rem0[p]) ;
            }
            if (!ok) continue ;
            int index = i1*ketb2.nbKey+i2 ;
            if (overbook ) {
                index &= MAX_HASH_ALLOC-1 ;
                if  (all_key.ketb.tb[index].key !=hC &&  all_key.ketb.tb[index].key != 0) {
                    index = indexConflicts + nb_conflits++ ;
                    if (nb_conflits >= MAX_CONFLITS) {
                        fprintf(stderr,"Limite de %d conflits atteinte. %d/%ld id faits \n "
                            ,MAX_CONFLITS,i1,ST1->nleaf) ;
                        exit(5);
                    }
                }
                // else {
                //    if  (all_key.ketb.tb[index].key !=hC  && all_key.ketb.tb[index].key != 0 ) {
                 //       nb_conflits++ ;
                 //   }
                // }
            }
            all_key.ketb.tb[index].key = hC ;
            all_key.ketb.tb[index].weight += ketb1.tb[i1].weight * ketb2.tb[i2].weight ;
        }
    }
    all_key.ketb.nbKey = nbKeyAlloc ;
    compute_key_tb_diff(&all_key.ketb) ;
    all_key.ketb.tb = realloc(all_key.ketb.tb,all_key.ketb.nbKey*sizeof(all_key.ketb.tb[0])) ;
    printf("%lds Total reduction %d->%d \n",(long)(time(NULL) - t0),ketb1.nbKey*ketb2.nbKey,all_key.ketb.nbKey);
    HashLeaf * tbKeyLeaf = calloc (all_key.ketb.nbKey, sizeof(HashLeaf));
    int nbKeyLeaf = 0;
    for (int id=0;id<all_key.ketb.nbKey;id++) {
        if (all_key.ketb.tb[id].key== 0 ) continue;
        keyToCov(all_key.ketb.tb[id].key,tbKeyLeaf[nbKeyLeaf].cov);
        tbKeyLeaf[nbKeyLeaf].weight = all_key.ketb.tb[id].weight;
        nbKeyLeaf++ ;
    }
    fflush(stdout);
    return all_key ;
}

key_tb compact_key(const Stage *ST) {
    key_tb ketb ;
    ketb.tb =  malloc(ST->nleaf*sizeof(ketb.tb[0])); ;
    typeRemCov cov1[nbPairs];
    for (int id = 0; id < ST->nleaf; id++) {
        const Leaf *r = &ST->leafs[id];
        memset(cov1, 0, sizeof(cov1));
        for (int i = 0; i < ST->nbT; i++) {  // contrib stage 1
            cov1[r->pairs[i]]++ ;
        }
        hashCov hC = 0 ;
        hC = covToKey(cov1);
        /*
        for (int p = 0; p < nbPairs; p++) {
            hC += pow3[p] * cov1[p] ;
        }
        */
        ketb.tb[id].key= hC ;
        ketb.tb[id].weight= r->weight ;
    }
    ketb.nbKey = ST->nleaf ;
    ketb.nbKey =  compute_key_tb_diff(&ketb) ;
    ketb.tb = realloc(ketb.tb,ketb.nbKey*sizeof(ketb.tb[0])); ;
    // printf("Size ketb=%d\n",ketb.nbKey) ;
    return ketb ;
}

void compute_jobs( const char *fileName,const typeRemCov *rem0,all_key_tb all_tb)   {
    // long nb_leaf = ST1->nleaf;
    char fname[64];
    // t0 = time(NULL);
    if (fileName != NULL) {
        strcpy(fname, fileName);
    } else {
        sprintf(fname, "resv2_%s_%d.txt",CIBLE_NAME,0);
    }
    FILE *fp = fopen(fname, "r");
    int begin_jid = 0 ;
    u64 Total = 0 ;
    if (fp) {
        long jid;
        char line[200] ;
        while (fgets(line,200,fp)) {
            if (line[0] != '#') {
                if (sscanf(line, "%ld %llu", &jid, &Total) == 2) {
                    if (jid >= begin_jid  && jid <= all_tb.ketb.nbKey) {
                        begin_jid = jid ;
                    }
                }
            }
        }
        fclose(fp);
    }
    print_info(stdout,all_tb);
    fp = fopen(fname, "a");
    if (!fp) { perror("res"); return ; }
    print_info(fp,all_tb);
    printf(" fichier %s, deja faits : %d\n", fname, begin_jid);

    tbHashElement = calloc ((1 << MEMO_BITS) , sizeof(tbHashElement[0]));
    int nh_period = all_tb.ketb.nbKey / 100 ;
    nh_period -= (nh_period % 10000) ;
    int nh= 0 ;
    int np = 0 ;
    if (begin_jid) {
        nh = begin_jid+1 ;
        np=1 ;
    }
    for (; nh < all_tb.ketb.nbKey ; np++, nh++) {
        hashCov key = all_tb.ketb.tb[nh].key ;
        keyToCov(key,cov);
        for (typeRemCov p =0; p< nbPairs;  p++) {
            rem_cov[p] = (typeRemCov)(cov[p] + rem0[p]) ;
        }
        cnt = 0;

        bt(0);
        u64 T = cnt ;

        Total +=  all_tb.ketb.tb[nh].weight * T;
        if (np == nh_period ) {
            printf("%lds %d->%lld\n",(long)(time(NULL) - t0),nh,Total); fflush(stdout);
            fprintf(fp, "%d %llu\n", nh,Total);
            np = 0 ;
        }
    }
    if (nh <= all_tb.ketb.nbKey) {
        fprintf(fp, "%d %llu\n", nh,Total);
        printf("%lds %d=> T=%lld \n",(long)(time(NULL) - t0),nh,Total);
        fprintf(fp,"#[%lds] %d jobs faits dans cette session %s (dernier id %d) Total=%lld\n",
           (long)(time(NULL) - t0), all_tb.ketb.nbKey,CIBLE_NAME,all_tb.ketb.nbKey,Total);
    }

    fclose(fp);
    check_sum(all_tb.ketb.nbKey,fname ) ;
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
    if (k == nbRemT) {
        if (MEMO_LEVEL < nbRemT) {
            printf("Error for level %d over %s must be bt_memo()\n") ;
            exit (4);
        }
        cnt++;
        return ;
    }
    const typeRemCov pa = remT_bt[k][0]; const typeRemCov  pb = remT_bt[k][1]; const typeRemCov  pc = remT_bt[k][2];
    if (rem_cov[pa] >= 2 && rem_cov[pb] >= 2 &&  rem_cov[pc] >= 2) {
        rem_cov[pa]-- ; rem_cov[pb]-- ; rem_cov[pc]-- ;
        if (cov[pa] < 2 && rem_cov[pb] >= 2 && rem_cov[pc] >= 2) { // pas encore de double attribution
                cov[pa]++;rem_cov[pa]++;
                if (k== MEMO_LEVEL-1) {
                    bt_memo(k + 1,0);
                } else {
                    bt(k + 1);
                }
                cov[pa]--;rem_cov[pa]--;
        }
         if (cov[pb] < 2 && rem_cov[pa] >= 2 && rem_cov[pc] >= 2) { // pas encore de double attribution
            cov[pb]++;rem_cov[pb]++;
             if (k== MEMO_LEVEL-1) {
                 bt_memo(k + 1,0);
             } else {
                 bt(k + 1);
             }
             cov[pb]--;rem_cov[pb]--;
        }
        if (cov[pc] < 2 && rem_cov[pa] >= 2 && rem_cov[pb] >= 2) { // pas encore de double attribution
            cov[pc]++;rem_cov[pc]++;
            if (k== MEMO_LEVEL-1) {
                bt_memo(k + 1,0);
            } else {
                bt(k + 1);
            }
            cov[pc]--;rem_cov[pc]--;
        }
        rem_cov[pa]++ ; rem_cov[pb]++ ; rem_cov[pc]++ ;
    }
}


void bt_memo(int k,int hash_key) {
    if ( k== MEMO_LEVEL) {
        hashCov key = covToKey(cov);
        hash_key = hashKey(key) ;
        if (tbHashElement[hash_key].key == key) {
            cnt += tbHashElement[hash_key].cnt ;
            return ;
        } else {
            tbHashElement[hash_key].key = key ;
            tbHashElement[hash_key].cnt = 0 ;
            // cnt_cache = 0 ;
        }
    } else if ( k == nbRemT) {
        tbHashElement[hash_key].cnt ++ ;
        cnt ++ ;
        return ;
    }
    const typeRemCov pa = remT_bt[k][0]; const typeRemCov  pb = remT_bt[k][1]; const typeRemCov  pc = remT_bt[k][2];
    if (rem_cov[pa] >= 2 && rem_cov[pb] >= 2 &&  rem_cov[pc] >= 2) {
        rem_cov[pa]-- ; rem_cov[pb]-- ; rem_cov[pc]-- ;
        if (cov[pa] < 2 && rem_cov[pb] >= 2 && rem_cov[pc] >= 2) { // pas encore de double attribution
            cov[pa]++;rem_cov[pa]++; bt_memo(k + 1,hash_key); cov[pa]--;rem_cov[pa]--;
        }
        if (cov[pb] < 2 && rem_cov[pa] >= 2 && rem_cov[pc] >= 2) { // pas encore de double attribution
            cov[pb]++;rem_cov[pb]++; bt_memo(k + 1,hash_key); cov[pb]--;rem_cov[pb]--;
        }
        if (cov[pc] < 2 && rem_cov[pa] >= 2 && rem_cov[pb] >= 2) { // pas encore de double attribution
            cov[pc]++;rem_cov[pc]++; bt_memo(k + 1,hash_key); cov[pc]--;rem_cov[pc]--;
        }
        rem_cov[pa]++ ; rem_cov[pb]++ ; rem_cov[pc]++ ;
    }

}
/***** END ************************* phase 2 (jobs) *****************************/

/*************************************** routines tritement options de main() **************/

static void print_info(FILE * fout,all_key_tb all_tb) {
    long  tw1 = 0, tw2 = 0, full = 0;
    for (long i = 0; i < all_tb.ST1->nleaf; i++) {
        tw1 += all_tb.ST1->leafs[i].weight;
        if (all_tb.ST1->leafs[i].weight == all_tb.ST1->nbPerm) full++;
    }
    for (long i = 0; i < all_tb.ST2->nleaf; i++) tw2 += all_tb.ST2->leafs[i].weight;
    fprintf(fout,"# Etage 1 : %ld affectations faisables de W1 (%d triplets), %ld orbites sous S6 "
           "(somme poids = %ld, orbites pleines de %d : %ld)\n"
           , all_tb.ST1->wAssign,all_tb.ST1->nbT, all_tb.ST1->nleaf, tw1,all_tb.ST1->nbPerm, full);
    fprintf(fout,"# Etage 2 : %ld affectations faisables de W2 (%d triplets), %ld orbites sous S2 "
           "(somme poids = %ld)\n", all_tb.ST2->wAssign,all_tb.ST2->nbT, all_tb.ST2->nleaf, tw2);
    fprintf(fout,"# Triplets restants a parcourir : %d\n", nbRemT);
    fprintf(fout,"# Etage 1 %ld->%d ,  Etage 2 %ld->%d  Total=%dx%d=%lld MAX=%d\n"
            ,all_tb.ST1->nleaf,all_tb.ktb1.nbKey,all_tb.ST2->nleaf,all_tb.ktb2.nbKey
            ,all_tb.ktb1.nbKey,all_tb.ktb2.nbKey,(int64_t)all_tb.ktb1.nbKey * all_tb.ktb2.nbKey,MAX_HASH_ALLOC) ;

    fprintf(fout,"# Nombre de developpement (apres le merge) = %ld \n",
       all_tb.ketb.nbKey);
    double wAssign = ((double)all_tb.ST1->wAssign)*all_tb.ST2->wAssign ;
    fprintf(fout,"# Gain nominal sur les %.0f configurations de depart : %.1f\n",
           wAssign,(double)wAssign / ((double)all_tb.ketb.nbKey));

    fflush(fout);
}

void check_sum(long nb_key,char *name_fp ) {
    unsigned long long total = 0;
    long  missing = 0;
    printf("CHECKSUM nb_key=%d\n", nb_key);
    FILE *fp = fopen(name_fp, "r");
    long jid;  u64 T;
    if (!fp) {
        fprintf(stderr, "ouverture impossible : %s\n", name_fp);
    }else {
        char line[200] ;
        while (fgets(line,200,fp)) {
            if (line[0] != '#') {
                if (sscanf(line, "%ld %llu", &jid, &T) == 2) {
                    if (jid >= 0 && jid <= nb_key ) {
                        total =  T;
                    }
                }
            }
        }
    }
    fclose(fp);

    missing = nb_key -jid  ;
    printf("jobs termines : %ld / %ld   manquants : %ld\n", jid, nb_key, missing);
    printf("TOTAL (partiel si manquants) = %llu\n", total);
    if (missing == 0)
        printf("Comparaison avec N attendu 125047355835960 : %s\n",
               total == 125047355835960ULL ? "IDENTIQUE" : "DIFFERENT");

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
static typeRemCov l_rem0[nbPairs] = {0};
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
        for (typeRemCov j = (typeRemCov)(i + 1); j < 8; j++) {
            pidx[i][j] = pidx[j][i] = np;
            P2[np][0] = i; P2[np][1] = j; np++;
        }
    for (int a = 0; a < 8; a++)
        for (int b = a + 1; b < 8; b++)
            for (int c = b + 1; c < 8; c++) {
                tri[nt][0] = (typeRemCov)a; tri[nt][1] = (typeRemCov)b; tri[nt][2] = (typeRemCov)c;
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
