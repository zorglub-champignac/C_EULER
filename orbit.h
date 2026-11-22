//
// Created by Jeannot on 22/11/2026.
//

#ifndef C_EULER_ORBIT_H
#define C_EULER_ORBIT_H
typedef int8_t typeRemCov ;
typedef unsigned long long u64;


#define nbPairs 28
#define nbTriplets 56
// 0x2000000=33 554 432

#define MAX_HASH_ALLOC 0x2000000
#define MAX_CONFLITS 60000000

typedef typeRemCov remT_Type[nbRemT][3];
typedef int64_t hashCov ;

/************************* structures **********************************/
typedef struct {
    typeRemCov pairs[MAXT] ;
    int16_t weight;
} Leaf;

typedef struct {
    typeRemCov cov[nbPairs];
    int64_t weight ;
} HashLeaf;

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

typedef struct {
    hashCov key ;
    u64 weight ;
} k_element ;
typedef struct {
    k_element *tb ;
    int nbKey ;
} key_tb ;
typedef struct {
    const Stage *ST1 ;
    key_tb ktb1 ;
    const Stage *ST2 ;
    key_tb ktb2 ;
    key_tb ketb ;
} all_key_tb ;


static void init_basic(void) ;
static Stage * build_ST1(void) ;
static Stage * build_ST2(void) ;
static Stage *stage_build(Stage * S,int perms[][8]) ;
static void stage_phase1(Stage *S) ;
remT_Type * build_remT(const Stage *ST1,const Stage *ST2) ;
static const typeRemCov *   build_rem0(const remT_Type * remT) ;

static void check_sum(long nb_key,char *name_fp ) ;
static void print_info(FILE * fout,all_key_tb all_tb) ;
static void compute_sample( int K,const typeRemCov *rem0,const Stage *ST1,const Stage *ST2) ;

static int gen_group(int lo, int m, int perms[][8]) ;
static all_key_tb compute_merge_dev(const typeRemCov *rem0,const Stage *ST1,const Stage *ST2) ;
static void compute_jobs( const char *fileName,const typeRemCov *rem0const,all_key_tb all_tb) ;

static void stage_dfs(Stage *S, int j ) ;
static void stage_leaf(Stage *S) ;


#endif //C_EULER_ORBIT_H
