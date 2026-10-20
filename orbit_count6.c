/*
 * orbit_count6.c
 *
 * Meme methode que orbit_count.c, avec des groupes plus gros :
 *
 * Etage 1 : W1 = les 20 triplets de {0,...,5}, groupe G1 = S6 (720 perms,
 *           6 et 7 fixes).
 * Etage 2 : W2 = les 6 triplets {a,6,7}, a=0..5, groupe G2 = S2 = {id,(6 7)}.
 *           G2 fixe point par point toute affectation de W1.
 *
 *   N = somme sur les representants r1 (orbites de G1 sur les affectations
 *       faisables de W1) de  |O1| *
 *       somme sur les representants r2 (orbites de G2 sur les affectations
 *       faisables de W2) de  |O2| * count(r1, r2)
 *
 * Facteur nominal 720 * 2 = 1440. Le comptage est exact (orbites courtes
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
 *   ./orbit_count6                         tous les jobs (fichier res6_0.txt)
 *   ./orbit_count6 PART NPARTS             jobs d'id % NPARTS == PART
 *                                          (fichier res6_PART.txt), reprise auto
 *   ./orbit_count6 --sum res6_0.txt ...    total et jobs manquants
 *
 * Un job = un representant r1 de l'etage 1 ; il boucle sur les representants
 * de l'etage 2 et ecrit une ligne "id poids1 T" avec T = somme poids2*count.
 * Resultat attendu : N = 125 047 355 835 960.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef unsigned long long u64;

/* ---------------------------------------------------------------- base */
static int pid_[8][8];
static int P2[28][2];
static int tri[56][3];
static int tidx[8][8][8];
static int tp[56][3];

static void init_basic(void)
{
    int np = 0, nt = 0, i, j, a, b, c;
    for (i = 0; i < 8; i++)
        for (j = i + 1; j < 8; j++) {
            pid_[i][j] = pid_[j][i] = np;
            P2[np][0] = i; P2[np][1] = j; np++;
        }
    for (a = 0; a < 8; a++)
        for (b = a + 1; b < 8; b++)
            for (c = b + 1; c < 8; c++) {
                tri[nt][0] = a; tri[nt][1] = b; tri[nt][2] = c;
                tidx[a][b][c] = nt;
                tp[nt][0] = pid_[a][b];
                tp[nt][1] = pid_[a][c];
                tp[nt][2] = pid_[b][c];
                nt++;
            }
}

static int tindex3(int a, int b, int c)
{
    int t;
    if (a > b) { t = a; a = b; b = t; }
    if (b > c) { t = b; b = c; c = t; }
    if (a > b) { t = a; a = b; b = t; }
    return tidx[a][b][c];
}

static int next_perm(int *a, int n)
{
    int i = n - 2, j, t;
    while (i >= 0 && a[i] > a[i + 1]) i--;
    if (i < 0) return 0;
    j = n - 1;
    while (a[j] < a[i]) j--;
    t = a[i]; a[i] = a[j]; a[j] = t;
    for (j = i + 1, i = n - 1; j < i; j++, i--) { t = a[j]; a[j] = a[i]; a[i] = t; }
    return 1;
}

static int gen_group(int lo, int m, int perms[][8])
{
    int a[8], ng = 0, i;
    for (i = 0; i < m; i++) a[i] = lo + i;
    do {
        for (i = 0; i < 8; i++) perms[ng][i] = i;
        for (i = 0; i < m; i++) perms[ng][lo + i] = a[i];
        ng++;
    } while (next_perm(a, m));
    return ng;
}

/* ------------------------------------------------------- etages */
#define MAXN 20
typedef struct {
    unsigned char ch[MAXN];
    short weight;
} Rep;

typedef struct {
    int n, ng;
    const int *W;
    int inv[720][MAXN];                    /* inv[g][j] = i tel que g envoie i en j */
    unsigned char dig[720][MAXN][3];       /* dig[g][j][c] = chiffre image en j si le chiffre en inv vaut c */
    int cov[28];
    unsigned char ch[MAXN];
    Rep *reps; long nreps, cap;
    long nassign;
    int killer;
} Stage;

static void stage_leaf(Stage *S)
{
    int g, j, stab = 0, n = S->n;
    const unsigned char *ch = S->ch;
    S->nassign++;
    /* tueuse d'abord */
    {
        g = S->killer;
        for (j = 0; j < n; j++) {
            int d = S->dig[g][j][ch[S->inv[g][j]]];
            if (d < ch[j]) return;
            if (d > ch[j]) break;
        }
    }
    for (g = 0; g < S->ng; g++) {
        int less = 0, greater = 0;
        for (j = 0; j < n; j++) {
            int d = S->dig[g][j][ch[S->inv[g][j]]];
            if (d < ch[j]) { less = 1; break; }
            if (d > ch[j]) { greater = 1; break; }
        }
        if (less) { S->killer = g; return; }
        if (!greater) stab++;           /* image identique : stabilisateur */
    }
    /* canonique */
    if (S->nreps == S->cap) {
        S->cap = S->cap ? S->cap * 2 : 1024;
        S->reps = (Rep *)realloc(S->reps, (size_t)S->cap * sizeof(Rep));
        if (!S->reps) { fprintf(stderr, "memoire insuffisante\n"); exit(1); }
    }
    {
        Rep *r = &S->reps[S->nreps++];
        memcpy(r->ch, ch, n);
        r->weight = (short)(S->ng / stab);
    }
}

static void stage_dfs(Stage *S, int i)
{
    int c, p;
    if (i == S->n) { stage_leaf(S); return; }
    for (c = 0; c < 3; c++) {
        p = tp[S->W[i]][c];
        if (S->cov[p] < 2) {
            S->cov[p]++; S->ch[i] = (unsigned char)c;
            stage_dfs(S, i + 1);
            S->cov[p]--;
        }
    }
}

static Stage *stage_build(const int *W, int n, int perms[][8], int ng)
{
    Stage *S = (Stage *)calloc(1, sizeof(Stage));
    int g, i, j, c, k;
    S->n = n; S->ng = ng; S->W = W; S->killer = 0;
    for (g = 0; g < ng; g++)
        for (i = 0; i < n; i++) {
            int *t = tri[W[i]];
            int gi = tindex3(perms[g][t[0]], perms[g][t[1]], perms[g][t[2]]);
            for (j = 0; j < n; j++) if (W[j] == gi) break;
            if (j == n) { fprintf(stderr, "W non invariant sous le groupe\n"); exit(1); }
            S->inv[g][j] = i;
            for (c = 0; c < 3; c++) {
                int p = tp[W[i]][c];
                int q = pid_[perms[g][P2[p][0]]][perms[g][P2[p][1]]];
                for (k = 0; k < 3; k++) if (tp[W[j]][k] == q) break;
                if (k == 3) { fprintf(stderr, "paire image absente\n"); exit(1); }
                S->dig[g][j][c] = (unsigned char)k;
            }
        }
    /* dig[g][j][c] est indexe par le chiffre c du triplet ANTECEDENT : voir stage_leaf.
       Il faut aussi que dig soit defini par l'antecedent i=inv[g][j] : c'est le cas ci-dessus. */
    stage_dfs(S, 0);
    {
        long sw = 0; long q;
        for (q = 0; q < S->nreps; q++) sw += S->reps[q].weight;
        if (sw != S->nassign) {
            fprintf(stderr, "ERREUR: somme des poids %ld != %ld affectations\n", sw, S->nassign);
            exit(1);
        }
    }
    return S;
}

/* ------------------------------------------------- donnees globales */
static int  W1[20], W2[6];
static Stage *ST1, *ST2;

static int  R[48][3];
static int  NR;

static void build_stages(void)
{
    static int perms1[720][8], perms2[2][8];
    int ng1, ng2, a, b, c, k = 0;
    for (a = 0; a < 6; a++)
        for (b = a + 1; b < 6; b++)
            for (c = b + 1; c < 6; c++) W1[k++] = tidx[a][b][c];
    for (a = 0; a < 6; a++) W2[a] = tidx[a][6][7];
    ng1 = gen_group(0, 6, perms1);
    ng2 = gen_group(6, 2, perms2);
    ST2 = stage_build(W2, 6, perms2, ng2);
    ST1 = stage_build(W1, 20, perms1, ng1);
}

static void build_order(void)
{
    int used[56], i, j, k, p, best, bestscore;
    int remtot[28], remleft[28];
    memset(used, 0, sizeof(used));
    for (i = 0; i < 20; i++) used[W1[i]] = 1;
    for (i = 0; i < 6; i++) used[W2[i]] = 1;
    memset(remtot, 0, sizeof(remtot));
    for (i = 0; i < 56; i++) if (!used[i])
        for (k = 0; k < 3; k++) remtot[tp[i][k]]++;
    memcpy(remleft, remtot, sizeof(remtot));
    NR = 56 - 26;
    for (j = 0; j < NR; j++) {
        best = -1; bestscore = -1000000;
        for (i = 0; i < 56; i++) if (!used[i]) {
            int sc = 0;
            for (k = 0; k < 3; k++) {
                p = tp[i][k];
                if (remleft[p] == 1) sc += 10;
                if (remleft[p] == remtot[p]) sc -= 1;
            }
            if (sc > bestscore) { bestscore = sc; best = i; }
        }
        used[best] = 1;
        for (k = 0; k < 3; k++) { R[j][k] = tp[best][k]; remleft[tp[best][k]]--; }
    }
}

/* ------------------------------------------------ backtracking */
static int  cov[28], rem[28];
static u64  cnt;

static void bt(int k)
{
    int a, b, c, j, p;
    if (k == NR) { cnt++; return; }
    a = R[k][0]; b = R[k][1]; c = R[k][2];
    rem[a]--; rem[b]--; rem[c]--;
    for (j = 0; j < 3; j++) {
        p = R[k][j];
        if (cov[p] < 2) {
            cov[p]++;
            if (cov[a] + rem[a] >= 2 && cov[b] + rem[b] >= 2 && cov[c] + rem[c] >= 2)
                bt(k + 1);
            cov[p]--;
        }
    }
    rem[a]++; rem[b]++; rem[c]++;
}

static int rem0[28];

/* un job = un representant de l'etage 1 ; retourne T = somme w2*count */
static u64 run_job(long id)
{
    int base[28], i, p, q;
    u64 T = 0;
    const Rep *r1 = &ST1->reps[id];
    memset(base, 0, sizeof(base));
    for (i = 0; i < 20; i++) base[tp[W1[i]][r1->ch[i]]]++;
    for (q = 0; q < ST2->nreps; q++) {
        const Rep *r2 = &ST2->reps[q];
        int ok = 1;
        memcpy(cov, base, sizeof(base));
        for (i = 0; i < 6; i++) cov[tp[W2[i]][r2->ch[i]]]++;
        memcpy(rem, rem0, sizeof(rem0));
        for (p = 0; p < 28; p++)
            if (cov[p] > 2 || cov[p] + rem[p] < 2) { ok = 0; break; }
        if (!ok) continue;
        cnt = 0;
        bt(0);
        T += (u64)r2->weight * cnt;
    }
    return T;
}

static void init_rem0(void)
{
    int i;
    memset(rem0, 0, sizeof(rem0));
    for (i = 0; i < NR; i++) { rem0[R[i][0]]++; rem0[R[i][1]]++; rem0[R[i][2]]++; }
}

/* --------------------------------------------------------- main */
static void print_info(void)
{
    long i, tw1 = 0, tw2 = 0, full = 0;
    for (i = 0; i < ST1->nreps; i++) { tw1 += ST1->reps[i].weight; if (ST1->reps[i].weight == 720) full++; }
    for (i = 0; i < ST2->nreps; i++) tw2 += ST2->reps[i].weight;
    printf("Etage 1 : %ld affectations faisables de W1 (20 triplets), %ld orbites sous S6 "
           "(somme poids = %ld, orbites pleines de 720 : %ld)\n", ST1->nassign, ST1->nreps, tw1, full);
    printf("Etage 2 : %ld affectations faisables de W2 (6 triplets), %ld orbites sous S2 "
           "(somme poids = %ld)\n", ST2->nassign, ST2->nreps, tw2);
    printf("Nombre de jobs (representants de l'etage 1) = %ld ; sous-comptages par job = %ld\n",
           ST1->nreps, ST2->nreps);
    printf("Gain nominal sur les configurations de depart : %.1f\n",
           (double)ST1->nassign * ST2->nassign / ((double)ST1->nreps * ST2->nreps));
    printf("Triplets restants a parcourir : %d\n", NR);
}

int main(int argc, char **argv)
{
    long J, id;
    clock_t c0 = clock();
    init_basic();
    build_stages();
    build_order();
    init_rem0();
    J = ST1->nreps;
    fprintf(stderr, "[orbites construites en %.1f s]\n", (double)(clock() - c0) / CLOCKS_PER_SEC);

    if (argc >= 2 && strcmp(argv[1], "--info") == 0) { print_info(); return 0; }

    if (argc >= 3 && strcmp(argv[1], "--sample") == 0) {
        int K = atoi(argv[2]), i;
        long step = J / K;
        double t0 = (double)clock() / CLOCKS_PER_SEC, dt;
        double sumT = 0.0;
        for (i = 0; i < K; i++) {
            id = (long)i * step + step / 2;
            sumT += (double)ST1->reps[id].weight * (double)run_job(id);
        }
        dt = (double)clock() / CLOCKS_PER_SEC - t0;
        printf("%d jobs echantillonnes sur %ld : %.2f s  (%.4f s/job)\n", K, J, dt, dt / K);
        printf("Temps total estime (1 coeur) : %.0f s = %.1f h\n", dt / K * J, dt / K * J / 3600.0);
        printf("Estimation de N (tres approximative) : %.3e\n", sumT / K * J);
        return 0;
    }

    if (argc >= 3 && strcmp(argv[1], "--sum") == 0) {
        char *done = (char *)calloc((size_t)J, 1);
        unsigned long long total = 0;
        long ndone = 0, missing = 0;
        int f;
        for (f = 2; f < argc; f++) {
            FILE *fp = fopen(argv[f], "r");
            long jid; int w; u64 T;
            if (!fp) { fprintf(stderr, "ouverture impossible : %s\n", argv[f]); continue; }
            while (fscanf(fp, "%ld %d %llu", &jid, &w, &T) == 3) {
                if (jid >= 0 && jid < J && !done[jid]) {
                    done[jid] = 1; ndone++;
                    total += (u64)w * T;
                }
            }
            fclose(fp);
        }
        for (id = 0; id < J; id++) if (!done[id]) missing++;
        printf("jobs termines : %ld / %ld   manquants : %ld\n", ndone, J, missing);
        printf("TOTAL (partiel si manquants) = %llu\n", total);
        if (missing == 0)
            printf("Comparaison avec N attendu 125047355835960 : %s\n",
                   total == 125047355835960ULL ? "IDENTIQUE" : "DIFFERENT");
        return 0;
    }

    {
        int part = 0, nparts = 1;
        char fname[64];
        char *done;
        FILE *fp;
        time_t t0 = time(NULL);
        long ndo = 0, nskip = 0;
        if (argc >= 3) { part = atoi(argv[1]); nparts = atoi(argv[2]); }
        sprintf(fname, "res6_%d.txt", part);
        done = (char *)calloc((size_t)J, 1);
        fp = fopen(fname, "r");
        if (fp) {
            long jid; int w; u64 T;
            while (fscanf(fp, "%ld %d %llu", &jid, &w, &T) == 3)
                if (jid >= 0 && jid < J) { done[jid] = 1; nskip++; }
            fclose(fp);
        }
        print_info();
        printf("partie %d/%d, fichier %s, deja faits : %ld\n", part, nparts, fname, nskip);
        fp = fopen(fname, "a");
        if (!fp) { perror("res"); return 1; }
        for (id = part; id < J; id += nparts) {
            u64 T;
            if (done[id]) continue;
            T = run_job(id);
            fprintf(fp, "%ld %d %llu\n", id, (int)ST1->reps[id].weight, T);
            ndo++;
            if (ndo % 2000 == 0) {
                fflush(fp);
                printf("[%lds] %ld jobs faits dans cette session (dernier id %ld)\n",
                       (long)(time(NULL) - t0), ndo, id);
            }
        }
        fclose(fp);
        printf("Termine. Utiliser --sum pour le total.\n");
    }
    return 0;
}