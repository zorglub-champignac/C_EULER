/*
 *
 *
 * orbit_count.c
 *
 * Comptage exact du nombre N de fonctions f : triplets(8) -> paires(8)
 * (f(t) inclus dans t, chaque paire image de exactement 2 triplets)
 * par la methode des orbites de configurations partielles.
 *
 * PRINCIPE (exact, sans hypothese de divisibilite)
 * ------------------------------------------------
 * Etage 1 : W1 = les 10 triplets de {0,1,2,3,4}, groupe G1 = S5 (120 perms
 *           de {0..4}, 5,6,7 fixes).
 * Etage 2 : W2 = {5,6,7}, {0,5,6}, {0,5,7}, {0,6,7}, groupe G2 = S3 sur
 *           {5,6,7} (6 perms, 0..4 fixes). G2 fixe point par point toute
 *           affectation de W1, donc G1 et G2 se composent.
 *
 *   N = somme sur les representants r1 des orbites de G1 sur les
 *       affectations faisables de W1, de  |O1| *
 *       somme sur les representants r2 des orbites de G2 sur les
 *       affectations faisables de W2, de  |O2| * count(r1, r2)
 *
 * ou count(r1,r2) = nombre de completions sur les 42 triplets restants
 * (backtracking avec elagage cov/remaining, test incremental).
 * Facteur nominal : 120 * 6 = 720 (les orbites courtes sont ponderees par
 * leur vraie taille : le resultat reste exact).
 *
 * UTILISATION
 *   gcc -O3 -march=native -o orbit_count orbit_count.c
 *   ./orbit_count --info                  structure des jobs
 *   ./orbit_count --bench ID SECONDES     essai de vitesse sur un job
 *   ./orbit_count                         tous les jobs (fichier res_0.txt)
 *   ./orbit_count PART NPARTS             jobs d'id % NPARTS == PART
 *                                         (fichier res_PART.txt), reprise
 *                                         automatique si le fichier existe
 *   ./orbit_count --sum res_0.txt res_1.txt ...   total et jobs manquants
 *
 * Le resultat final doit valoir N = 125 047 355 835 960 (valeur attendue).
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

// table des paires et des triplets
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
// cherche l'index en remettant dans l'ordre a<b<c
static int tindex3(int a, int b, int c)
{
    int t;
    if (a > b) { t = a; a = b; b = t; }
    if (b > c) { t = b; b = c; c = t; }
    if (a > b) { t = a; a = b; b = t; }
    return tidx[a][b][c];
}

/* ------------------------------------------------ permutations */
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

/* genere les permutations de {lo..lo+m-1}, autres elements fixes */
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
typedef struct {
    unsigned char ch[10];   /* choix (0,1,2) pour chaque triplet de W */
    int weight;             /* taille de l'orbite */
} Rep;

static int  S_n, S_ng;
static const int *S_W;
static int  S_pos[720][10];
static unsigned char S_cmap[720][10][3];
static int  pow3[11];
static int  S_cov[28];
static unsigned char S_ch[10];
static Rep *S_reps;
static int  S_nreps;
static long S_nassign;

static void stage_leaf(void)
{
    int g, i, j, code = 0, mincode, nmin = 0;
    for (i = 0; i < S_n; i++) code += S_ch[i] * pow3[i];
    mincode = code;
    for (g = 0; g < S_ng; g++) {
        int c2 = 0;
        unsigned char img[10];
        for (i = 0; i < S_n; i++) {
            j = S_pos[g][i];
            img[j] = S_cmap[g][i][S_ch[i]];
        }
        for (i = 0; i < S_n; i++) c2 += img[i] * pow3[i];
        if (c2 < mincode) { mincode = c2; }
        if (c2 == code) nmin++;
    }
    S_nassign++;
    if (code == mincode) {
        Rep *r = &S_reps[S_nreps++];
        memcpy(r->ch, S_ch, S_n);
        r->weight = S_ng / nmin;      /* orbite = |G| / |stabilisateur| */
    }
}

static void stage_dfs(int i)
{
    int c, p;
    if (i == S_n) { stage_leaf(); return; }
    for (c = 0; c < 3; c++) {
        p = tp[S_W[i]][c];
        if (S_cov[p] < 2) {
            S_cov[p]++; S_ch[i] = (unsigned char)c;
            stage_dfs(i + 1);
            S_cov[p]--;
        }
    }
}

static int stage_build(const int *W, int n, int perms[][8], int ng,
                       Rep *reps, long *nassign)
{
    int g, i, j, c, k;
    S_n = n; S_ng = ng; S_W = W; S_reps = reps; S_nreps = 0; S_nassign = 0;
    pow3[0] = 1;
    for (i = 1; i <= 10; i++) pow3[i] = pow3[i - 1] * 3;
    for (g = 0; g < ng; g++)
        for (i = 0; i < n; i++) {
            int *t = tri[W[i]];
            int gi = tindex3(perms[g][t[0]], perms[g][t[1]], perms[g][t[2]]);
            for (j = 0; j < n; j++) if (W[j] == gi) break;
            if (j == n) { fprintf(stderr, "W non invariant sous le groupe\n"); exit(1); }
            S_pos[g][i] = j;
            for (c = 0; c < 3; c++) {
                int p = tp[W[i]][c];
                int q = pid_[perms[g][P2[p][0]]][perms[g][P2[p][1]]];
                for (k = 0; k < 3; k++) if (tp[W[j]][k] == q) break;
                if (k == 3) { fprintf(stderr, "paire image absente\n"); exit(1); }
                S_cmap[g][i][c] = (unsigned char)k;
            }
        }
    memset(S_cov, 0, sizeof(S_cov));
    stage_dfs(0);
    {
        long sw = 0;
        for (i = 0; i < S_nreps; i++) sw += reps[i].weight;
        if (sw != S_nassign) {
            fprintf(stderr, "ERREUR: somme des poids %ld != %ld affectations\n", sw, S_nassign);
            exit(1);
        }
    }
    *nassign = S_nassign;
    return S_nreps;
}

/* ------------------------------------------------- donnees globales */
static int  W1[10], W2[4];
static Rep  reps1[40000], reps2[200];
static int  n1, n2;
static long nass1, nass2;

static int  R[48][3];     /* paires des triplets restants (ordre de parcours) */
static int  NR;
static int  R_tr[48];

static void build_stages(void)
{
    static int perms1[720][8], perms2[720][8];
    int ng1, ng2, a, b, c, k = 0;
    for (a = 0; a < 5; a++)
        for (b = a + 1; b < 5; b++)
            for (c = b + 1; c < 5; c++) W1[k++] = tidx[a][b][c];
    W2[0] = tidx[5][6][7];
    W2[1] = tidx[0][5][6];
    W2[2] = tidx[0][5][7];
    W2[3] = tidx[0][6][7];
    ng1 = gen_group(0, 5, perms1);
    ng2 = gen_group(5, 3, perms2);
    n1 = stage_build(W1, 10, perms1, ng1, reps1, &nass1);
    n2 = stage_build(W2, 4, perms2, ng2, reps2, &nass2);
}

/* ordre des triplets restants : glouton (ferme les paires au plus tot) */
static void build_order(void)
{
    int used[56], i, j, k, p, best, bestscore;
    int remtot[28], remleft[28];
    memset(used, 0, sizeof(used));
    for (i = 0; i < 10; i++) used[W1[i]] = 1;
    for (i = 0; i < 4; i++) used[W2[i]] = 1;
    memset(remtot, 0, sizeof(remtot));
    for (i = 0; i < 56; i++) if (!used[i])
        for (k = 0; k < 3; k++) remtot[tp[i][k]]++;
    memcpy(remleft, remtot, sizeof(remtot));
    NR = 56 - 14;
    for (j = 0; j < NR; j++) {
        best = -1; bestscore = -1000000;
        for (i = 0; i < 56; i++) if (!used[i]) {
            int sc = 0;
            for (k = 0; k < 3; k++) {
                p = tp[i][k];
                if (remleft[p] == 1) sc += 10;          /* ferme la paire */
                if (remleft[p] == remtot[p]) sc -= 1;   /* ouvre la paire */
            }
            if (sc > bestscore) { bestscore = sc; best = i; }
        }
        used[best] = 1;
        R_tr[j] = best;
        for (k = 0; k < 3; k++) { R[j][k] = tp[best][k]; remleft[tp[best][k]]--; }
    }
}

/* ------------------------------------------------ backtracking */
static int  cov[28], rem[28];
static u64  cnt, nodes;
static int  abort_flag;
static double t_limit, t_start;

static void bt(int k)
{
    int a, b, c, j, p;
    if (k == NR) { cnt++; return; }
    if (abort_flag) return;
    a = R[k][0]; b = R[k][1]; c = R[k][2];
    rem[a]--; rem[b]--; rem[c]--;
    for (j = 0; j < 3; j++) {
        p = R[k][j];
        if (cov[p] < 2) {
            cov[p]++;
            /* test incremental : seules les 3 paires du triplet ont change */
            if (cov[a] + rem[a] >= 2 && cov[b] + rem[b] >= 2 && cov[c] + rem[c] >= 2)
                bt(k + 1);
            cov[p]--;
        }
    }
    rem[a]++; rem[b]++; rem[c]++;
    if ((++nodes & 0xFFFFFFULL) == 0 && t_limit > 0) {
        if ((double)(clock() - 0) / CLOCKS_PER_SEC - t_start > t_limit) abort_flag = 1;
    }
}

/* prepare cov/rem pour le job id ; retourne 0 si infaisable d'emblee */
static int setup_job(long id, int *w1, int *w2)
{
    int r1 = (int)(id / n2), r2 = (int)(id % n2), i, p;
    memset(cov, 0, sizeof(cov));
    for (i = 0; i < 10; i++) cov[tp[W1[i]][reps1[r1].ch[i]]]++;
    for (i = 0; i < 4; i++) cov[tp[W2[i]][reps2[r2].ch[i]]]++;
    memset(rem, 0, sizeof(rem));
    for (i = 0; i < NR; i++) { rem[R[i][0]]++; rem[R[i][1]]++; rem[R[i][2]]++; }
    for (p = 0; p < 28; p++)
        if (cov[p] > 2 || cov[p] + rem[p] < 2) return 0;
    *w1 = reps1[r1].weight; *w2 = reps2[r2].weight;
    return 1;
}

static u64 run_job(long id, int *weight, double maxsec)
{
    int w1 = 0, w2 = 0;
    cnt = 0; nodes = 0; abort_flag = 0;
    t_limit = maxsec; t_start = (double)clock() / CLOCKS_PER_SEC;
    if (!setup_job(id, &w1, &w2)) { *weight = 0; return 0; }
    *weight = w1 * w2;
    bt(0);
    return cnt;
}

/* --------------------------------------------------------- main */
static void print_info(void)
{
    long totw1 = 0, totw2 = 0;
    int i, c1[130] = {0};
    for (i = 0; i < n1; i++) { totw1 += reps1[i].weight; c1[reps1[i].weight]++; }
    for (i = 0; i < n2; i++) totw2 += reps2[i].weight;
    printf("Etage 1 : %ld affectations faisables de W1, %d orbites sous S5 (somme poids = %ld)\n",
           nass1, n1, totw1);
    printf("Etage 2 : %ld affectations faisables de W2, %d orbites sous S3 (somme poids = %ld)\n",
           nass2, n2, totw2);
    printf("Nombre de jobs = %d x %d = %ld\n", n1, n2, (long)n1 * n2);
    printf("Gain nominal sur le nombre de configurations de depart : %.1f\n",
           (double)nass1 * nass2 / ((double)n1 * n2));
    printf("Triplets restants a parcourir : %d\n", NR);
}

int main(int argc, char **argv)
{
    long J, id;
    init_basic();
    build_stages();
    build_order();
    J = (long)n1 * n2;

    if (argc >= 2 && strcmp(argv[1], "--info") == 0) { print_info(); return 0; }

    if (argc >= 4 && strcmp(argv[1], "--bench") == 0) {
        int w; u64 c; double t0, dt;
        id = atol(argv[2]);
        t0 = (double)clock() / CLOCKS_PER_SEC;
        c = run_job(id, &w, atof(argv[3]));
        dt = (double)clock() / CLOCKS_PER_SEC - t0;
        printf("job %ld : poids=%d  count=%llu  noeuds=%llu  temps=%.1fs  %s  (%.2e sol/s)\n",
               id, w, c, nodes, dt, abort_flag ? "INTERROMPU" : "termine", dt > 0 ? c / dt : 0.0);
        return 0;
    }

    if (argc >= 3 && strcmp(argv[1], "--sum") == 0) {
        char *done = (char *)calloc((size_t)J, 1);
        unsigned long long total = 0;
        long ndone = 0, missing = 0;
        int f;
        for (f = 2; f < argc; f++) {
            FILE *fp = fopen(argv[f], "r");
            long jid; int w; u64 c;
            if (!fp) { fprintf(stderr, "ouverture impossible : %s\n", argv[f]); continue; }
            while (fscanf(fp, "%ld %d %llu", &jid, &w, &c) == 3) {
                if (jid >= 0 && jid < J && !done[jid]) {
                    done[jid] = 1; ndone++;
                    total += (u64)w * c;
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
        sprintf(fname, "res_%d.txt", part);
        done = (char *)calloc((size_t)J, 1);
        fp = fopen(fname, "r");
        if (fp) {
            long jid; int w; u64 c;
            while (fscanf(fp, "%ld %d %llu", &jid, &w, &c) == 3)
                if (jid >= 0 && jid < J) { done[jid] = 1; nskip++; }
            fclose(fp);
        }
        print_info();
        printf("partie %d/%d, fichier %s, deja faits : %ld\n", part, nparts, fname, nskip);
        fp = fopen(fname, "a");
        if (!fp) { perror("res"); return 1; }
        for (id = part; id < J; id += nparts) {
            int w; u64 c;
            if (done[id]) continue;
            c = run_job(id, &w, 0.0);
            fprintf(fp, "%ld %d %llu\n", id, w, c);
            fflush(fp);
            ndo++;
            if (ndo % 20 == 0)
                printf("[%lds] %ld jobs faits dans cette session (dernier id %ld)\n",
                       (long)(time(NULL) - t0), ndo, id);
        }
        fclose(fp);
        printf("Termine. Utiliser --sum pour le total.\n");
    }
    return 0;
}
