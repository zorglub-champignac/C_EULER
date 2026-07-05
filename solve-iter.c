#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <time.h>
/*
 * solve_final.c - Backtracking itératif EXACT de solve2.c
 * avec sauvegarde/reprise sur disque.
 *
 * solve2.c de référence :
 *   bt(ti):
 *     if ti==56: count++; return
 *     rem[p0,p1,p2]--                     (pc=0→1, SANS feasible ici)
 *     for k=0,1,2:                         (pc=1,2,3)
 *       pi = tp[ti][k]
 *       if cov[pi]<2:
 *         cov[pi]++
 *         if feasible(): bt(ti+1)
 *         cov[pi]--
 *     rem[p0,p1,p2]++                     (pc=4→leave)
 */

int tp[56][3];

void init_tp() {
    int pi=0, idx[8][8];
    memset(idx,-1,sizeof(idx));
    for(int i=0;i<8;i++) for(int j=i+1;j<8;j++) idx[i][j]=pi++;
    int ti=0;
    for(int a=0;a<8;a++) for(int b=a+1;b<8;b++) for(int c=b+1;c<8;c++){
        tp[ti][0]=idx[a][b]; tp[ti][1]=idx[a][c]; tp[ti][2]=idx[b][c]; ti++;
    }
}

static inline int feasible(int cov[28], int rem[28]) {
    for(int p=0;p<28;p++) if(2-cov[p]>rem[p]) return 0;
    return 1;
}

/* --- Frame : état d'un appel bt(ti) suspendu ---
 *   pc=0 : rem-- pas encore fait
 *   pc=1 : rem-- fait, tenter k=0
 *   pc=2 : tenter k=1
 *   pc=3 : tenter k=2
 *   pc=4 : rem++ à faire puis remonter
 *   chosen : paire pushée par ce niveau (-1=aucune)
 */
typedef struct { int ti, pc, chosen; } Frame;

#define MAXD 58
Frame     S[MAXD];
int       top;
int       cov[28], rem[28];
long long cnt;
long long nodes   = 0;
int bestUp = 0 ;

/* --- Sauvegarde / chargement --- */
#define SAVEFILE "solve_state.bin"
#define MAGIC    0xF1A17E54UL
#define MAGIC2    0xF3A37F53UL
char *  load_name = SAVEFILE ;
void save_state() {
    FILE *f = fopen(load_name,"wb");
    if(!f) { perror("save"); return; }
    unsigned long m = MAGIC2;
    fwrite(&m,   sizeof m,          1, f);
    fwrite(&cnt, sizeof cnt,        1, f);
    fwrite(&top, sizeof top,        1, f);
    fwrite(cov,  sizeof cov,        1, f);
    fwrite(rem,  sizeof rem,        1, f);
    fwrite(S,    sizeof(Frame), top+1, f);
	fwrite(&bestUp,sizeof(bestUp),   1, f);
	fwrite(&nodes,sizeof(nodes),	1, f);
    fclose(f);
    //fprintf(stderr,"[saved]  top=%d  count=%lld\n", top, cnt);
}

int load_state() {
    FILE *f = fopen(load_name,"rb");
    if(!f) return 0;
    unsigned long m;
    if(fread(&m,sizeof m,1,f)<1 || m!=MAGIC2) {
        fprintf(stderr,"Sauvegarde invalide.\n"); fclose(f); return 0;
    }
    fread(&cnt, sizeof cnt,        1, f);
    fread(&top, sizeof top,        1, f);
    fread(cov,  sizeof cov,        1, f);
    fread(rem,  sizeof rem,        1, f);
    fread(S,    sizeof(Frame), top+1, f);
	fread(&bestUp,sizeof(bestUp), 1,f);
	fread(&nodes,sizeof(nodes),	1, f);
    fclose(f);
    fprintf(stderr,"[loaded] top=%d  count=%lld\n", top, cnt);
	cnt--;
	// nodes ;
    return 1;
}

volatile int interrupted = 0;
void on_sig(int s) { (void)s; interrupted = 1; }

void run(long long next_sv,long long save_every) {

    time_t    t0      = time(NULL);

    for(;;) {
        Frame *f = &S[top];
		if(f->ti < bestUp) {
			clock_t end = clock();
			printf("BestUp=%d count=%lld in %llds\n",bestUp,cnt,(long long)(time(NULL)-t0));
			fflush(stdout);
			bestUp = f->ti ;
		}

        /* Feuille */
        if(f->ti == 56) {
			if (cnt == 0) bestUp = 56 ;
            cnt++;
			if(cnt >= next_sv) {
                        fprintf(stderr,"[%llds] nodes=%lld  count=%lld  depth=%d\n",
                                (long long)(time(NULL)-t0), nodes, cnt, top);
                        save_state();
                        next_sv += save_every;
            }

            if(top==0) break;
            top--;
            if(S[top].chosen>=0) { cov[S[top].chosen]--; S[top].chosen=-1; }
            S[top].pc++;
            continue;
        }

        int ti=f->ti, p0=tp[ti][0], p1=tp[ti][1], p2=tp[ti][2];

        /* PC=0 : ENTER - rem-- sans feasible (exactement solve2.c) */
        if(f->pc == 0) {
            rem[p0]--; rem[p1]--; rem[p2]--;
            f->pc = 1;
            continue;
        }

        /* PC=1,2,3 : tenter k=pc-1 */
        if(f->pc <= 3) {
            //int k=f->pc-1 ; // , pi=tp[ti][k];
 			switch(f->pc) {
				case 1:
    				if ((cov[p0] < 2) &&(rem[p1]+cov[p1] >= 2) && (rem[p2]+cov[p2] >= 2)) {
						f->chosen = p0;
						cov[p0]++ ;
						top++;
						S[top] = (Frame){ ti+1, 0, -1 };
						nodes++;
						if(interrupted) { save_state(); return; }
						continue;
					}
					break ;
				case 2:
					if ((cov[p1] < 2) && (rem[p0]+cov[p0] >= 2) && (rem[p2]+cov[p2] >= 2)) {
						f->chosen = p1;
						cov[p1]++ ;
						top++;
						S[top] = (Frame){ ti+1, 0, -1 };
						nodes++;
						if(interrupted) { save_state(); return; }
						continue;
					}
					break ;
				case 3:
					if ((cov[p2] < 2) && (rem[p0]+cov[p0] >= 2) && (rem[p1]+cov[p1] >= 2)) {
						f->chosen = p2;
						cov[p2]++ ;
						top++;
						S[top] = (Frame){ ti+1, 0, -1 };
						nodes++;
						if(interrupted) { save_state(); return; }
						continue;
					}
					break ;
            } // end switch
            f->pc++;
            continue;
        }

        /* PC=4 : LEAVE - rem++ puis remonter */
        rem[p0]++; rem[p1]++; rem[p2]++;
        if(top==0) break;
        top--;
        if(S[top].chosen>=0) { cov[S[top].chosen]--; S[top].chosen=-1; }
        S[top].pc++;
    }

    printf("Solutions: %lld\n", cnt);
    remove(SAVEFILE);
}

int main(int argc, char **argv) {
    init_tp();
    signal(SIGINT, on_sig);
    long long save_every = 100000000LL;
    if(argc >= 2) save_every = atoll(argv[1]);
	if(argc >= 3) load_name = argv[2] ;
	long long next_sv = save_every ;

    if(!load_state()) {
        memset(cov, 0, sizeof cov);
        for(int p=0;p<28;p++) rem[p]=6;
        top    = 0;
        S[0]   = (Frame){ 0, 0, -1 };
        cnt    = 0;
        fprintf(stderr,"Démarrage from scratch.\n");
	}
	if (cnt != 0) {
		fprintf(stderr,"Démarrage from count=%lld.\n",cnt);
		next_sv = cnt + save_every - (cnt % save_every) ;
	}
    run(next_sv,save_every);
    return 0;
}