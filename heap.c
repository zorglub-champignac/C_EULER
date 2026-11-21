/* -----------------------------------------------------------------------------
 * Fonction de tri suivant la methode DHEAPSORT.
 *
 * La structure HEAP_DATA permet de gerer une pile d'elements ou le meilleur
 * element est toujours sur le sommet de la pile. On peut a tout moment retirer
 * le meilleur element de la pile (HEAP_GetFirst) ou inserer un nouvel element
 * (HEAP_Insert)
 * Attention : si la pile est pleine l'insertion du nouvel element n'est pas
 * faite.
 * La structure DHEAP_DATA offre les memes fonctionalites mais en plus si la
 * pile est pleine DHEAP_Insert insere le nouvel element et retourne le plus
 * mauvais contenu dans la pile, ce qui garantit que la pile contient toujours
 * les meilleurs elements.
 * Remarque : les elements ne sont pas physiquement inseres dans la pile,
 * on garde un pointeur sur eux. Le programme appelant proprietaire des
 * elements ne doit donc pas les modifier a l'insu  des routines de gestion
 * des piles.
 * Pour DHEAP_DATA en plus du classement up, on gere un classement down.
 * Ce classement est gere grace a la correspondance upToDown et sa reciproque
 * downToUp.
 * On a donc toujours :
 *		UpToDown(DownToUp()) = Identite()
 *		DownToUp(UpToDown()) = Identite()
 *		comp( hp->data[i] , hp->data[j] ) meme signe que (j - i) (class up)
 *		comp( hp->data[DownToUp(i)] , hp->data[DownToUp(j)] ) de signe
 *		contraire de (j - i) (classement down ).
 *
 * HEAP_INDEX a un fonctionnement identiqu mais les donnees sont gerees par index (> 0)  au lieu
 * de pointeur 
 * -----------------------------------------------------------------------------
 */

/* ----------------------
 * Directives d'inclusion
 * ----------------------
 */
#include <stdio.h>		/* Fcts systeme */
#include <stdlib.h>

#include "heap.h"



struct HEAP_DATA {
	void        		**data;			/* pointeur sur les elements 	*/
	uint32_t            sizeMax;	    /* taille maximum				*/
	uint32_t	        size;			/* taille effective				*/
	int32_t        		isTopFirst;		/* Premier element present	?	*/
	int				(*comp) (const void *el1,const void *el2, void * ptUser) ;
	int				(*compSimple) (const void *el1,const void *el2) ;
	void * ptUser;
										/* Routine de comparaison		*/
} ;

struct DHEAP_DATA {
	void        		**data ;		/* pointeur sur les elements 	*/
	uint32_t	        sizeMax;        /* taille maximum				*/
	uint32_t	        size;			/* taille effective				*/
	int32_t        		isTopFirst;		/* Premier element present	?	*/
	int         	(*comp) (const void *el1,const void *el2, void * ptUser) ;
										/* Routine de comparaison		*/
	int				(*compSimple) (const void *el1,const void *el2);
	uint32_t	        freeDown;		/* indice de l'el. libre decr.	*/
	uint32_t	        *upToDown;		/* Corresp. croiss. decroiss.	*/
	uint32_t	        *downToUp;		/* Corresp. reciproque			*/
	void * ptUser;
} ;



static void HEAP_InsertFromTop(HEAP_DATA *hp,const void *element) ;
static void DHEAP_InsertFromTop(DHEAP_DATA *hp,const void *element) ;
static void DHEAP_InsertFromDown(DHEAP_DATA *hp) ;
static int HEAP_LocalComp(const void * el1, const void * el2, void * ptUser);
static int DHEAP_LocalComp(const void * el1, const void * el2, void * ptUser);


/* ---------------------------
 * Code des Fonctions internes
 * ---------------------------
 */

/* -----------------------------------------------------------------------------
 * Fonction	: HEAP_InsertFromTop
 *
 * Role		:
 *  Insere l'element a partir du sommet (rg 0 ) en le faisant
 *  descendre jusqu'a sa place
 *
 * Retour	:
 *	- neant
 * -----------------------------------------------------------------------------
 */
static void HEAP_InsertFromTop(
	HEAP_DATA	*hp,
	const void	*element)
{
	uint32_t i=0 , j = 1 ;
	void * ptUser;
	/* i est la position courante, j est le premier fils de i */
	if (hp->compSimple == NULL) {
		ptUser = hp->ptUser;
	} else {
		ptUser = (void *) hp;
	}
	while(j < hp->size) {
	
		/* on regarde le plus grand des 2 fils */
		if( ( j+1 <  hp->size ) &&
			(hp->comp (hp->data[j],hp->data[j+1], ptUser) < 0) ) {
			j++ ;
		}
		
		/* on regarde si l'element est plus grand que les fils */
		if(hp->comp(element,hp->data[j] ,ptUser) < 0 ) {
			/* non , on continue la descente en faisant remonter le plus
			grand des fils  */
			hp->data[i] = hp->data[j] ;
			i = j ;
			j += j + 1 ;
		} else {
			/* oui on arrete la */
			break ;
		}
	}
	
	/* on insere reellement l'element a sa place */
	hp->data[i] = (void *) element ;
	return ;
}

/* -----------------------------------------------------------------------------
 * Fonction	: DHEAP_InsertFromTop
 *
 * Role		:
 *  Insere l'element a partir du sommet (rg 0 ) en le faisant
 *  descendre jusqu'a sa place dans le classement up.
 *
 * Retour	:
 *	- neant
 * -----------------------------------------------------------------------------
 */
static void DHEAP_InsertFromTop(
	DHEAP_DATA *dhp,
	const void *element)
{
	uint32_t iu = 0 , ju = 1 ;
	void * ptUser;
	/* iu est la position courante, ju est le premier fils de iu dans le
	classement up. iud est la position initiale dans le classement down */
	uint32_t iud = dhp->upToDown[iu] ;
	
	if (dhp->compSimple == NULL){
		ptUser = dhp->ptUser;
	} else {
		ptUser = (void *) dhp;
	}
	while(ju < dhp->size) {
	
		/* on regarde le plus grand des 2 fils */
		if( ( ju + 1 <  dhp->size ) &&
			(dhp->comp (dhp->data[ju],dhp->data[ju + 1] , ptUser) < 0) ) {
			ju++ ;
		}
		
		/* on regarde si l'element est plus grand que les fils */
		if(dhp->comp(element,dhp->data[ju],ptUser) < 0 ) {
			/* non , on continue la descente en faisant remonter le plus
			grand des de fils  */
			dhp->data[iu] = dhp->data[ju] ;
			/* mise a jour de la correspondance up <-> down */
			dhp->upToDown[iu] = dhp->upToDown[ju] ;
			dhp->downToUp[dhp->upToDown[iu]] = iu ;
			iu = ju ;
			ju += ju + 1 ;
		} else {
			/* oui on arrete la */
			break ;
		}
	}
	
	/* on insere reellement l'element */
	dhp->data[iu ] = (void *) element ;
	
	/* mise a jour de la correspondance up <-> down */
	dhp->upToDown[iu] = iud ;
	dhp->downToUp[iud] = iu ;
	return ;
}


/* -----------------------------------------------------------------------------
 * Fonction	: DHEAP_InsertFromDown
 *
 * Role		:
 *  Met a sa place dans le classement down, l'element de rang dhp->freeDown
 *  On ne sait pas a priori si l'element doit remonter ou redescendre
 *
 * Retour	:
 *	- neant
 * -----------------------------------------------------------------------------
 */
static void DHEAP_InsertFromDown(DHEAP_DATA *dhp)
{
	uint32_t id = dhp->freeDown ;
	uint32_t idu = dhp->downToUp[id] ;
	uint32_t jd = ( id -1 ) >> 1 ;
	void * ptUser;
	if (dhp->compSimple == NULL){
		ptUser = dhp->ptUser;
	} else {
		ptUser = dhp;
	}
	/* id position courante , jd position du pere , idu position
	initiale dans le classement up */
	if(  (id > 0) &&
		 (dhp->comp(dhp->data[idu],dhp->data[dhp->downToUp[jd]], ptUser) < 0) ) {
		
		/* il faut faire remonter l'lement */
		do {
			/* on remonte l'element tant qu'il est inferieur au pere */
			dhp->downToUp[id] = dhp->downToUp[jd] ;
			dhp->upToDown[dhp->downToUp[id]] = id ;
			id = jd ;
			jd = (id -1) >> 1 ;
		} while ((id > 0) &&
				(dhp->comp(dhp->data[idu],dhp->data[dhp->downToUp[jd]], ptUser) < 0) );
		
		/* on insere l'element a sa place finale */
		dhp->downToUp[id] = idu ;
		dhp->upToDown[idu] = id ;
	
	} else {
		/* il faut faire descendre l'element, jd est le premier fils */
		jd = ( id << 1 ) + 1 ;
	
		while(jd < dhp->size) {
		
			/* on regarde le plus petit des 2 fils */
			if( ( jd + 1 <  dhp->size ) &&
				(dhp->comp (dhp->data[dhp->downToUp[jd]],
						   dhp->data[dhp->downToUp[jd + 1]],ptUser) > 0) ) {
				jd++ ;
			}
			
			/* on regarde si l'element est plus petit que les fils */
			if(dhp->comp(dhp->data[idu],dhp->data[dhp->downToUp[jd]],ptUser) > 0 ) {
				/* non , on contine la descente */
				dhp->downToUp[id] = dhp->downToUp[jd] ;
				dhp->upToDown[dhp->downToUp[id]] = id ;
				id = jd ;
				jd += jd + 1 ;
			} else {
				/* oui on arrete la */
				break ;
			}
		}
		
		/* on insere l'element a sa place finale */
		dhp->downToUp[id] = idu ;
		dhp->upToDown[idu] = id ;
	}
	return ;
}



static int HEAP_LocalComp(const void * el1, const void * el2, void * ptUser)
{
	HEAP_DATA * hp;
	hp = (HEAP_DATA *) ptUser;
	return (hp->compSimple(el1,el2));
}

static int DHEAP_LocalComp(const void * el1, const void * el2, void * ptUser)
{
	DHEAP_DATA * dhp;
	dhp = (DHEAP_DATA *) ptUser;
	return (dhp->compSimple(el1,el2));
}


/* ----------------------------
 * Code des Fonctions exportees
 * ----------------------------
 */

/* -----------------------------------------------------------------------------
 * Fonction	: HEAP_Alloc
 *
 * Role		:
 *	Allocation d'une pile simple de taille sizeMax
 *  La pile est initialisee
 *
 * Retour	:
 *	- NULL en cas d'echec
 *  - ptr sur la pile sinon
 * -----------------------------------------------------------------------------
 */
HEAP_DATA * HEAP_Alloc(
	const uint32_t sizeMax,
	int (*comp) (const void *el1,const void *el2) )
{
	HEAP_DATA *newHeap = (HEAP_DATA *) calloc(1, sizeof(HEAP_DATA));

	if(newHeap == NULL) return newHeap ;
	
	newHeap->data = (void **) malloc(sizeMax*sizeof(void *));
	
	if(newHeap->data == NULL) {
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->sizeMax = sizeMax ;
	newHeap->comp = HEAP_LocalComp ;
	newHeap->compSimple = comp;
	newHeap->isTopFirst = 0 ;

	return newHeap ;
}


HEAP_DATA * HEAP_AllocU(
	const uint32_t sizeMax,
	int (*comp) (const void *el1,const void *el2, void * ptUser),
	void * ptUser)
{
	HEAP_DATA *newHeap = (HEAP_DATA *) calloc(1, sizeof(HEAP_DATA));

	if(newHeap == NULL) return newHeap ;
	
	newHeap->data = (void **) malloc(sizeMax*sizeof(void *));
	
	if(newHeap->data == NULL) {
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->sizeMax = sizeMax ;
	newHeap->comp = comp ;
	newHeap->compSimple = NULL;
	newHeap->isTopFirst = 0 ;
	newHeap->ptUser = ptUser;

	return newHeap ;
}


/* -----------------------------------------------------------------------------
 * Fonction	: HEAP_Realloc
 *
 * Role		:
 *	Reallocation d'une nouvelle taille
 *
 * Retour	:
 *	- NULL en cas d'echec
 * -----------------------------------------------------------------------------
 */
HEAP_DATA * HEAP_Realloc(
	HEAP_DATA *hp,			/* Structure a desallouer	*/
	const int32_t newSize)	/* Nouvelle taille			*/
{
	if(hp == NULL) return hp;
	{
		void * newData = realloc(hp->data,sizeof(void *) * newSize) ;

		if(newData == NULL) {
			free(hp->data);
			free(hp) ;
			return NULL ;
		}
		
		hp->data = (void **) newData ;
		hp->sizeMax = newSize ;
		return hp ;
	}
}

/* -----------------------------------------------------------------------------
 * Fonction	: DHEAP_Alloc
 *
 * Role		:
 *	Allocation d'une pile double de taille sizeMax
 *  La pile est initialisee
 *
 * Retour	:
 *	- NULL en cas d'echec
 *  - ptr sur la pile sinon
 * -----------------------------------------------------------------------------
 */
DHEAP_DATA * DHEAP_Alloc(
	const uint32_t sizeMax,
	int (*comp) (const void *el1,const void *el2) )
{
	DHEAP_DATA *newHeap = (DHEAP_DATA *) calloc(1, sizeof(DHEAP_DATA));

	if(newHeap == NULL) return newHeap ;
	
//	newHeap->data = (void **) calloc(sizeMax, sizeof(void *));
	newHeap->data = (void **) malloc(sizeMax*sizeof(void *));
	
	if(newHeap->data == NULL) {
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->upToDown = (uint32_t *) calloc(sizeMax, sizeof(newHeap->upToDown[0])) ;
	
	if(newHeap->upToDown == NULL) {
		free(newHeap->data);
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->downToUp = (uint32_t *) calloc(sizeMax, sizeof(newHeap->downToUp[0])) ;
	
	if(newHeap->downToUp == NULL) {
		free(newHeap->upToDown) ;
		free(newHeap->data);
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->sizeMax = sizeMax ;
	newHeap->comp = DHEAP_LocalComp ;
	newHeap->compSimple = comp;
	newHeap->isTopFirst = 0 ;
	newHeap->freeDown = 0 ;
	return newHeap ;
}


DHEAP_DATA * DHEAP_AllocU(
	const uint32_t sizeMax,
	int (*comp) (const void *el1,const void *el2, void * ptUser),
	void * ptUser)
{
	DHEAP_DATA *newHeap = (DHEAP_DATA *) calloc(1, sizeof(DHEAP_DATA));

	if(newHeap == NULL) return newHeap ;
	
	newHeap->data = (void **) malloc(sizeMax*sizeof(void *));
	
	if(newHeap->data == NULL) {
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->upToDown = (uint32_t *) calloc(sizeMax, sizeof(newHeap->upToDown[0])) ;
	
	if(newHeap->upToDown == NULL) {
		free(newHeap->data);
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->downToUp = (uint32_t *) calloc(sizeMax, sizeof(newHeap->downToUp[0])) ;
	
	if(newHeap->downToUp == NULL) {
		free(newHeap->upToDown) ;
		free(newHeap->data);
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->sizeMax = sizeMax ;
	newHeap->comp = comp ;
	newHeap->compSimple = NULL;
	newHeap->isTopFirst = 0 ;
	newHeap->freeDown = 0 ;
	newHeap->ptUser = ptUser;
	return newHeap ;
}

/* -----------------------------------------------------------------------------
 * Fonction	: DHEAP_Realloc
 *
 * Role		:
 *	Reallocation d'une nouvelle taille
 *
 * Retour	:
 *	- NULL en cas d'echec
 * -----------------------------------------------------------------------------
 */
DHEAP_DATA * DHEAP_Realloc(
	DHEAP_DATA *dhp,		/* Structure a desallouer	*/
	const int32_t newSize)	/* Nouvelle taille			*/
{
	void *newData ;

	if(dhp == NULL) return NULL ;
	
	newData = realloc(dhp->data,sizeof(void *) * newSize);
	
	if(newData == NULL) {
		DHEAP_Free(dhp) ;
		return NULL ;
	} else {
		dhp->data = (void **) newData ;
	}
	
	newData = realloc(dhp->upToDown,sizeof(dhp->upToDown[0])* newSize) ;
	
	if(newData == NULL) {
		DHEAP_Free(dhp) ;
		return NULL ;
	} else {
		dhp->upToDown = (uint32_t *) newData ;
	}

	newData = realloc(dhp->downToUp,sizeof(dhp->downToUp[0]) * newSize) ;
	
	if(newData == NULL) {
		DHEAP_Free(dhp) ;
		return NULL ;
	} else {
		dhp->downToUp = (uint32_t *) newData ;
	}
	
	dhp->sizeMax = newSize ;
	
	return dhp ;
}

/* -----------------------------------------------------------------------------
 * Fonction	: HEAP_GetFirst
 *
 * Role		:
 *	Retourne un pointeur sur le meilleur element dans la pile
 *
 * Retour	:
 *	- NULL si pile vide
 * -----------------------------------------------------------------------------
 */
void * HEAP_GetFirst(HEAP_DATA *hp)
{
	if(hp->size == 0) return NULL ;  /* pile vide */

	if(hp->isTopFirst) {
		hp->isTopFirst = 0 ;   /* le meilleure est en tete on le retourne */
	} else {
		/* on va faire remonter le dernier element pour recupere la tete */
		HEAP_InsertFromTop(hp,hp->data[hp->size]) ;
	}
	
	hp->size-- ;  /* mise a jour de la taille de la pile */
	
	return hp->data[0] ;
}



/* -----------------------------------------------------------------------------
 * Fonction	: DHEAP_GetFirst
 *
 * Role		:
 *	Retourne un pointeur sur le meilleur element dans la pile
 *
 * Retour	:
 *	- NULL si pile vide
 * -----------------------------------------------------------------------------
 */
void * DHEAP_GetFirst(DHEAP_DATA *dhp)
{
	
	if(dhp->size == 0) return NULL; /* pile vide */

	if(dhp->isTopFirst == 0) {
		/* on va faire remonter le dernier element pour recupere la tete */
		/*puisque l'element de rg dhp->size passe en zero dans le classement
		  up, on met a jour les correspondances entre le classement
		  up et down */
		dhp->upToDown[0] = dhp->upToDown[dhp->size] ;
		dhp->downToUp[dhp->upToDown[0]] = 0 ;
		
		/* on met a jour le classement up en remettant l'element en
		tete a sa place */
		DHEAP_InsertFromTop(dhp,dhp->data[dhp->size]) ;
		
		/* on reorganise le classement down en faisant passer l'element
		de rg dhp->size (ds le classement down) au rang dhp->freeDown
		(toujours ds le class. down), puisque il s'agit de l'element libere
		au coup precedent */
		if(dhp->freeDown < dhp->size) {
			/* on met a jour la correspondance entre up et down */
			dhp->downToUp[dhp->freeDown] = dhp->downToUp[dhp->size] ;
			dhp->upToDown[dhp->downToUp[dhp->freeDown]] = dhp->freeDown ;
			/* on reorganise le classement down */
			DHEAP_InsertFromDown(dhp);
		}
	} else {
		dhp->isTopFirst = 0 ;
	}
	
	dhp->size-- ;	/* mise a jour de la taille */
	
	/* on stocke l'emplacment libre */
	dhp->freeDown = dhp->upToDown[0] ;
	
	return dhp->data[0] ;
}

/* -----------------------------------------------------------------------------
 * Fonction	: HEAP_Insert
 *
 * Role		:
 *	Insere un element
 *  ATTENTION! si la pile est pleine, l'element n'est pas insere
 *
 * Retour	:
 *	- 0 si OK
 *  -1 si erreur
 * -----------------------------------------------------------------------------
 */
int  HEAP_Insert(
	HEAP_DATA *hp,
	const void *element)
{
	void * ptUser;
	if(hp->compSimple == NULL) {
		ptUser = hp->ptUser;
	} else {
		ptUser = hp;
	}
	if(hp->isTopFirst) {
		
		uint32_t j = hp->size++ ;

		if(j >= hp->sizeMax) {
			hp->size-- ; /* pile pleine on ne fait rien */
			return -1 ;
		}

		/* on a insere l'element en queue (hp->size), on le fait remonter
		a sa place en permuttant avec ses peres qui lui sont inferieur */
		while(j > 0) {
			/* rang du pere */
			uint32_t i = (j - 1) >> 1 ;
			if(hp->comp(element,hp->data[i], ptUser) > 0 ) {
				/* on permute */
				hp->data[j] = hp->data[i] ;
				j = i ;
			} else {
				break ;
			}
		}
		/* on insere reellement l'element a sa place */
		hp->data[j] = (void *) element ;
	} else { /* la place en tete est libre, on insere a partir de la tete */

		hp->isTopFirst = 1 ;
		hp->size++;
		HEAP_InsertFromTop(hp,element) ;
	}
    return 0 ;
}

/* -----------------------------------------------------------------------------
 * Fonction	: HEAP_Reorg
 *
 * Role		:
 *	Reordonnancement de la pile
 *
 * Retour	:
 *	- neant
 * -----------------------------------------------------------------------------
 */
void HEAP_Reorg(
	HEAP_DATA *hp,
	int (*newComp) (const void *el1,const void *el2))
{
	uint32_t nbOk ;

	hp->compSimple = newComp ;
	hp->comp = HEAP_LocalComp;
	if(hp->size <= 0) return;
	
	if(hp->isTopFirst == 0) {
		hp->data[0] = hp->data[hp->size] ;
		hp->isTopFirst = 1 ;
	}
	
	for(nbOk = 1 ; nbOk < hp->size; nbOk ++) {
		
		uint32_t j	= nbOk ;
		void *element		= hp->data[j] ;
		
		/* on insere l'element de rang nbOk a sa place en le faisant remonter
		a sa place en permuttant avec ses peres qui lui sont inferieur */
		while(j > 0) {
			/* rang du pere */
			uint32_t i = (j - 1) >> 1 ;

			if(newComp(element,hp->data[i]) > 0 ) {
				/* on permute */
				hp->data[j] = hp->data[i] ;
				j = i ;
			} else {
				break ;
			}
		}
		
		/* on insere reellement l'element a sa place */
		hp->data[j] = (void *) element ;
	}
}


void  HEAP_ReorgU(
	HEAP_DATA *hp,
	int (*newComp) (const void *el1,const void *el2, void * ptUser),
	void * ptUser)
{
	uint32_t nbOk ;

	hp->comp = newComp ;
	hp->ptUser = ptUser;
	hp->compSimple = NULL;
	if(hp->size <= 0) return;
	
	if(hp->isTopFirst == 0) {
		hp->data[0] = hp->data[hp->size] ;
		hp->isTopFirst = 1 ;
	}
	
	for(nbOk = 1 ; nbOk < hp->size; nbOk ++) {
		
		uint32_t j	= nbOk ;
		void *element		= hp->data[j] ;
		
		/* on insere l'element de rang nbOk a sa place en le faisant remonter
		a sa place en permuttant avec ses peres qui lui sont inferieur */
		while(j > 0) {
			/* rang du pere */
			uint32_t i = (j - 1) >> 1 ;

			if(newComp(element,hp->data[i],ptUser) > 0 ) {
				/* on permute */
				hp->data[j] = hp->data[i] ;
				j = i ;
			} else {
				break ;
			}
		}
		
		/* on insere reellement l'element a sa place */
		hp->data[j] = (void *) element ;
	}
}

/* -----------------------------------------------------------------------------
 * Fonction	: DHEAP_Insert
 *
 * Role		:
 *  Insere un element
 *  Si la pile est pleine, l'element est insere a la place du plus mauvais
 *  des elements qui est retourne
 *
 * Retour	:
 *	- NULL si la pile n'etait pas pleine
 * -----------------------------------------------------------------------------
 */
void * DHEAP_Insert(
	DHEAP_DATA *dhp,
	void *element)
{
	void * ptUser;
	if (dhp->compSimple == NULL) {
		ptUser = dhp->ptUser;
	} else {
		ptUser = dhp;
	}
	if(dhp->isTopFirst) {

		uint32_t ju = dhp->size++ ;
		uint32_t jd ;

		if(ju >= dhp->sizeMax) {
		
			/* la pile es pleine */
			dhp->size-- ;
			
			/* on compare l'element au pire dans la pile */
			if(dhp->comp(element,dhp->data[dhp->downToUp[0]],ptUser) < 0 ) {
				return element ; /* l'element est le pire */
			} else {
			
				/* on expulse l'element le pire */
				void * elementExpulse = dhp->data[dhp->downToUp[0]] ;
				ju = dhp->downToUp[0] ;
				
				/* on insere l'element a sa place, en le faisant remonter
				pusique le pire element etait forcement une feuille */
				while(ju > 0) {
					uint32_t iu = (ju - 1) >> 1 ;

					if(dhp->comp(element,dhp->data[iu],ptUser) > 0 ) {
						
						/* l'element est meilleurs que le pere */
						dhp->data[ju] = dhp->data[iu] ;
					
						/* mise a jour des corresp. up <-> down */
						dhp->upToDown[ju] = dhp->upToDown[iu] ;
						dhp->downToUp[dhp->upToDown[ju]] = ju ;
						ju = iu ;
					} else {
						break ;
					}
				}
				
				/* on insere reellement l'element a sa place */
				dhp->data[ju] = (void *) element ;
				dhp->downToUp[0] = ju ;
				dhp->upToDown[ju] = 0 ;
				
				/* il faut maintenant l'inserer dans la file down */
				dhp->freeDown = 0 ;
				DHEAP_InsertFromDown(dhp);
				return elementExpulse ;
			}
		}
		while(ju > 0) {
			
			/* la pile n'est pas pleine, on insere l'element a partir
			du rang dhp->size en le faisant remonter jusqu'a sa place */
			uint32_t iu = (ju - 1) >> 1 ;

			if(dhp->comp(element,dhp->data[iu],ptUser) > 0 ) {
			
				/* l'element est meilleur que le pere */
				dhp->data[ju] = dhp->data[iu] ;
				dhp->upToDown[ju] = dhp->upToDown[iu] ;
				dhp->downToUp[dhp->upToDown[ju]] = ju ;
				ju = iu ;
			} else {
				break ;
			}
		}
		
		/* on insere reellement l'element a sa place */
		dhp->data[ju] = (void *) element ;
		
		/* il reste a le mettre a sa place dans la file down */
		jd = dhp->size - 1 ;
		
		while(jd > 0) {
			uint32_t id = (jd - 1) >> 1 ;

			if(dhp->comp(element,dhp->data[dhp->downToUp[id]],ptUser) < 0 ) {
			
				/* l'element est meilleur que le pere pour le classement
				down, on permute */
				dhp->downToUp[jd] = dhp->downToUp[id] ;
				dhp->upToDown[dhp->downToUp[jd]] = jd ;
				jd = id ;
			} else {
				break ;
			}
		}
		
		/* on insere l'element dans le classement down */
		dhp->downToUp[jd] = ju ;
		dhp->upToDown[ju] = jd ;
	
	} else {
	
		/* l'emplacement en tete est libre, on insere a partir de la tete */
		dhp->isTopFirst = 1 ;
		dhp->size++ ;
		DHEAP_InsertFromTop(dhp,element) ;
		
		/* on insere dans le classement down */
		DHEAP_InsertFromDown(dhp);
	}
	return NULL ;
}
/* -----------------------------------------------------------------------------
 * Fonction	: HEAP_Free
 *
 * Role		:
 *	Liberation de la pile
 *
 * Retour	:
 *	- NULL
 * -----------------------------------------------------------------------------
 */
HEAP_DATA * HEAP_Free(HEAP_DATA *hp)
{
	if(hp != NULL) {
		if(hp->data != NULL) free(hp->data);
		free(hp) ;
	}
	return NULL ;
}

/* -----------------------------------------------------------------------------
 * Fonction	: DHEAP_Free
 *
 * Role		:
 *	Liberation de la pile
 *
 * Retour	:
 *	- NULL
 * -----------------------------------------------------------------------------
 */
DHEAP_DATA * DHEAP_Free(DHEAP_DATA *dhp)
{
	if(dhp != NULL) {
		if(dhp->data != NULL) free(dhp->data);
		if(dhp->upToDown != NULL) free(dhp->upToDown);
		if(dhp->downToUp != NULL) free(dhp->downToUp);
		free(dhp) ;
	}
	return NULL ;
}


#ifdef HEAP_TEST

/* -----------------------------------------------------------------------------
 * Programme de test qui va inserer sizHeap elements dans une pile
 * simple Heap de taille sizHeap, et inserer les memes elements dans une
 * pile double DHeap de taille sizDHeap.
 * Par defaut sizHeap = 10000 et sizDHeap = 5000. Ces valeurs peuvent
 * etre changees en argument du programme, sizDHeap comme 1er argument
 * sizHeap comme 2ieme argument.
 * L'insertion des elements dans Heap se fait en inserant 3 elements puis
 * en retirant le meilleur element et ainsi de suite. Chaque element retire
 * est reinsere apres epuisement des autres elements.
 * L'insertion dans DHeap se fait de la meme facon, mais de plus les elements
 * en debordement (pile pleine) sont stockes dans une queue secondaire qui
 * sera inseree quand la pile sera videe.
 * On retire ensuite les lements de chaque pile Heap et Dheap en verifiant
 * que les valeurs sont identiques. Si la pile Dheap est vide, on reinsere
 * tous les elements de la queue secondaire.
 * -----------------------------------------------------------------------------
 */

typedef struct element ELEMENT ;

struct element {
	int32_t	val;
	ELEMENT *nxt;
};

/* routine de comparaison pour heap */
static int32_t nbComp = 0L ;

int comp(const void *el1,const void *el2)
{
	int32_t v1 = ((ELEMENT *) el1)->val ;
    int32_t v2 = ((ELEMENT *) el2)->val ;
    nbComp++;
    return v1 - v2 ;
}



/* routine de comparaison pour Dheap */
static int32_t nbCompD = 0L ;

int compD(const void *el1,const void *el2)
{
	int32_t v1 = ((ELEMENT *) el1)->val ;
    int32_t v2 = ((ELEMENT *) el2)->val ;
	nbCompD++;
    return v1 - v2 ;
}


/* routines locales */
static void InsertInHeap(ELEMENT *ptElIn, HEAP_DATA *hp) ;
static ELEMENT *InsertInDHeap(ELEMENT *ptElIn, DHEAP_DATA *dhp) ;



/* -----------------------------------------------------------------------------
 * Fonction	: HEAP_TEST
 *
 * Role		:
 *	Programme de test
 *
 * Retour	:
 *	- neant
 * -----------------------------------------------------------------------------
 */
int HEAP_Test(int argc,char **argv)
/* argv[1] = sizDHeap;  argv[2] = sizHeap;  */
{
	ELEMENT		*tabNonTrie = NULL;
	ELEMENT		*nxtSerie	= NULL;
	DHEAP_DATA	*dhp		= NULL;
	HEAP_DATA	*hp			= NULL;

	static int32_t sizHeap  = 10000;
	static int32_t sizDHeap = 5000;

	int32_t i,nbExtr ;

	/* prise en compte des arguments */
	if(argc > 1) {
		sizDHeap = atoi(argv[1]) ;
		if(argc > 2) {
			sizHeap = atoi(argv[2]) ;
		}
	}

	/* on alloue le tableau des elements non tries ainsi que les 2 piles */
	tabNonTrie = (ELEMENT *) calloc(sizHeap, sizeof(tabNonTrie[0])) ;
	if(tabNonTrie == NULL){
		printf("Manque de place pour tabNonTrie\n");
		exit(0) ;
	}

	dhp = DHEAP_Alloc(sizDHeap ,compD) ;
	if(dhp == NULL) {
		printf("Manque de place pour DHEAP_Alloc\n");
		free( tabNonTrie );
		exit(0) ;
	}

	hp  = HEAP_Alloc(sizHeap ,comp) ;
	if(hp == NULL) {
		printf("Manque de place pour HEAP_Alloc\n");
		free( tabNonTrie );
		DHEAP_Free( dhp );
		exit(0) ;
	}

    srand(1);
	/* on initialise les valeurs du tableau non trie */
	for(i=0;i<sizHeap;i++) {
		tabNonTrie[i].val = ( rand() * 32768L + rand())  & 0x7fffffff  ;
	}
	
	/* on chaine les elements en boucle pour preparer une insertion */
	for(i=0;i<sizHeap-1;i++) {
		tabNonTrie[i].nxt = tabNonTrie + i + 1 ;
	}
	
	tabNonTrie[sizHeap-1].nxt = tabNonTrie ;   /* boucle */
	
	/* on insere dans la pile simple */
	InsertInHeap(tabNonTrie,hp);
	
	/* on chaine les elements en boucle pour preparer une insertion */
	for(i=0;i<sizHeap-1;i++) {
		tabNonTrie[i].nxt = tabNonTrie + i + 1 ;
	}
	
	tabNonTrie[sizHeap-1].nxt = tabNonTrie ;   /* boucle */

	/* on insere dans la pile double en recuperant la chaine des elements
	en debordement (si sizDHeap < sizHeap ). */
	nxtSerie = InsertInDHeap(tabNonTrie,dhp);

	/* on va maintenant vider les piles en comparant */
	{
		int32_t valAnt = 0 ;

		for(nbExtr = 0 ; hp->size > 0 ; nbExtr++ ) {
		
			/* on extrait un element de chaque pile */
			ELEMENT *el  = (ELEMENT *) HEAP_GetFirst(hp);
			ELEMENT *del  = (ELEMENT *) DHEAP_GetFirst(dhp);
			if(del == NULL) {
				/* dhp vide, on insere les elements en debordement contenus
				   dans nxtSerie */
				nxtSerie = InsertInDHeap(nxtSerie,dhp);
				del = (ELEMENT *) DHEAP_GetFirst(dhp) ;
				if(del == NULL) printf(" Dheap anormalement vide\n");
			}
			
			/* comparaison de la valeur des 2 elements. Remarque on ne
			peut comparer les 2 elements car s'il existe 2 elements ayant
			la meme la valeur il ne sont pas necessairement classes dans
			le meme ordre dans les 2 piles */
			if(del->val != el->val ) {
				printf(" erreur : ( HP[%d] = %x ) !=  ( DHP[%d] = %x )\n"
				,nbExtr,el->val,nbExtr,del->val );
			} ;
			
			/* on verifie que l'ordre est bien decroissant */
			if((nbExtr > 0) && (el->val > valAnt)) {
				printf(" erreur : ( HP[%d] = %x ) > ( HP[%d] = %x )\n"
				,nbExtr,el->val,nbExtr-1,valAnt );
			}
			valAnt = el->val ;
		}
	}
	
	free( tabNonTrie );
	DHEAP_Free(dhp);
	HEAP_Free(hp);
		
	printf(" %d comparaison pour extraire %d el dans Heap\n",
					nbComp,nbExtr);
	printf(" %d comparaison pour extraire %d el dans DHeap\n",
					nbCompD,nbExtr);
	return (0);
}



/*  Insere les elements contenus dans la boucle ptElIn dans la pile hp
	Un fois sur 4 on retire le meilleur element que l'on reinsere dans la
	boucle ptElIn
*/
static void InsertInHeap(ELEMENT *ptElIn, HEAP_DATA *hp)
{
	int32_t j = 0 ;
	int32_t nbInser = 0 ;
	int32_t oldNbComp = nbComp ;

	while(ptElIn != NULL) {
		if( j++ == 4 ) {
			ELEMENT *el = (ELEMENT *) HEAP_GetFirst(hp) ;
			j = 0 ;
			/* on remet l'element retire dans la boucle */
			if(el != NULL) {
				el->nxt = ptElIn->nxt  ;
				ptElIn->nxt = el ;
				ptElIn = el ;
			}
		} else {
			/* on insere le prochain element */
			ELEMENT *el = ptElIn->nxt ;
			nbInser++ ;
			if(el != ptElIn) {
				ptElIn->nxt = el->nxt ;
			} else {
				ptElIn = NULL ;
			}
            if(HEAP_Insert(hp,el) != 0) { printf("ERR insert %u element\n",el->val); }
		}
	}
	printf(" insertion ds HP de %d elements -> %.1f comp pour %d insertions\n",
		hp->size,( nbComp - oldNbComp) / (float) nbInser ,nbInser );
	nbComp = oldNbComp ;
	return ;
}

/*  insere les elements contenus dans la boucle ptElIn dans la pile dhp
	Un fois sur 4 on retire le meilleur element que l'on reinsere dans la
	boucle ptElIn
	Les elements en debordement de pile (pile pleine) sont retournees dans
	une boucle ptElOut
*/

static ELEMENT *InsertInDHeap(ELEMENT *ptElIn, DHEAP_DATA *dhp)
{
	int32_t	 j = 0 ;
	ELEMENT *ptElOut = NULL ;
	int32_t nbInser = 0 ;
	int32_t oldNbCompD = nbCompD ;

	while(ptElIn != NULL) {
		if( j++ == 4 ) {
			ELEMENT *el = (ELEMENT *) DHEAP_GetFirst(dhp) ;
			j = 0 ;
			/* remet l'element retire dans la boucle a inserer */
			if(el != NULL) {
				el->nxt = ptElIn->nxt  ;
				ptElIn->nxt = el ;
				ptElIn = el ;
			}
		} else {
			ELEMENT *el = ptElIn->nxt ;
			if(el != ptElIn) {
				ptElIn->nxt = el->nxt ;
			} else {
				ptElIn = NULL ;
			}
			nbInser++ ;
			/* on insere le prochain element */
			if( ( el = (ELEMENT *) DHEAP_Insert(dhp,el) ) != NULL) {
				/* on recupere un element en debordement */
				if(ptElOut == NULL) {
					ptElOut = el ;
					el->nxt = el ;
				} else {
					el->nxt = ptElOut->nxt ;
					ptElOut->nxt = el ;
				}
			}
		}
	}

	printf(" insertion ds DHP de %d elements -> %.1f comp pour %d insertions\n",
		dhp->size,( nbCompD - oldNbCompD) / (float) nbInser ,nbInser );
	nbCompD = oldNbCompD ;
	return ptElOut ;
}




#endif	/* HEAP_TEST */

/********************** version HEAP_INDEX ******************************************/


struct HEAP_INDEX {
	uint32_t        *   index;			/* index des elements 	        */
	uint32_t            sizeMax;	    /* taille maximum				*/
	uint32_t	        size;			/* taille effective				*/
	int32_t        		isTopFirst;		/* Premier element present	?	*/
	int				(*comp) (uint32_t ind1,uint32_t ind2, void * ptUser) ;
	int				(*compSimple) (uint32_t ind1,uint32_t ind2) ;
	void * ptUser;
										/* Routine de comparaison		*/
} ;

struct DHEAP_INDEX {
	uint32_t        *   index;			/* index des elements 	        */
	uint32_t	        sizeMax;        /* taille maximum				*/
	uint32_t	        size;			/* taille effective				*/
	int32_t        		isTopFirst;		/* Premier element present	?	*/
	int        	(*comp) (uint32_t ind1,uint32_t ind2, void * ptUser) ;
										/* Routine de comparaison		*/
	int				(*compSimple) (uint32_t ind1,uint32_t ind2);
	uint32_t	        freeDown;		/* indice de l'el. libre decr.	*/
	uint32_t	        *upToDown;		/* Corresp. croiss. decroiss.	*/
	uint32_t	        *downToUp;		/* Corresp. reciproque			*/
	void * ptUser;
} ;

static void HEAPIND_InsertFromTop(HEAP_INDEX *hp,H_INDEX ind) ;
static void DHEAPIND_InsertFromTop(DHEAP_INDEX *dhp,H_INDEX ind) ;
static void DHEAPIND_InsertFromDown(DHEAP_INDEX *dhp) ;
static int HEAPIND_LocalComp(H_INDEX ind1, H_INDEX ind2, void * ptUser);
static int DHEAPIND_LocalComp(H_INDEX ind1, H_INDEX ind2, void * ptUser);

static void HEAPIND_InsertFromTop(HEAP_INDEX *hp,H_INDEX indEl) {
	uint32_t i=0 , j = 1 ;
	void * ptUser;
	/* i est la position courante, j est le premier fils de i */
	if (hp->compSimple == NULL) {
		ptUser = hp->ptUser;
	} else {
		ptUser = (void *) hp;
	}
	while(j < hp->size) {
		/* on regarde le plus grand des 2 fils */
		if( ( j+1 <  hp->size ) &&
			(hp->comp (hp->index[j],hp->index[j+1], ptUser) < 0) ) {
			j++ ;
		}
		/* on regarde si l'element est plus grand que les fils */
		if(hp->comp(indEl,hp->index[j] ,ptUser) < 0 ) {
			/* non , on continue la descente en faisant remonter le plus
			grand des fils  */
			hp->index[i] = hp->index[j] ;
			i = j ;
			j += j + 1 ;
		} else {
			/* oui on arrete la */
			break ;
		}
	}
	/* on insere reellement l'element a sa place */
	hp->index[i] = indEl ;
	return ;
}

static void DHEAPIND_InsertFromTop(DHEAP_INDEX *dhp,H_INDEX indEl) 
{
	uint32_t iu = 0 , ju = 1 ;
	void * ptUser;
	/* iu est la position courante, ju est le premier fils de iu dans le
	classement up. iud est la position initiale dans le classement down */
	uint32_t iud = dhp->upToDown[iu] ;
	
	if (dhp->compSimple == NULL){
		ptUser = dhp->ptUser;
	} else {
		ptUser = (void *) dhp;
	}
	while(ju < dhp->size) {
	
		/* on regarde le plus grand des 2 fils */
		if( ( ju + 1 <  dhp->size ) &&
			(dhp->comp (dhp->index[ju],dhp->index[ju + 1] , ptUser) < 0) ) {
			ju++ ;
		}
		
		/* on regarde si l'element est plus grand que les fils */
		if(dhp->comp(indEl,dhp->index[ju],ptUser) < 0 ) {
			/* non , on continue la descente en faisant remonter le plus
			grand des de fils  */
			dhp->index[iu] = dhp->index[ju] ;
			/* mise a jour de la correspondance up <-> down */
			dhp->upToDown[iu] = dhp->upToDown[ju] ;
			dhp->downToUp[dhp->upToDown[iu]] = iu ;
			iu = ju ;
			ju += ju + 1 ;
		} else {
			/* oui on arrete la */
			break ;
		}
	}
	
	/* on insere reellement l'element */
	dhp->index[iu ] = indEl ;
	
	/* mise a jour de la correspondance up <-> down */
	dhp->upToDown[iu] = iud ;
	dhp->downToUp[iud] = iu ;
	return ;
}
static void DHEAPIND_InsertFromDown(DHEAP_INDEX *dhp) {
	uint32_t id = dhp->freeDown ;
	uint32_t idu = dhp->downToUp[id] ;
	uint32_t jd = ( id -1 ) >> 1 ;
	void * ptUser;
	if (dhp->compSimple == NULL){
		ptUser = dhp->ptUser;
	} else {
		ptUser = dhp;
	}
	/* id position courante , jd position du pere , idu position
	initiale dans le classement up */
	if(  (id > 0) &&
        (dhp->comp(dhp->index[idu],dhp->index[dhp->downToUp[jd]], ptUser) < 0) ) {
		
		/* il faut faire remonter l'lement */
		do {
			/* on remonte l'element tant qu'il est inferieur au pere */
			dhp->downToUp[id] = dhp->downToUp[jd] ;
			dhp->upToDown[dhp->downToUp[id]] = id ;
			id = jd ;
			jd = (id -1) >> 1 ;
		} while ((id > 0) &&
				(dhp->comp(dhp->index[idu],dhp->index[dhp->downToUp[jd]], ptUser) < 0) );
		
		/* on insere l'element a sa place finale */
		dhp->downToUp[id] = idu ;
		dhp->upToDown[idu] = id ;
	
	} else {
		/* il faut faire descendre l'element, jd est le premier fils */
		jd = ( id << 1 ) + 1 ;
	
		while(jd < dhp->size) {
		
			/* on regarde le plus petit des 2 fils */
			if( ( jd + 1 <  dhp->size ) &&
				(dhp->comp (dhp->index[dhp->downToUp[jd]],
						   dhp->index[dhp->downToUp[jd + 1]],ptUser) > 0) ) {
				jd++ ;
			}
			
			/* on regarde si l'element est plus petit que les fils */
            if(dhp->comp(dhp->index[idu],dhp->index[dhp->downToUp[jd]],ptUser) > 0 ) {
				/* non , on contine la descente */
				dhp->downToUp[id] = dhp->downToUp[jd] ;
				dhp->upToDown[dhp->downToUp[id]] = id ;
				id = jd ;
				jd += jd + 1 ;
			} else {
				/* oui on arrete la */
				break ;
			}
		}
		
		/* on insere l'element a sa place finale */
		dhp->downToUp[id] = idu ;
		dhp->upToDown[idu] = id ;
	}
	return ;
}

static int HEAPIND_LocalComp(H_INDEX ind1, H_INDEX ind2, void * ptUser)
{
	HEAP_INDEX * hp = (HEAP_INDEX *) ptUser;
	return (hp->compSimple(ind1,ind2));
}

static int DHEAPIND_LocalComp(H_INDEX ind1, H_INDEX ind2, void * ptUser)
{
	DHEAP_INDEX * dhp = (DHEAP_INDEX *) ptUser;
	return (dhp->compSimple(ind1,ind2));
}
HEAP_INDEX * HEAPIND_Alloc(
	const uint32_t sizeMax,
	int (*comp) (H_INDEX ind1,H_INDEX ind2) )
{
	HEAP_INDEX *newHeap = (HEAP_INDEX *) calloc(1, sizeof(HEAP_INDEX));

	if(newHeap == NULL) return newHeap ;
	
	newHeap->index = (H_INDEX *) malloc(sizeMax*sizeof(newHeap->index[0]));
	
	if(newHeap->index == NULL) {
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->sizeMax = sizeMax ;
	newHeap->comp = HEAPIND_LocalComp ;
	newHeap->compSimple = comp;
	newHeap->isTopFirst = 0 ;

	return newHeap ;
}

HEAP_INDEX * HEAPIND_AllocU(
	const uint32_t sizeMax,
	int (*comp) (H_INDEX ind1,H_INDEX ind2, void * ptUser),
	void * ptUser)
{
	HEAP_INDEX *newHeap = (HEAP_INDEX *) calloc(1, sizeof(HEAP_INDEX));

	if(newHeap == NULL) return newHeap ;
	
	newHeap->index = (H_INDEX *) malloc(sizeMax*sizeof(newHeap->index[0]));
	
	if(newHeap->index == NULL) {
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->sizeMax = sizeMax ;
	newHeap->comp = comp ;
	newHeap->compSimple = NULL;
	newHeap->isTopFirst = 0 ;
	newHeap->ptUser = ptUser;

	return newHeap ;
}

HEAP_INDEX * HEAPIND_Realloc(
	HEAP_INDEX *hp,			/* Structure a desallouer	*/
	const int32_t newSize)	/* Nouvelle taille			*/
{
	if(hp == NULL) return hp;
	{
		void * newIndex = realloc(hp->index,sizeof(hp->index[0]) * newSize) ;

		if(newIndex == NULL) {
			free(hp->index);
			free(hp) ;
			return NULL ;
		}
		
		hp->index = (H_INDEX *) newIndex ;
		hp->sizeMax = newSize ;
		return hp ;
	}
}

DHEAP_INDEX * DHEAPIND_Alloc(
	const uint32_t sizeMax,
	int (*comp) (H_INDEX ind1,H_INDEX ind2) )
{
	DHEAP_INDEX *newHeap = (DHEAP_INDEX *) calloc(1, sizeof(DHEAP_INDEX));

	if(newHeap == NULL) return newHeap ;
	
    newHeap->index = (H_INDEX *) malloc(sizeMax*sizeof(newHeap->index[0]));
	
	if(newHeap->index == NULL) {
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->upToDown = (uint32_t *) calloc(sizeMax, sizeof(newHeap->upToDown[0])) ;
	
	if(newHeap->upToDown == NULL) {
        free(newHeap->index);
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->downToUp = (uint32_t *) calloc(sizeMax, sizeof(newHeap->downToUp[0])) ;
	
	if(newHeap->downToUp == NULL) {
		free(newHeap->upToDown) ;
        free(newHeap->index);
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->sizeMax = sizeMax ;
	newHeap->comp = DHEAPIND_LocalComp ;
	newHeap->compSimple = comp;
	newHeap->isTopFirst = 0 ;
	newHeap->freeDown = 0 ;
	return newHeap ;
}
DHEAP_INDEX * DHEAPIND_AllocU(
	const uint32_t sizeMax,
	int (*comp) (H_INDEX ind1,H_INDEX ind2, void * ptUser),
	void * ptUser)
{
	DHEAP_INDEX *newHeap = (DHEAP_INDEX *) calloc(1, sizeof(DHEAP_INDEX));

	if(newHeap == NULL) return newHeap ;
	
    newHeap->index = (H_INDEX *) malloc(sizeMax*sizeof(newHeap->index[0]));
	
	if(newHeap->index == NULL) {
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->upToDown = (uint32_t *) calloc(sizeMax, sizeof(newHeap->upToDown[0])) ;
	
	if(newHeap->upToDown == NULL) {
		free(newHeap->index);
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->downToUp = (uint32_t *) calloc(sizeMax, sizeof(newHeap->downToUp[0])) ;
	
	if(newHeap->downToUp == NULL) {
		free(newHeap->upToDown) ;
		free(newHeap->index);
		free(newHeap);
		return newHeap = NULL ;
	} ;
	
	newHeap->sizeMax = sizeMax ;
	newHeap->comp = comp ;
	newHeap->compSimple = NULL;
	newHeap->isTopFirst = 0 ;
	newHeap->freeDown = 0 ;
	newHeap->ptUser = ptUser;
	return newHeap ;
}

DHEAP_INDEX * DHEAPIND_Realloc(
	DHEAP_INDEX *dhp,		/* Structure a desallouer	*/
	const int32_t newSize)	/* Nouvelle taille			*/
{
	void *newIndex ;

	if(dhp == NULL) return NULL ;
	
    newIndex = realloc(dhp->index,sizeof(dhp->index[0]) * newSize);
	
	if(newIndex == NULL) {
		return DHEAPIND_Free(dhp) ;
	} else {
		dhp->index = (H_INDEX *) newIndex ;
	}
	
	newIndex = realloc(dhp->upToDown,sizeof(dhp->upToDown[0])* newSize) ;
	
	if(newIndex == NULL) {
		return DHEAPIND_Free(dhp) ;
	} else {
		dhp->upToDown = (uint32_t *) newIndex ;
	}

	newIndex = realloc(dhp->downToUp,sizeof(dhp->downToUp[0]) * newSize) ;
	
	if(newIndex == NULL) {
		return DHEAPIND_Free(dhp) ;
	} else {
		dhp->downToUp = (uint32_t *) newIndex ;
	}
	
	dhp->sizeMax = newSize ;
	
	return dhp ;
}

H_INDEX HEAPIND_GetFirst(HEAP_INDEX *hp)
{
	if(hp->size == 0) return H_INDEX_NULL ;  /* pile vide */

	if(hp->isTopFirst) {
		hp->isTopFirst = 0 ;   /* le meilleure est en tete on le retourne */
	} else {
		/* on va faire remonter le dernier element pour recupere la tete */
		HEAPIND_InsertFromTop(hp,hp->index[hp->size]) ;
	}
	
	hp->size-- ;  /* mise a jour de la taille de la pile */
	
	return hp->index[0] ;
}
H_INDEX DHEAPIND_GetFirst(DHEAP_INDEX *dhp)
{
	
	if(dhp->size == 0) return H_INDEX_NULL; /* pile vide */

	if(dhp->isTopFirst == 0) {
		/* on va faire remonter le dernier element pour recupere la tete */
		/*puisque l'element de rg dhp->size passe en zero dans le classement
		  up, on met a jour les correspondances entre le classement
		  up et down */
		dhp->upToDown[0] = dhp->upToDown[dhp->size] ;
		dhp->downToUp[dhp->upToDown[0]] = 0 ;
		
		/* on met a jour le classement up en remettant l'element en
		tete a sa place */
		DHEAPIND_InsertFromTop(dhp,dhp->index[dhp->size]) ;
		
		/* on reorganise le classement down en faisant passer l'element
		de rg dhp->size (ds le classement down) au rang dhp->freeDown
		(toujours ds le class. down), puisque il s'agit de l'element libere
		au coup precedent */
		if(dhp->freeDown < dhp->size) {
			/* on met a jour la correspondance entre up et down */
			dhp->downToUp[dhp->freeDown] = dhp->downToUp[dhp->size] ;
			dhp->upToDown[dhp->downToUp[dhp->freeDown]] = dhp->freeDown ;
			/* on reorganise le classement down */
			DHEAPIND_InsertFromDown(dhp);
		}
	} else {
		dhp->isTopFirst = 0 ;
	}
	
	dhp->size-- ;	/* mise a jour de la taille */
	
	/* on stocke l'emplacment libre */
	dhp->freeDown = dhp->upToDown[0] ;
	
	return dhp->index[0] ;
}

int  HEAPIND_Insert(
	HEAP_INDEX *hp,
	H_INDEX indEl)
{
	void * ptUser;
	if(hp->compSimple == NULL) {
		ptUser = hp->ptUser;
	} else {
		ptUser = hp;
	}
	if(hp->isTopFirst) {
		
		uint32_t j = hp->size++ ;

		if(j >= hp->sizeMax) {
			hp->size-- ; /* pile pleine on ne fait rien */
			return -1 ;
		}

		/* on a insere l'element en queue (hp->size), on le fait remonter
		a sa place en permuttant avec ses peres qui lui sont inferieur */
		while(j > 0) {
			/* rang du pere */
			uint32_t i = (j - 1) >> 1 ;
			if(hp->comp(indEl,hp->index[i], ptUser) > 0 ) {
				/* on permute */
				hp->index[j] = hp->index[i] ;
				j = i ;
			} else {
				break ;
			}
		}
		/* on insere reellement l'element a sa place */
		hp->index[j] = indEl ;
	} else { /* la place en tete est libre, on insere a partir de la tete */

		hp->isTopFirst = 1 ;
		hp->size++;
		HEAPIND_InsertFromTop(hp,indEl) ;
	}
    return 0 ;
}

void HEAPIND_Reorg(
	HEAP_INDEX *hp,
	int (*newComp) (H_INDEX ind1,H_INDEX ind2))
{
	uint32_t nbOk ;

	hp->compSimple = newComp ;
	hp->comp = HEAPIND_LocalComp;
	if(hp->size <= 0) return;
	
	if(hp->isTopFirst == 0) {
		hp->index[0] = hp->index[hp->size] ;
		hp->isTopFirst = 1 ;
	}
	
	for(nbOk = 1 ; nbOk < hp->size; nbOk ++) {
		
		uint32_t j	= nbOk ;
		H_INDEX  indEl		= hp->index[j] ;
		
		/* on insere l'element de rang nbOk a sa place en le faisant remonter
		a sa place en permuttant avec ses peres qui lui sont inferieur */
		while(j > 0) {
			/* rang du pere */
			uint32_t i = (j - 1) >> 1 ;

			if(newComp(indEl,hp->index[i]) > 0 ) {
				/* on permute */
				hp->index[j] = hp->index[i] ;
				j = i ;
			} else {
				break ;
			}
		}
		
		/* on insere reellement l'element a sa place */
		hp->index[j] = indEl ;
	}
}


void  HEAPIND_ReorgU(
	HEAP_INDEX *hp,
	int (*newComp) (H_INDEX ind1,H_INDEX ind2, void * ptUser),
	void * ptUser)

{
	uint32_t nbOk ;

	hp->comp = newComp ;
	hp->ptUser = ptUser;
	hp->compSimple = NULL;
	if(hp->size <= 0) return;
	
	if(hp->isTopFirst == 0) {
		hp->index[0] = hp->index[hp->size] ;
		hp->isTopFirst = 1 ;
	}
	
	for(nbOk = 1 ; nbOk < hp->size; nbOk ++) {
		
		uint32_t j	= nbOk ;
		H_INDEX indEl		= hp->index[j] ;
		
		/* on insere l'element de rang nbOk a sa place en le faisant remonter
		a sa place en permuttant avec ses peres qui lui sont inferieur */
		while(j > 0) {
			/* rang du pere */
			uint32_t i = (j - 1) >> 1 ;

			if(newComp(indEl,hp->index[i],ptUser) > 0 ) {
				/* on permute */
				hp->index[j] = hp->index[i] ;
				j = i ;
			} else {
				break ;
			}
		}
		
		/* on insere reellement l'element a sa place */
		hp->index[j] = indEl ;
	}
}


H_INDEX DHEAPIND_Insert(
	DHEAP_INDEX *dhp,
	H_INDEX indEl)
{
	void * ptUser;
	if (dhp->compSimple == NULL) {
		ptUser = dhp->ptUser;
	} else {
		ptUser = dhp;
	}
	if(dhp->isTopFirst) {

		uint32_t ju = dhp->size++ ;
		uint32_t jd ;

		if(ju >= dhp->sizeMax) {
		
			/* la pile es pleine */
			dhp->size-- ;
			
			/* on compare l'element au pire dans la pile */
			if(dhp->comp(indEl,dhp->index[dhp->downToUp[0]],ptUser) < 0 ) {
				return indEl ; /* l'element est le pire */
			} else {
			
				/* on expulse l'element le pire */
				H_INDEX indExpulse = dhp->index[dhp->downToUp[0]] ;
				ju = dhp->downToUp[0] ;
				
				/* on insere l'element a sa place, en le faisant remonter
				pusique le pire element etait forcement une feuille */
				while(ju > 0) {
					uint32_t iu = (ju - 1) >> 1 ;

					if(dhp->comp(indEl,dhp->index[iu],ptUser) > 0 ) {
						
						/* l'element est meilleurs que le pere */
						dhp->index[ju] = dhp->index[iu] ;
					
						/* mise a jour des corresp. up <-> down */
						dhp->upToDown[ju] = dhp->upToDown[iu] ;
						dhp->downToUp[dhp->upToDown[ju]] = ju ;
						ju = iu ;
					} else {
						break ;
					}
				}
				
				/* on insere reellement l'element a sa place */
				dhp->index[ju] = indEl ;
				dhp->downToUp[0] = ju ;
				dhp->upToDown[ju] = 0 ;
				
				/* il faut maintenant l'inserer dans la file down */
				dhp->freeDown = 0 ;
				DHEAPIND_InsertFromDown(dhp);
				return indExpulse ;
			}
		}
		while(ju > 0) {
			
			/* la pile n'est pas pleine, on insere l'element a partir
			du rang dhp->size en le faisant remonter jusqu'a sa place */
			uint32_t iu = (ju - 1) >> 1 ;

			if(dhp->comp(indEl,dhp->index[iu],ptUser) > 0 ) {
			
				/* l'element est meilleur que le pere */
				dhp->index[ju] = dhp->index[iu] ;
				dhp->upToDown[ju] = dhp->upToDown[iu] ;
				dhp->downToUp[dhp->upToDown[ju]] = ju ;
				ju = iu ;
			} else {
				break ;
			}
		}
		
		/* on insere reellement l'element a sa place */
		dhp->index[ju] = indEl ;
		
		/* il reste a le mettre a sa place dans la file down */
		jd = dhp->size - 1 ;
		
		while(jd > 0) {
			uint32_t id = (jd - 1) >> 1 ;

			if(dhp->comp(indEl,dhp->index[dhp->downToUp[id]],ptUser) < 0 ) {
			
				/* l'element est meilleur que le pere pour le classement
				down, on permute */
				dhp->downToUp[jd] = dhp->downToUp[id] ;
				dhp->upToDown[dhp->downToUp[jd]] = jd ;
				jd = id ;
			} else {
				break ;
			}
		}
		
		/* on insere l'element dans le classement down */
		dhp->downToUp[jd] = ju ;
		dhp->upToDown[ju] = jd ;
	
	} else {
	
		/* l'emplacement en tete est libre, on insere a partir de la tete */
		dhp->isTopFirst = 1 ;
		dhp->size++ ;
		DHEAPIND_InsertFromTop(dhp,indEl) ;
		
		/* on insere dans le classement down */
		DHEAPIND_InsertFromDown(dhp);
	}
	return H_INDEX_NULL ;
}

HEAP_INDEX * HEAPIND_Free(HEAP_INDEX *hp)
{
	if(hp != NULL) {
		if(hp->index != NULL) free(hp->index);
		free(hp) ;
	}
	return NULL ;
}

DHEAP_INDEX * DHEAPIND_Free(DHEAP_INDEX *dhp)
{
	if(dhp != NULL) {
		if(dhp->index != NULL) free(dhp->index);
		if(dhp->upToDown != NULL) free(dhp->upToDown);
		if(dhp->downToUp != NULL) free(dhp->downToUp);
		free(dhp) ;
	}
	return NULL ;
}

#ifdef HEAP_TEST

static int32_t nbCompInd = 0L ;

int compInd(H_INDEX ind1,H_INDEX ind2,void * ptUser)
{
	int32_t * val = (int32_t *) ptUser ;
    nbCompInd++;
	return val[ind1] - val[ind2] ;
}

static int32_t nbCompIndD = 0L ;

int compIndD(H_INDEX ind1,H_INDEX ind2,void * ptUser)
{
  	int32_t * val = (int32_t *) ptUser ;
	nbCompIndD++;
	return val[ind1] - val[ind2] ;
}
static void InsertInHeapInd(HEAP_INDEX *hp,H_INDEX maxVal) ;
static int InsertInDHeapInd(DHEAP_INDEX *dhp,H_INDEX maxVal,H_INDEX * outIndex) ;


int HEAPIND_Test(int argc,char **argv)
/* argv[1] = sizDHeap;  argv[2] = sizHeap;  */
{
	int32_t		*valNonTrie = NULL;
	H_INDEX		*outIndex	= NULL;
	DHEAP_INDEX	*dhp		= NULL;
	HEAP_INDEX	*hp			= NULL;

	static int32_t sizHeap  = 10000;
	static int32_t sizDHeap = 5000;
    int32_t nbOut ;

	int32_t i,nbExtr ;

	/* prise en compte des arguments */
	if(argc > 1) {
		sizDHeap = atoi(argv[1]) ;
		if(argc > 2) {
			sizHeap = atoi(argv[2]) ;
		}
	}

	/* on alloue le tableau des elements non tries ainsi que les 2 piles */
	valNonTrie = (int32_t *) malloc((sizHeap+1)* sizeof(valNonTrie[0])) ;
	outIndex = (H_INDEX *) malloc(sizHeap * sizeof(outIndex[0])) ;
    
    if(valNonTrie == NULL){
		printf("Manque de place pour valNonTrie\n");
		exit(0) ;
	}

	dhp = DHEAPIND_AllocU(sizDHeap ,compIndD,valNonTrie) ;
	if(dhp == NULL) {
		printf("Manque de place pour DHEAPIND_Alloc\n");
		free( valNonTrie );
		exit(0) ;
	}

	hp  = HEAPIND_AllocU(sizHeap ,compInd,valNonTrie) ;
	if(hp == NULL) {
		printf("Manque de place pour HEAP_Alloc\n");
		free( valNonTrie );
		DHEAPIND_Free( dhp );
		exit(0) ;
	}

    srand(1);
	/* on initialise les valeurs du tableau non trie */
	for(i=0;i<sizHeap;i++) {
		valNonTrie[i+1] = ( rand() * 32768L + rand())  & 0x7fffffff ;
        outIndex[i] = i+1 ; /* car index == 0 interdit */
	}
	
	/* on insere dans la pile simple */
	InsertInHeapInd(hp,sizHeap);
    nbOut = InsertInDHeapInd(dhp,sizHeap,outIndex) ;

	/* on va maintenant vider les piles en comparant */
	{
		int32_t valAnt = 0 ;
		for(nbExtr = 0 ; hp->size > 0 ; nbExtr++ ) {
		
			/* on extrait un element de chaque pile */
			H_INDEX indEl  = HEAPIND_GetFirst(hp);
			H_INDEX dindEl  = DHEAPIND_GetFirst(dhp);
			if(dindEl == H_INDEX_NULL) {
				/* dhp vide, on insere les elements en debordement contenus
				   dans nxtSerie */
				nbOut = InsertInDHeapInd(dhp,nbOut,outIndex);
				dindEl = DHEAPIND_GetFirst(dhp) ;
				if(dindEl == H_INDEX_NULL) printf(" Dheap anormalement vide\n");
			}
			
			/* comparaison de la valeur des 2 elements. Remarque on ne
			peut comparer les 2 elements car s'il existe 2 elements ayant
			la meme la valeur il ne sont pas necessairement classes dans
			le meme ordre dans les 2 piles */
			if(valNonTrie[indEl] != valNonTrie[dindEl] ) {
				printf(" erreur : ( HP[%d] = %d ) !=  ( DHP[%d] = %d )\n"
				,nbExtr,valNonTrie[indEl],nbExtr,valNonTrie[dindEl] );
			} ;
			
			/* on verifie que l'ordre est bien decroissant */
			if((nbExtr > 0) && (valNonTrie[indEl] > valAnt)) {
				printf(" erreur : ( HP[%d] = %d ) > ( HP[%d] = %d )\n"
				,nbExtr,valNonTrie[indEl],nbExtr-1,valAnt );
			}
			valAnt = valNonTrie[indEl] ;
		}
	}
	
	free( valNonTrie );
	DHEAPIND_Free(dhp);
	HEAPIND_Free(hp);
		
	printf(" %d comparaison pour extraire %d el dans Heap\n",
					nbCompInd,nbExtr);
	printf(" %d comparaison pour extraire %d el dans DHeap\n",
					nbCompIndD,nbExtr);
	return (0);
}




/* attention pour simuler a l'identique le comportement chaine de la version sans ind
la gestion des index est un peu etrange */
static void InsertInHeapInd( HEAP_INDEX *hp,H_INDEX maxVal)
{
	int32_t j = 0 ;
	int32_t nbInser = 0 ;
	int32_t oldNbCompInd = nbCompInd ;
    H_INDEX * nxtIndex ;
    uint32_t curInd = 1 ; /* pour simuler le chainage */ 
    uint32_t lastInd = 1 ;
    nxtIndex = malloc(maxVal * sizeof(nxtIndex[0])) ;
    {
        uint32_t i ;
        for(i=0;i<maxVal;i++) {
            nxtIndex[i] = i+1 ;
        }
    }
	do {
		if( j++ == 4 ) {
			/* on remet l'element retire dans la boucle */
			 nxtIndex[lastInd++] = HEAPIND_GetFirst(hp) ;
			j = 0 ;
		} else {
			/* on insere le prochain element */
            if(HEAPIND_Insert(hp,nxtIndex[curInd]) != 0) { printf("ERR insert %u element\n",nxtIndex[curInd]); }
            curInd++ ;
            if(curInd == maxVal) curInd = 0 ;
            nbInser++ ;
 		}
	} while (curInd != lastInd) ;
    free(nxtIndex) ;
	printf(" insertion ds HP de %d elements -> %.1f comp pour %d insertions\n",
		hp->size,( nbCompInd - oldNbCompInd) / (float) nbInser ,nbInser );
	nbCompInd = oldNbCompInd ;
	return ;
}

static int InsertInDHeapInd(DHEAP_INDEX *dhp,H_INDEX maxVal,H_INDEX * outIndex) {
	int32_t	 j = 0 ;
	int32_t nbInser = 0 ;
	int32_t oldNbCompIndD = nbCompIndD ;

    uint32_t curInd = 1 ;
    uint32_t lastInd = 1 ;
    int32_t nbIndOut = 0 ;
    H_INDEX * nxtIndex ;
    nxtIndex = malloc(maxVal * sizeof(nxtIndex[0])) ;

	do  {
		if( j++ == 4 ) {
			/* on remet l'element retire dans la boucle */
			 outIndex[lastInd++] = DHEAPIND_GetFirst(dhp) ;
			j = 0 ;
		} else {
			/* on insere le prochain element */
            H_INDEX supInd = DHEAPIND_Insert(dhp,outIndex[curInd]) ;
            if(supInd != H_INDEX_NULL) { 
                nxtIndex[nbIndOut++] = supInd ;
            }
            curInd++ ;
            if(curInd == maxVal) curInd = 0 ;
			nbInser++ ;
        }
	} while (curInd != lastInd) ;
    /* on va maintenant recopier en ordre inverse */
    outIndex[0] = nxtIndex[0] ;
    for(j=1;j<nbIndOut;j++) {
        outIndex[j] = nxtIndex[nbIndOut-j] ;
    }
    free(nxtIndex); 
	printf(" insertion ds DHP de %d elements -> %.1f comp pour %d insertions\n",
		dhp->size,( nbCompIndD - oldNbCompIndD) / (float) nbInser ,nbInser );
	nbCompIndD = oldNbCompIndD ;
	return nbIndOut ;
}

#endif	/* HEAP_TEST */
