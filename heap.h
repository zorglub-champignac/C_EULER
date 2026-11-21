#ifndef _HEAP_H_
#define _HEAP_H_

/* -----------------------------------------------------------------------------
 * Structure des piles et tris .....
 *
 * La structure HEAP_DATA permet de gerer une pile d'elements ou le meilleur
 * element est toujours sur le sommet de la pile. On peut a tout moment retirer
 * le meilleur element de la pile (HeapGetFirst) ou inserer un nouvel element
 * (HeapInsert)
 * Attention : si la pile est pleine l'insertion du nouvel element n'est pas
 * faite.
 * La structure DHEAP_DATA offre les meme fonctionalites mais en plus si la
 * pile est pleine DHeapInsert insere le nouvel element et retourne le plus
 * mauvais contenu dans la pile, ce qui garantit que la pile contient toujours
 * les meilleurs elements.
 * Remarque : les elements ne sont pas physiquement inseres dans la pile,
 * on garde un pointeur sur eux. Le programme appelant proprietaire des
 * elements ne doit donc pas les modifier a l'insu des routines de gestion
 * des piles.
 * La structure HEAP_INDEX a exactement le meme fonctionnet, mails les donnees inserees
 * sont gerees par un index (> 0) au lieu d'un pointeur. C'est a l'appelant de gerer lui meme 
 * les donnees
  * -----------------------------------------------------------------------------
 */

/* ----------------------
 * Directives d'inclusion
 * ----------------------
 */
#include <stdint.h>

/* ------------------------------------
 * Declaration des Constantes exportees
 * ------------------------------------
 */
 
 /* ---------------------------------------------
 * Declaration des Types et Structures exportees
 * ---------------------------------------------
 */
typedef struct HEAP_DATA HEAP_DATA ;


typedef struct DHEAP_DATA DHEAP_DATA ;


/* ----------------------------------
 * Prototypes des Fonctions exportees
 * ----------------------------------
 */

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------------------
 * Fonction	: HEAP_Alloc et DHEAP_Alloc
 *
 * Role		:
 *	Allocation d'une pile
 *
 * Retour	:
 *	-
 * -----------------------------------------------------------------------------
 */
 HEAP_DATA * HEAP_Alloc(
	const uint32_t uneSizeMax, /* Taille demandee			*/
	int (*pUneFctComp) (const void *el1,const void *el2));
									/* Routine de comparaison	*/
 HEAP_DATA *  HEAP_AllocU(
	const uint32_t uneSizeMax, /* Taille demandee			*/
	int (*pUneFctComp) (const void *el1,const void *el2, void * ptUser),
	void * ptUser);
									/* Routine de comparaison	*/
 DHEAP_DATA *  DHEAP_Alloc(
	const uint32_t uneSizeMax,
	int (*pUneFctComp) (const void *el1,const void *el2));

DHEAP_DATA *  DHEAP_AllocU(
	const uint32_t sizeMax,
	int (*comp) (const void *el1,const void *el2, void * ptUser),
	void * ptUser);									

/* -----------------------------------------------------------------------------
 * Fonction	: HEAP_GetFirst et DHEAP_GetFirst
 *
 * Role		:
 *	Retourne le meilleur element
 *
 * Retour	:
 *	-
 * -----------------------------------------------------------------------------
 */
 void *  HEAP_GetFirst(HEAP_DATA *pUneHeap);

 void *  DHEAP_GetFirst(DHEAP_DATA *pUneDHeap);


/* -----------------------------------------------------------------------------
 * Fonction	: HEAP_Insert et DHEAP_Insert
 *
 * Role		:
 *	Insere un element
 *
 * Retour	:
 *	- cas HEAP_Insert  : neant
 *	- cas DHEAP_Insert :
 * -----------------------------------------------------------------------------
 */
 int  HEAP_Insert(
	HEAP_DATA	*pUneHeap,		/* La pile				*/
	const void	*pUnElement);   /* L'element a inserer	*/

 void *  DHEAP_Insert(
	DHEAP_DATA	*pUneDHeap,
	void		*pUnElement);


/* -----------------------------------------------------------------------------
 * Fonction	: HEAP_Reorg
 *
 * Role		:
 *	Reordonnancement d'une pile
 *
 * Retour	:
 *	- neant
 * -----------------------------------------------------------------------------
 */
 void  HEAP_Reorg(
	HEAP_DATA	*pUneHeap,		/* La pile					*/
	int (*pUneFctComp) (const void *el1,const void *el2));
								/* routine de comparaison	*/

 void  HEAP_ReorgU(
	HEAP_DATA	*pUneHeap,		/* La pile					*/
	int (*pUneFctComp) (const void *el1,const void *el2, void * ptUser),
	void * ptUser);

/* -----------------------------------------------------------------------------
 * Fonction	: HEAP_Free et DHEAP_Free
 *
 * Role		:
 *	Liberation de la pile
 *
 * Retour	:
 *	-
 * -----------------------------------------------------------------------------
 */
 HEAP_DATA *  HEAP_Free(HEAP_DATA *pUneHeap);

 DHEAP_DATA *  DHEAP_Free(DHEAP_DATA *pUneDHeap);


/* -----------------------------------------------------------------------------
 * Fonction	: HEAP_Realloc et DHEAP_Realloc
 *
 * Role		:
 *	Reallocation de la pile
 *
 * Retour	:
 *	-
 * -----------------------------------------------------------------------------
 */
 HEAP_DATA *  HEAP_Realloc(
	HEAP_DATA *pUneHeap,
	const int32_t uneNewSize);

 DHEAP_DATA *  DHEAP_Realloc(
	DHEAP_DATA *pUneDHeap,
	const int32_t uneNewSize);


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
#ifdef	HEAP_TEST
 int  HEAP_Test(
    int		argc,
	char	*argv[]);
#endif


#ifdef __cplusplus
}
#endif



typedef struct HEAP_INDEX HEAP_INDEX ;
typedef  uint32_t H_INDEX ;
#define H_INDEX_NULL  ((uint32_t) 0 )
typedef struct DHEAP_INDEX DHEAP_INDEX ;

#ifdef __cplusplus
extern "C" {
#endif

HEAP_INDEX * HEAPIND_Alloc(
	const uint32_t uneSizeMax, /* Taille demandee			*/
	int (*pUneFctComp) (H_INDEX ind1,H_INDEX ind2));
									/* Routine de comparaison	*/
 HEAP_INDEX *  HEAPIND_AllocU(
	const uint32_t uneSizeMax, /* Taille demandee			*/
	int (*pUneFctComp) (H_INDEX ind1,H_INDEX ind2, void * ptUser),
	void * ptUser);
									/* Routine de comparaison	*/
 DHEAP_INDEX *  DHEAPIND_Alloc(
	const uint32_t uneSizeMax,
	int (*pUneFctComp) (H_INDEX ind1,H_INDEX ind2));

DHEAP_INDEX *  DHEAPIND_AllocU(
	const uint32_t sizeMax,
	int (*comp) (H_INDEX ind1,H_INDEX ind2, void * ptUser),
	void * ptUser);									

 H_INDEX  HEAPIND_GetFirst(HEAP_INDEX *pUneHeapInd);

 H_INDEX  DHEAPIND_GetFirst(DHEAP_INDEX *pUneDHeapInd);

 int  HEAPIND_Insert(
	HEAP_INDEX	*pUneHeap,		/* La pile				*/
    H_INDEX   ind );   /* Lindex de l'element a inserer	*/

 H_INDEX  DHEAPIND_Insert(
	DHEAP_INDEX	*pUneDHeap,
	H_INDEX   ind);

  HEAP_INDEX *  HEAPIND_Free(HEAP_INDEX *pUneHeap);

 DHEAP_INDEX *  DHEAPIND_Free(DHEAP_INDEX *pUneDHeap);
 HEAP_INDEX *  HEAPIND_Realloc(
	HEAP_INDEX *pUneHeap,
	const int32_t uneNewSize);

 DHEAP_INDEX *  DHEAPIND_Realloc(
	DHEAP_INDEX *pUneDHeap,
	const int32_t uneNewSize);

#ifdef	HEAP_TEST
 int HEAPIND_Test(int argc,char **argv);
#endif


#ifdef __cplusplus
}
#endif

#endif	/* _HEAP_H_ */
