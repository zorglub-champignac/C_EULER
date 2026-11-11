//
// Created by Jeannot on 09/11/2026.
//

#ifndef C_EULER_CIBLE_H
#define C_EULER_CIBLE_H
#define CIBLE 21

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
#elif CIBLE == 23
#   define CIBLE_NAME  "perm6_W7"
#   define PERM1   6  // permutation de [012..6[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 7
#   define LIST_W2 {{0,6,7},{1,6,7},{2,6,7},{0,1,6},{0,1,7},{0,2,6},{0,2,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 29  // NbTriplets - MAXT - SIZE_W2  56-(20+7)
#elif CIBLE == 24
#   define CIBLE_NAME  "perm6_W10e"
#   define PERM1   6  // permutation de [012..6[
#   define NBPERM1 720 // 6!
#   define MAXT    20 // 6x5x4/6
#   define SIZE_W2 10
#   define LIST_W2 {{0,1,6},{0,1,7},{0,2,6},{0,2,7},{1,2,6},{1,2,7},{0,3,6},{0,3,7},{1,3,6},{1,3,7}}
#   define PERM2   6
#   define G_PERM2 2
#   define NBPERM2 2
#   define nbRemT 26  // NbTriplets - MAXT - SIZE_W2  56-(20+10)
#endif

#endif //C_EULER_CIBLE_H
