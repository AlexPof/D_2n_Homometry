#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
* EnsD2N is the structure definition for subsets of D_2n
*
* Each element in the D_2n dihedral group can be written as (g,h) where
* g is an element of Z_n and h is an element of Z_12.
* We represent a subset of D_2n by its decomposition on the two cosets
* (Z_n,1) and (Z_n,h), which are here represented as tables A0 and A1.
* In addition, we represent a subset of D_2n as a single value valid
* which is the concatenation of two words of length n in {0,1}, one form
* each coset.
*
*/
typedef struct EnsD2N {
  int* A0;
  int* A1;
  unsigned long long int val;
  int N;
} EnsD2N;

void initialize_EnsD2N(EnsD2N* ensemble, unsigned long long int a, int N);
void niceprint_EnsD2N(EnsD2N* ensemble);
int is_left_homometric(EnsD2N* X,EnsD2N* Y);
int is_right_homometric(EnsD2N* X,EnsD2N* Y);
int is_left_translated(EnsD2N* X,EnsD2N* Y);
int is_right_translated(EnsD2N* X,EnsD2N* Y);
int is_pureZn(EnsD2N* X);
int mod(int a, int b);
unsigned long long int power2(unsigned long long int  a);
unsigned long long int next_samebits_number(unsigned long long int x);
void binaryprint(unsigned int a,int N);

/*
* main takes 2 arguments
*
*
* N: order of the D_2N dihedral group
* P: cardinality of the subsets to be considered
*
*/

int main(int argc, char *argv[]) {
  EnsD2N* ensembles;
  int total_ensembles=0;
  EnsD2N X;
  unsigned long long int x;
  int i,j,flag;
  int count;
  int homomtype;

  int N,P;

  if (argc<4) {
    printf("Not enough arguments - Need at least 3\n");
    exit(1);
  } else {
    if (strcmp("left",argv[1])==0)
      homomtype=0;
    else if (strcmp("right",argv[1])==0)
      homomtype=1;
    else {
      printf("Homometry type should be either left or right\n");
      exit(1);
    }
    N = atoi(argv[2]);
    P = atoi(argv[3]);
  }

  ensembles = (EnsD2N*)malloc(sizeof(EnsD2N));

  /*
  * The enumeration proceeds in two parts:
  *
  * First, we determine a collection of D_2n ensembles that are
  * not related by a translate transform of D_2n and which are not ensembles
  * of pure Z_n
  *
  */

  printf("Building the collection of D_%d ensembles of cardinality %d...\n",2*N,P);

  x=power2(P-1)-1;
  while(x<power2(2*N-1)) {
    initialize_EnsD2N(&X, (1+2*x), N);
    if (!is_pureZn(&X)) {
      if (total_ensembles==0) {
        total_ensembles++;
        ensembles = (EnsD2N*)malloc(sizeof(EnsD2N));
        initialize_EnsD2N(ensembles, (1+2*x), N);
      } else {
        flag=0;

        for (i=0;i<total_ensembles;i++) {
          switch(homomtype) {
            case 0:
              flag=flag || is_right_translated(&X,&ensembles[i]);
              break;
            case 1:
              flag=flag || is_left_translated(&X,&ensembles[i]);
              break;
          }
        }
        if (flag==0) {
          total_ensembles++;
          ensembles = (EnsD2N*)realloc(ensembles,total_ensembles*sizeof(EnsD2N));
          initialize_EnsD2N(&ensembles[total_ensembles-1], (1+2*x), N);
        }
      }
    }
    x = next_samebits_number(x);
  }
  printf("Found %d D_%d ensembles of cardinality %d...\n",total_ensembles,2*N,P);

  /*
  * Then we determine the homometric sets by examining
  * all pair-wise combinations
  */

  switch(homomtype) {
    case 0:
      printf("Determining left-homometric D_%d ensembles of cardinality %d...\n",2*N,P);
      break;
    case 1:
      printf("Determining right-homometric D_%d ensembles of cardinality %d...\n",2*N,P);
      break;
  }

  count=0;
  for(i=0;i<total_ensembles;i++) {
    for(j=i+1;j<total_ensembles;j++) {
      flag=0;
      switch(homomtype) {
        case 0:
          flag = is_left_homometric(&ensembles[i],&ensembles[j]);
          break;
        case 1:
          flag = is_right_homometric(&ensembles[i],&ensembles[j]);
          break;
      }
      if (flag) {
        printf("%llu-%llu\n",ensembles[i].val,ensembles[j].val);
        niceprint_EnsD2N(&ensembles[i]);
        niceprint_EnsD2N(&ensembles[j]);
        printf("*******************\n");
        count++;
      }
    }
  }
  exit(0);
}

///////////////////////////////////////////////////////

int is_left_homometric(EnsD2N* X,EnsD2N* Y){
  /*
  * Function:  is_left_homometric
  * --------------------
  * Checks if two D_2n subsets X and Y are left homometric. We check This
  * by calculating the left interval vector of X and Y, and comparing them, since
  * left homometric sets have identical left interval vectors.
  *
  * X,Y: the D_2n subsets to be checked
  *
  *  returns: True if the subsets X and Y are left homometric,
  *            False otherwise.
  */
  int i,j,N;
  int countX,countY;

  N=X->N;
  for(i=0;i<N;i++) {
    countX=0;
    countY=0;

    for(j=0;j<N;j++) {
      countX+= (X->A0[j]) & (X->A0[ mod((j+i),N) ]);
      countX+= (X->A1[j]) & (X->A1[ mod((j+i),N) ]);

      countY+= (Y->A0[j]) & (Y->A0[ mod((j+i),N) ]);
      countY+= (Y->A1[j]) & (Y->A1[ mod((j+i),N) ]);
    }
    if (countX!=countY)
      return 0;
  }

  for(i=0;i<N;i++) {
    countX=0;
    countY=0;

    for(j=0;j<N;j++) {
      countX+= (X->A0[j]) & (X->A1[ mod((i-j),N) ]);

      countY+= (Y->A0[j]) & (Y->A1[ mod((i-j),N) ]);
    }
    if (countX!=countY)
      return 0;
  }


  return 1;
}

int is_right_homometric(EnsD2N* X,EnsD2N* Y){
  /*
  * Function:  is_right_homometric
  * --------------------
  * Checks if two D_2n subsets X and Y are right homometric. We check This
  * by calculating the right interval vector of X and Y, and comparing them,
  * since right homometric sets have identical right interval vectors.
  *
  * X,Y: the D_2n subsets to be checked
  *
  *  returns: True if the subsets X and Y are right homometric,
  *            False otherwise.
  */

  int i,j,N;
  int countX,countY;

  N=X->N;
  for(i=0;i<N;i++) {
    countX=0;
    countY=0;

    for(j=0;j<N;j++) {
      countX+= (X->A0[j]) & (X->A0[ mod((j+i),N) ]);
      countX+= (X->A1[j]) & (X->A1[ mod((j+i),N) ]);

      countY+= (Y->A0[j]) & (Y->A0[ mod((j+i),N) ]);
      countY+= (Y->A1[j]) & (Y->A1[ mod((j+i),N) ]);
    }
    if (countX!=countY)
      return 0;
  }

  for(i=0;i<N;i++) {
    countX=0;
    countY=0;

    for(j=0;j<N;j++) {
      countX+= (X->A0[j]) & (X->A1[ mod((j+i),N) ]);

      countY+= (Y->A0[j]) & (Y->A1[ mod((j+i),N) ]);
    }
    if (countX!=countY)
      return 0;
  }


  return 1;
}

int is_left_translated(EnsD2N* X,EnsD2N* Y){
  /*
  * Function:  is_left_translated
  * --------------------
  * Checks if two D_2n subsets X and Y are left translates of each other, i.e.
  * there exists (p,q) in D_2n such that the set {(p,q)(g,h) with (g,h) in X}
  * is the subset Y.
  *
  * X,Y: the D_2n subsets to be checked
  *
  *  returns: True if the subsets X and Y are left translates of each other,
  *            False otherwise.
  */
  int p,i,N,c,d;

  N=X->N;
  for(p=0;p<N;p++) {
    c=0;
    d=0;
    for(i=0;i<N;i++) {
      if (X->A0[i] == Y->A0[ mod((i+p),N)])
        c++;
      if (X->A1[i] == Y->A1[ mod((i+p),N)])
        c++;

      if (X->A0[i] == Y->A1[ mod((p-i),N)])
        d++;
      if (X->A1[i] == Y->A0[ mod((p-i),N)])
        d++;
    }
    if (c==(2*N) || d==(2*N))
      return 1;
  }

  return 0;
}


int is_right_translated(EnsD2N* X,EnsD2N* Y){
  /*
  * Function:  is_right_translated
  * --------------------
  * Checks if two D_2n subsets X and Y are right translates of each other, i.e.
  * there exists (p,q) in D_2n such that the set {(g,h)(p,q) with (g,h) in X}
  * is the subset Y.
  *
  * X,Y: the D_2n subsets to be checked
  *
  *  returns: True if the subsets X and Y are right translates of each other,
  *            False otherwise.
  */
  int p,i,N,c,d;

  N=X->N;
  for(p=0;p<N;p++) {
    c=0;
    d=0;
    for(i=0;i<N;i++) {
      // This checks the case where (p,q) is such that q is the identity in Z_2
      if (X->A0[i] == Y->A0[ mod((i+p),N)])
        c++;
      if (X->A1[i] == Y->A1[ mod((i-p),N)])
        c++;

      // This checks the case where (p,q) is such that q is non-trivial in Z_2
      if (X->A0[i] == Y->A1[ mod((i+p),N)])
        d++;
      if (X->A1[i] == Y->A0[ mod((i-p),N)])
        d++;
    }
    if (c==(2*N) || d==(2*N))
    return 1;
  }

  return 0;
}

int is_pureZn(EnsD2N* X) {
  /*
  * Function:  is_pureZn
  * --------------------
  * Checks if a D_2n subset is a pure Z_n subset, i.e. all elements of the
  * subset belong to a single coset (Z_n,1) or (Z_n,h).
  *
  * X: the D_2n subset to be checked
  *
  *  returns: True if the subset is a pure Z_n subset, False otherwise
  */
  int i,N,c,d;

  N=X->N;
  c=0;
  d=0;
  for (i=0;i<N;i++) {
    // We count the number of elements of the form (g,1) or (g,h)
    c+=X->A0[i];
    d+=X->A1[i];
  }
  // If one of them is equal to zero, this means that all elements are
  // either of the form (g,1) or (g,h)
  return (c==0 || d==0);
}

void initialize_EnsD2N(EnsD2N* ensemble, unsigned long long int a, int N) {
  /*
  * Function:  initialize_EnsD2N
  * --------------------
  * Initialize a D_2n subset with the given values
  *
  * ensemble: the D_2n subset to be initialized
  * a: the integer representation of the D_2n subset
  * N: the order n of D_2n
  *
  *  returns: None
  */
  int i,j;
  int countX;

  ensemble->N = N;
  ensemble->val = a;
  ensemble->A0 = calloc(N,sizeof(int));
  ensemble->A1 = calloc(N,sizeof(int));

  for(i=0;i<ensemble->N;i++) {
    ensemble->A0[i] = (a>>i) & 1;
    ensemble->A1[i] = (a>>(N+i)) & 1;
  }
}

void niceprint_EnsD2N(EnsD2N* ensemble) {
  /*
  * Function:  niceprint_EnsD2N
  * --------------------
  * Print an interpretable version of a D_2n subset
  *
  * ensemble: the D_2n subset to be printed
  *
  *
  *  returns: None. Prints an element (g,h) of the subset
  *           as g+ if h=1, g- otherwise.
  */
  int i;

  printf("{");
  for(i=0;i<ensemble->N;i++) {
    if (ensemble->A0[i])
      printf("%d+,",i);
  }
  for(i=0;i<ensemble->N;i++) {
    if (ensemble->A1[i])
      printf("%d-,",i);
  }
    printf("}\n");
}

unsigned long long int next_samebits_number(unsigned long long int x) {
  /*
  * Function:  next_samebits_number
  * --------------------
  * Given a set of n elements,
  * produce a new set of n elements. If you start with the
  * result of next_samebits_number(k)/ and then at each
  * step apply next_samebits_number to the previous result,
  * and keep going until a set is obtained containing m as a
  * member, you will have obtained a set representing all
  * possible ways of choosing k things from m things.
  */
  unsigned long long int  smallest, ripple, new_smallest, ones;

  smallest     = (x & -x);
  ripple       = x + smallest;
  new_smallest = (ripple & -ripple);
  ones         = ((new_smallest/smallest) >> 1) - 1;
  return ripple | ones;
}


int mod(int a, int b) {
  /*
  * Function:  mod
  * --------------------
  * Arithmetic modulo of a by b
  *
  *  returns: the value of mod(a,b).
  */
  return (a >= 0 ? a % b : b - (-a) % b);
}

unsigned long long int power2(unsigned long long int  a) {
  unsigned long long int v=1;
  return (v << a);
}

void binaryprint(unsigned int a,int N) {
  /*
  * Function:  binaryprint
  * --------------------
  * Print the binary representation of an integer a considered as a word
  * of length N. For debug purposes.
  */
  int i;
  for(i=0;i<N;i++)
  printf("%d",(a>>i)&1);
  printf("\n");
}
