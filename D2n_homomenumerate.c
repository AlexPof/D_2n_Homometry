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
* More precisely, if D_2n is the group given by the presentation
* D_2n = < s,t | t^n=1=s^2, sts=t^-1>
* then A0 represents the coset of elements of the form t^p, p=0...n
* and A1 represents the coset of elements of the form (t^p)s, p=0...n
*/
typedef struct EnsD2N {
  int* A0;
  int* A1;
  char* hash;
  int N;
} EnsD2N;

void enumerate_homometric(int N, int P, FILE* output_file);
void initialize_EnsD2N(EnsD2N* ensemble, int* sequence, int N);
void get_left_IV(EnsD2N* X,EnsD2N* IV);
void get_right_IV(EnsD2N* X,EnsD2N* IV);
void right_multiply(EnsD2N* X,EnsD2N* Y,int p,int sign);
int right_homom_by_translation(EnsD2N* X,EnsD2N* Y,int* a,int* b);
void niceprint_EnsD2N(EnsD2N* ensemble, char* str_rep);
int is_equal(EnsD2N* X,EnsD2N* Y);

int is_left_translated(EnsD2N* X,EnsD2N* Y);
int is_right_translated(EnsD2N* X,EnsD2N* Y);
int is_pureZn(EnsD2N* X);
int mod(int a, int b);
void next_kbit_seq(int* seq,int* nextseq,int N);
void binaryprint(unsigned int a,int N);

/*
* main takes 3 arguments
*
* homometry_type: 'left' or 'right'
* N: order of the D_2N dihedral group
* P: cardinality of the subsets to be considered
* output_file: the output for the enumeration
*/

int main(int argc, char *argv[]) {
  int N,P,homomtype;
  FILE* output_file=NULL;

  if (argc<3) {
    printf("Not enough arguments - Need at least 3\n");
    exit(1);
  }
  N = atoi(argv[1]);
  P = atoi(argv[2]);
  output_file=fopen(argv[3],"w");
  if (output_file==NULL) {
    printf("Error creating the output file...\n");
    exit(1);
  }
  printf("%d %d \n",N,P);

  enumerate_homometric(N, P, output_file);

  fclose(output_file);
  exit(0);
}

void enumerate_homometric(int N, int P, FILE* output_file) {
    int total_ensembles=0;
    int i,j,flag,c;
    EnsD2N *ensembles,*ensembles_left_IV;
    EnsD2N X,X_left_IV;
    int sequence[2*N],next_sequence[2*N];
    char str_rep[200]="";
    int p,sign,endflag=0;

    printf("Building the collection of D_%d subsets of cardinality %d...\n",2*N,P);
    for(i=0;i<2*N;i++) {
      if(i<P)
        sequence[i]=1;
      else
        sequence[i]=0;
    }

    while(!endflag) {
      initialize_EnsD2N(&X, sequence, N);
      if (!is_pureZn(&X)) {
        get_left_IV(&X,&X_left_IV);
        if (total_ensembles==0) {
          // If this is the first subset found, we initialize the ensemble table
          total_ensembles++;
          ensembles = (EnsD2N*)malloc(sizeof(EnsD2N));
          initialize_EnsD2N(ensembles, sequence, N);

          ensembles_left_IV = (EnsD2N*)malloc(sizeof(EnsD2N));
          get_left_IV(ensembles,ensembles_left_IV);
        } else {
          // Otherwise, we check if the IV already exists in our enumeration

          flag=0;
          for (i=0;i<total_ensembles;i++) {
            if(is_equal(&X_left_IV,&ensembles_left_IV[i])) {
              // Most commonly, the interval vector already exists because it is
              // a right translate of a previously enumerated ensemble, in which
              // case we flag it so we now we don't need it.
              if (is_right_translated(&X,&ensembles[i])) {
                flag=1;
              }
            }
          }

          if (flag==0) {
            // If we have not flagged it before, it means that we have an ensemble
            // with a common IV, but which is not trivially related to the previous
            // ones so we add it to the enumeration.
            total_ensembles++;
            ensembles = (EnsD2N*)realloc(ensembles,total_ensembles*sizeof(EnsD2N));
            ensembles_left_IV = (EnsD2N*)realloc(ensembles_left_IV,total_ensembles*sizeof(EnsD2N));

            initialize_EnsD2N(&ensembles[total_ensembles-1], sequence, N);
            get_left_IV(&ensembles[total_ensembles-1],&ensembles_left_IV[total_ensembles-1]);
          }
        }
      }
      endflag=1;
      for(i=2*N-1;i>2*N-P-1;i--) {
        endflag &= sequence[i];
      }
      next_kbit_seq(sequence,next_sequence,2*N);
      for(i=0;i<2*N;i++)
        sequence[i]=next_sequence[i];
    }

    // Now, we need to test all ensembles to output pairs with similar left IV

    c=0;
    for (i=0;i<total_ensembles;i++) {
      for (j=i+1;j<total_ensembles;j++) {
        if(is_equal(&ensembles_left_IV[i],&ensembles_left_IV[j]) && !is_equal(&ensembles[i],&ensembles[j])) {
          fprintf(output_file,"===== %d ======\n",c);
          fprintf(output_file,(&ensembles[i])->hash);
          fprintf(output_file," - ");
          fprintf(output_file,(&ensembles[j])->hash);
          fprintf(output_file,"\n");
          niceprint_EnsD2N(&ensembles[i], str_rep);
          fprintf(output_file,str_rep);
          fprintf(output_file,"\n");
          niceprint_EnsD2N(&ensembles[j], str_rep);
          fprintf(output_file,str_rep);

          // By the above enumeration, we know that the two homometric ensembles A and B
          // are not right translates of each other. However, it may be that
          // a right translation of B gives an ensemble B' such that A and B' also
          // have the same right Interval Vector, which would mean that A and B' are both
          // left and right homometric.
          if(right_homom_by_translation(&ensembles[i],&ensembles[j],&p,&sign)) {
            if(sign==1)
              fprintf(output_file,"\n<<<< Right homometric by right multiplication of B by (%d,+) >>>>\n",p);
            else
              fprintf(output_file,"\n<<<< Right homometric by right multiplication of B by (%d,-) >>>>\n",p);
          } else {
            fprintf(output_file,"\n<<<< No right homometry obtained by right multiplication of B >>>>\n");
          }
          c++;
        }
      }
    }
}

void get_left_IV(EnsD2N* X,EnsD2N* IV) {
  /*
  * Function:  get_left_IV
  * --------------------
  * Calculates the left Interval Vector (IV) of a subset of D_2n
  *
  * X: the D_2n subset, whose left IV is calculated
  *
  *  returns: the left IV as a subset 'IV' of D_2n.
  */
  int i,j;

  IV->N = X->N;
  IV->A0 = calloc(X->N,sizeof(int));
  IV->A1 = calloc(X->N,sizeof(int));

  for(i=0;i<X->N;i++) {
    IV->A0[i] = 0;
    IV->A1[i] = 0;
  }

  for(i=0;i<X->N;i++) {
    for(j=0;j<X->N;j++) {
      IV->A0[mod((i-j),X->N)] += X->A0[i]*X->A0[j];
      IV->A0[mod((i-j),X->N)] += X->A1[i]*X->A1[j];

      IV->A1[mod((i+j),X->N)] += X->A1[i]*X->A0[j];
      IV->A1[mod((i+j),X->N)] += X->A0[i]*X->A1[j];
    }
  }
}


void get_right_IV(EnsD2N* X,EnsD2N* IV) {
  /*
  * Function:  get_right_IV
  * --------------------
  * Calculates the right Interval Vector (IV) of a subset of D_2n
  *
  * X: the D_2n subset, whose right IV is calculated
  *
  *  returns: the right IV as a subset 'IV' of D_2n.
  */
  int i,j;

  IV->N = X->N;
  IV->A0 = calloc(X->N,sizeof(int));
  IV->A1 = calloc(X->N,sizeof(int));

  for(i=0;i<X->N;i++) {
    IV->A0[i] = 0;
    IV->A1[i] = 0;
  }

  for(i=0;i<X->N;i++) {
    for(j=0;j<X->N;j++) {
      IV->A0[mod((j-i),X->N)] += X->A0[i]*X->A0[j];
      IV->A0[mod((i-j),X->N)] += X->A1[i]*X->A1[j];

      IV->A1[mod((j-i),X->N)] += X->A0[i]*X->A1[j];
      IV->A1[mod((i-j),X->N)] += X->A1[i]*X->A0[j];
    }
  }
}


int right_homom_by_translation(EnsD2N* X,EnsD2N* Y,int* a,int* b) {
  /*
  * Function:  right_homom_by_translation
  * --------------------
  * Checks if two subsets X and Y of D_2n are right homometric if Y is right
  * translated.
  *
  * X:,Y the D_2n subsets
  *
  *  returns: 0 if no right translation gives right homometric subsets,
  *          1 otherwise, in which (a,b) is the element (t^a)s^b of D_2n,
  *          by which Y must be right-translated.
  */
  EnsD2N Y_mult,X_right_IV,Y_right_IV;
  int p,sign;

  get_right_IV(X,&X_right_IV);
  for (p=0;p<X->N;p++) {
    for (sign=0;sign<2;sign++) {
        right_multiply(Y,&Y_mult,p,sign);
        get_right_IV(&Y_mult,&Y_right_IV);
        if (is_equal(&X_right_IV,&Y_right_IV)) {
          *a=p;
          *b=sign;
          return 1;
        }
    }
  }
  return 0;
}

void right_multiply(EnsD2N* X,EnsD2N* Y,int p,int sign) {
  /*
  * Function:  right_multiply
  * --------------------
  * Right translates a subset X of D_2n by the element (t^p)s^sign of D_2n
  *
  */
  int i;

  Y->N = X->N;
  Y->A0 = calloc(X->N,sizeof(int));
  Y->A1 = calloc(X->N,sizeof(int));

  for(i=0;i<X->N;i++) {
    if (sign==1) {
      Y->A0[mod(i+p,X->N)] = X->A0[i];
      Y->A1[mod(i-p,X->N)] = X->A1[i];
    } else {
      Y->A1[mod(i+p,X->N)] = X->A0[i];
      Y->A0[mod(i-p,X->N)] = X->A1[i];
    }
  }
}

int is_equal(EnsD2N* X,EnsD2N* Y){
  /*
  * Function:  is_equal
  * --------------------
  * Checks if two subsets X and Y of D_2n are equal
  * Returns 1 if this is the case, 0 otherwise.
  */
  int i;

  for(i=0;i<X->N;i++) {
    if (X->A0[i]!=Y->A0[i])
      return 0;
    if (X->A1[i]!=Y->A1[i])
      return 0;
  }
  return 1;
}

///////////////////////////////////////////////////////

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

void initialize_EnsD2N(EnsD2N* ensemble, int* sequence, int N) {
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
  int numbytes=(2*N/4)+2;

  ensemble->N = N;
  ensemble->A0 = calloc(N,sizeof(int));
  ensemble->A1 = calloc(N,sizeof(int));
  ensemble->hash = calloc(numbytes,sizeof(char));

  for(i=0;i<numbytes;i++) {
    ensemble->hash[i]=65;
  }
  ensemble->hash[numbytes-1]=0;
  for(i=0;i<2*N;i++) {
    ensemble->hash[i/4] = ensemble->hash[i/4] | (sequence[i]<<((i%4)+1));
  }

  for(i=0;i<N;i++) {
    ensemble->A0[i] = sequence[i];
    ensemble->A1[i] = sequence[i+N];
  }
}

void niceprint_EnsD2N(EnsD2N* ensemble, char* str_rep) {
  /*
  * Function:  niceprint_EnsD2N
  * --------------------
  * Print an interpretable version of a D_2n subset
  *
  * ensemble: the D_2n subset to be printed
  * str_rep: the destination string. Prints an element (g,h) of the subset
  *          as g+ if h=1, g- otherwise.
  *
  *
  *  returns: None.
  *
  */
  int i;
  char* temp_str;

  strcpy(str_rep,"");
  // 200 is the max size of the string
  // Not the better way to code it, though
  snprintf(str_rep, 200, "%s%s", str_rep, "{");
  for(i=0;i<ensemble->N;i++) {
    if (ensemble->A0[i])
      snprintf(str_rep, 200, "%s%d+,", str_rep, i);
  }
  for(i=0;i<ensemble->N;i++) {
    if (ensemble->A1[i])
      snprintf(str_rep, 200, "%s%d-,", str_rep, i);
  }
  snprintf(str_rep, 200, "%s%s", str_rep, "}");
}

void next_kbit_seq(int* seq,int* nextseq,int N) {
  /*
  * Function:  next_kbit_seq
  * --------------------
  * From a sequence of k bits in N, returns the next sequence of k bits in N
  *
  * seq: array of ints of size N
  * nextseq: array of ints of size N, in which the next sequence is written
  * N: int, size of the sequence
  *
  */
  int i,smallest,next_smallest,c,nc;
  int ripples[N],ones[N],final[N];

  for(i=0;i<N;i++) {
    ripples[i]=seq[i];
    ones[i]=0;
    final[i]=0;
  }

  smallest=0;
  while(seq[smallest]==0)
    smallest+=1;

  ripples[smallest]=0;
  c=1;
  for(i=smallest+1;i<N;i++) {
    nc = c & ripples[i];
    ripples[i] = (ripples[i]+c)%2;
    c=nc;
  }

  next_smallest=0;
  while(ripples[next_smallest]==0)
    next_smallest+=1;

  for(i=0;i<(next_smallest-smallest-1);i++)
    ones[i]=1;
  for(i=0;i<N;i++)
    nextseq[i] = ripples[i] | ones[i];
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
