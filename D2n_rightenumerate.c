#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


typedef struct EnsD2N {
	int* A0;
	int* A1;
  unsigned long long int  val;
	int N;
} EnsD2N;

void initialize_EnsD2N(EnsD2N* ensemble, unsigned long long int  a, int N);
void niceprint_EnsD2N(EnsD2N* ensemble);

int is_right_homometric(EnsD2N* X,EnsD2N* Y);
int is_left_translated(EnsD2N* X,EnsD2N* Y);
int is_pureZn(EnsD2N* X);
int mod(int a, int b);
unsigned long long int  power2(unsigned long long int  a);

unsigned long long int  next_samebits_number(unsigned long long int  x);
void binaryprint(unsigned int a,int N);


int main(int argc, char *argv[]) {
    EnsD2N* ensembles;
    int total_ensembles=0;
    EnsD2N X,Y;
    unsigned long long int  x;
    int i,j,flag;
    int count;

    int N,P;

   	if (argc<3) {
        printf("Not enough arguments - Need at least 2\n");
        exit(1);
    } else {
		    N = atoi(argv[1]);
        P = atoi(argv[2]);
    }

    ensembles = (EnsD2N*)malloc(sizeof(EnsD2N));

    /////////////////////////////////
    //  The enumeration proceeds in two parts:
    //  - First, we determine a collection of D_2n ensembles that are not related by a right transform of D_2n
    //    and which are not ensembles of pure Z_n

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
                            if (is_left_translated(&X,&ensembles[i]))
                                flag=1;
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

    /////////////////////////////////
    //  - Then, we determine the left-homometric sets

    printf("Found %d D_%d ensembles of cardinality %d...\n",total_ensembles,2*N,P);
    /////////////////////////////////
    //  - Then, we determine the left-homometric sets

    printf("Determining right-homometric D_%d ensembles of cardinality %d...\n",2*N,P);
    count=0;
    for(i=0;i<total_ensembles;i++) {
        for(j=i+1;j<total_ensembles;j++) {
            if (is_right_homometric(&ensembles[i],&ensembles[j])) {
                printf("%llu-%llu\n",ensembles[i].val,ensembles[j].val);
                niceprint_EnsD2N(&ensembles[i]);
                niceprint_EnsD2N(&ensembles[j]);
                printf("*******************\n");
                count++;
            }
        }
    }

    printf("%d right homometric ensembles found\n", count);
}

///////////////////////////////////////////////////////

int is_right_homometric(EnsD2N* X,EnsD2N* Y){
/*
 * Function:  is_right_homometric 
 * --------------------
 * XXX:
 *    XXXX
 *
 *  X: XXX
 *
 *  returns: XXX
 *           XXX
 *           XXX
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
    int shift,i,N,c,d;

    N=X->N;
    for(shift=0;shift<N;shift++) {
        c=0;
        d=0;
        for(i=0;i<N;i++) {
            if ( (X->A0[i] == Y->A0[ mod((i+shift),N)]) )
                c++;
            if ( (X->A1[i] == Y->A1[ mod((i+shift),N)]) )
                c++;

            if ( (X->A0[i] == Y->A1[ mod((shift-i),N)]) )
                d++;
            if ( (X->A1[i] == Y->A0[ mod((shift-i),N)]) )
                d++;
        }
        if (c==(2*N) || d==(2*N))
            return 1;
    }

    return 0;
}

int is_pureZn(EnsD2N* X) {
    int i,N,c,d;

    N=X->N;
    c=0;
    d=0;
    for (i=0;i<N;i++) {
        c+=X->A0[i];
        d+=X->A1[i];
    }
    if (c==0 || d==0)
        return 1;
    return 0;
}


///////////////////////////////////////////////////////

void initialize_EnsD2N(EnsD2N* ensemble, unsigned long long int  a, int N) {
    int i;

    ensemble->N = N;
    ensemble->val = a;
    ensemble->A0 = calloc(N,sizeof(int));
    ensemble->A1 = calloc(N,sizeof(int));

    for(i=0;i<ensemble->N;i++) {
        ensemble->A0[i] = (a>>i) & 1;
    }
    for(i=0;i<ensemble->N;i++) {
        ensemble->A1[i] = (a>>(N+i)) & 1;
    }
}

void niceprint_EnsD2N(EnsD2N* ensemble) {
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

//////////////

unsigned long long int  next_samebits_number(unsigned long long int  x) {
    unsigned long long int  smallest, ripple, new_smallest, ones;

    smallest     = (x & -x);
    ripple       = x + smallest;
    new_smallest = (ripple & -ripple);
    ones         = ((new_smallest/smallest) >> 1) - 1;
    return ripple | ones;
}


int mod(int a, int b) {
    return (a >= 0 ? a % b : b - (-a) % b);
}

unsigned long long int  power2(unsigned long long int  a) {
    unsigned long long int v=1;
    return (v << a);
}

void binaryprint(unsigned int a,int N) {
    int i;
    for(i=0;i<N;i++)
        printf("%d",(a>>i)&1);
    printf("\n");
}
