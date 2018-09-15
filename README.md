# D_2n Homometry

C algorithm for brute-force enumeration of homometric sets in the dihedral groups D_2n

See the paper *Genuys G., Popoff A. (2017) Homometry in the Dihedral Groups: Lifting Sets from  ℤn to  Dn. In: Agustín-Aquino O., Lluis-Puebla E., Montiel M. (eds) Mathematics and Computation in Music. MCM 2017. Lecture Notes in Computer Science, vol 10527. Springer, Cham* for more theoretical information.

==========

Compile the C program with gcc and start:

    >>> ./D2n_homomenumerate n p output_file

where

  * *n* is an integer defining the order of the D_2n dihedral group
  * *p* is an integer defining the cardinality of the subsets of D_2n to be
        examined.
  * *output_file* is the name of the output file to be written

For example:

    >>> ./D2n_homomenumerate 12 5 output.txt

Use the python script to count the unique homometric n-uples:

    >>> python D2n_homomcounts.py output.txt

With the above example:

    >>> python D2n_homomcounts.py
    # of left homometric subsets
    2-uples: 8
    3-uples: 2
    # of simultaneous left and right homometric subsets
    2-uples: 8
    3-uples: 2

The zip file *partial_results.zip* contains a partial enumeration of subsets of cardinality *p* in the dihedral groups D_2n, for *p=4,5,6,7* and n from 8 to 18.
