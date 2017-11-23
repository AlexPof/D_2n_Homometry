# D_2n Homometry

C algorithms for brute-force enumeration of homometric sets in the dihedral groups D_2n

==========

Start the C program:

    >>> ./D2n_homomenumerate homometry_type n p output_file

where

  * *homometry_type* is either *left* or *right*
  * *n* is an integer defining the order of the D_2n dihedral group
  * *p* is an integer defining the cardinality of the subsets of D_2n to be
        examined.
  * *output_file* is the name of the output file to be written

For example:

    >>> ./D2n_homomenumerate left 12 5 output.txt
    Building the collection of D_24 ensembles of cardinality 5...
    Found 1705 D_24 ensembles of cardinality 5...
    Determining left-homometric D_24 ensembles of cardinality 5...

Use the python script to count the unique homometric n-uples:

    >>> python D2n_homomcounts.py

With the above example:

    >>> python D2n_homomcounts.py
    Found 10 unique homometric sets:
    [[20753, 82193], [41233, 164113], [86081, 267269, 267521], [143425, 268417], [172097, 533509, 533761], [200769, 268305], [266499, 266625], [532739, 532865], [1065219, 1065345], [2130179, 2130305]]
    Counts by n-uples:
      Number of 2-uples: 8
      Number of 3-uples: 2
