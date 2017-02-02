#Hilbert order#

This problem is similar to the one given in PS2. Here, however, one must
compute the Hilbert key for a set of points. The construction of the grid
points proceeds along the lines described in PS2, but with some differences.
Let us consider a square of side length ``L``. Write a program which computes 
the coordinates ``xpos[i]`` and ``ypos[i]`` of a set of points ``N = M × M``. 
The
points must be arranged into the square with a uniform spacing. However,
the program must compute the final list of points hierarchically, that is the
program must contain a recursive function which has as input the root lattice
(``N = 4 = 2 × 2``) of points and it returns the ``4N`` points which fill the four
sub-squares. The order in which the sub-quadrants are filled now must be :
left-bottom, left-top , right-top, right-bottom.
An important difference with respect PS2 is that here the coordinates of
the points filling the sub-squares must not be a simple replica of the parent
square, but their coordinates must be transformed as follows. 
* Left-bottom: rotate −90 deg and reverse order. 
* Left-top and right-top: identical. 
* Right-bottom rotate +90 deg and reverse order.
The function must be recursive so that for a given recursion depth levscan
it proceeds until there are 4 levscan points. Finally write the final points to-
gether with the corresponding Hilbert (or H−) key values and make a plot
of the points with a line joining them. Verify that for each point ``i > 1``
it is satisfied ``key(i) > key(i − 1)``. Compute the H−keys up to the order
``levhilbert = levscan``.
Here are given the corresponding [pseudocodes](exhpc_3.pdf).