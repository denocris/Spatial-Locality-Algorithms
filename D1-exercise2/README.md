#Morton order#
Let us consider a square of side length ``L``. Write a program which computes 
the coordinates ``xpos[i]`` and ``ypos[i]`` of a set of points ``N = M × M``. 
The points must be arranged into the square with a uniform spacing. However,
the program must compute the final list of points hierarchically, that is the
program must contain a recursive function which has as input the root lattice 
(``N = 4 = 2 × 2``) of points and it returns the 4N points which fill the
four sub-squares. The order in which the sub-quadrants are filled must be:
left-bottom, right-bottom, left-top, right top.
The function must be recursive so that for a given recursion depth levscan
it proceeds until there are 4 levscan points. 

Finally write the final points together 
with the corresponding Morton (or Z−) key values and make a plot
of the points with a line joining them. Verify that for each point ``i > 1``
it is satisfied ``key(i) > key(i − 1)``. 
Compute the Z−keys up to the order levmorton.
Assume as input ``L = 16``, ``levscan = 4``, ``levmorton = 8``.
Here are given the corresponding [pseudocodes](exhpc_2.pdf).
