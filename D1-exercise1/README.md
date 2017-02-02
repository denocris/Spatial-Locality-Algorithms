#Bit manipulation routines#

In what follows to help the students some exercises are presented with basic
pseudo-codes. Herafter a pseudocode are several lines of istructions with the aim
of illustrating the general procedure. The solution of the problems requires the
conversion of these pseudocodes into a computer program , it is left to the student
the choice of the programming language. The pseudocodes presented here will be
written in Fortran-style language.
There is not a conventional pseudocode language for bit manipulation. Here
for simplicity we make use of the Fortran names for the functions

The aim of this simple problem is to acquire some familiarity with the bit
manipulations routines which will be used later in the forthcoming exercises.

Write a program in which a variable K is declared as a 64-bit integer long.
On this variable perform the following operations (assuming that bit positions start from the left and the range is 0 − 63):

- set bit in position 3 to one (IBSET)
- check the non-zero bits (BTEST)
- now set to one the first three bits (IBSET)
- move them to the right by two positions (ISHFT)
- is bit 4 = 1 ? (BTEST)
- now from position 2 extract three bits (IBITS)
- set bit at position 3 to be zero (IBCLR)
- now perform a logical mask with another integer (K2 = 16) (IAND)
- move to the left by three positions (ISHFT)

Print the value of K after each operation.