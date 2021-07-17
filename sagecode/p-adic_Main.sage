"""To initialize, specify:
p: prime for the p-adic calculations.
N: Defines the unit used for computation. In this code we implemented the units
    N=4        - [4]-3[2]+2[1]
    N is prime - [N]-N[1]
M: p-adic precision. Ideally this value is at least 100.
"""
p = 3
N = 2
M = 100

load("Miscellaneous.sage")
#Defines the padic ring of rationals, with specified precision.
S = Qp(p, prec = 2*M+1, type = 'capped-rel', print_mode = 'series')
#Allows us to write polynomials in y with coefficients in S.
R.<y> = PowerSeriesRing(S, sparse=True)

"""Checks to see if a file exists with the necessary dictionaries of integrals given p, N, M
If so, loads the integrals from those files to save time.
Otherwise, runs algorithms 4_2a and 4_2b to calculate those integrals
Note that especially depending on M, this may take several hours."""

"""WARNING: Do not terminate process when integrals are being loaded. Else the file storing the integrals will be erased."""
integral_1, integral_2, integral_3, integral_4 = getIntegrals(p,N,M)

"""Find discriminants D such that (D/p)=-1, (D/N)=1, D=1(mod 4) when N is prime and =1(mod 8) when N=4"""
Ds = getGoodDs(N,p,1000)

"""For each fundamental discriminant, calculate polynomial"""
for D in Ds:
    getPoly(D, N, integral_1, integral_2, integral_3, integral_4)

