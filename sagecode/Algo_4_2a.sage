"""
This section creates some necessary helper methods and calculates integrals i and iii for algorithm 4.2
Integrals ii and iv will be calculated in Algo_4_2b
Paper referenced: https://services.math.duke.edu/~dasgupta/papers/comp.pdf
"""
import sympy,math, gmpy2, time
load("Miscellaneous.sage")

S = Qp(p, prec = M*2+1, type = 'capped-rel', print_mode = 'series')	#Defines the padic ring of rationals, with specified precision.
R.<y> = PowerSeriesRing(S, sparse = True)	#Allows us to write polynomials in y with coefficients in S.  


"""These methods are used in several different places in the code"""

def bernPolyList(p,M):
    """Get the list of Bernoulli Polynomials. Load from pre-existing file for any polynomial already calculated."""
    bernPolys = []
    
    """Tried storing the bernoulli polynomials to improve run time. 
    This takes up too much space, so it has been commented out. """
#     name = "bernPoly.txt"
#     file = open(name, 'w+')
#     try:
#         content = file.readlines()
#     except:
#         content = []
#     for bernPolyString in content:
#         bernPolys.append(sage_eval(bernPolyString, locals={'z':z}))
        
#     for i in range(len(content),p*(M+math.ceil(log(M)))):
#         bernPolys.append(bernoulli_polynomial(z,i))
#         file.write(str(bernPolys[-1]) + "\n")
        
#     file.close()

    """Add each Bernoulli polynomial to the list using SageMath's bernoulli_polynomial function"""
    for i in range(p*(M+math.ceil(log(M)))):
        bernPolys.append(bernoulli_polynomial(y,i))
    return bernPolys
    
def bernPolyTilde(s, x):  
    """Returns: Tilde bernoulli Polynomial (Pg.6)"""
    if s == 1 and round(x) == x:
        return 0
    else:
        return bernoulli_polynomial(x-math.floor(x),s)
    
def bernPolyTildeZp(s, x, bernPolys):  
    """Returns: Tilde bernoulli Polynomial in Zp (Pg.6)"""
    if s == 1 and round(x) == x:
        return 0
    else:
        return bernPolys[s](y=x-math.floor(x))

def genDedekind(s,t,a,c):
    """Returns the generalized Dedekind sum (Pg.7)"""
    output = 0
    if c < 0:
        a, c = -a, -c
    div = gcd(a,c)
    a, c = a/div, c/div

    for i in range(1,c+1):
        output += bernPolyTilde(s,i/c)*bernPolyTilde(t,i*a/c)
    output = output*(c^(s-1))/(s*t)
    return output

def genDedekindFast(s,t,a,c,bernPolys):
    """Returns the p-adic of generalized Dedekind sum (Pg.7)"""
    output = 0
    if c < 0:
        a, c = -a, -c
    div = gcd(a,c)
    a, c = a/div, c/div

    for i in range(1,c+1):
        output += S(bernPolyTildeZp(s,i/c,bernPolys)*bernPolyTildeZp(t,i*a/c,bernPolys))
    output = output*(c^(s-1))/(s*t)
    return output

def eqn14(p,a,N,c,n,k,bernPolys):
    """Equation 14 of the paper. prime p is fixed from start. Used to calculate integral 1 in algorithm 4.2"""
    result = 0
    for l in range(n+1): 
        """
        sumofDede should iterate over all d dividing N. 
        """
        if N != 4:
            """We assume N is prime here and the given unit is [N]-N[1].
            Hence only have to consider d=1 and N. n_N = 1 and n_1 = -N"""
            sumofDede = S((genDedekindFast(k-l-1,l+1,a,c,bernPolys)-p^(k-l-2)*genDedekindFast(k-l-1,l+1,p*a,c,bernPolys))/(N^l)
                         -N*(genDedekindFast(k-l-1,l+1,a,N*c,bernPolys)-p^(k-l-2)*genDedekindFast(k-l-1,l+1,p*a,N*c,bernPolys)))
            """If N = 2 or 3, we need to multiply by 4 and 3 respectively to get integral values."""
            if N == 2:
                sumofDede *= 4
            elif N == 3:
                sumofDede *= 3
        else:
            """If N = 4 and unit is 2[1]-3[2]+[4], or n_1 = 2, n_2 = -3, n_4 = 1"""
            sumofDede = 2*S((genDedekindFast(k-l-1,l+1,a,c,bernPolys)-p^(k-l-2)*genDedekindFast(k-l-1,l+1,p*a,c,bernPolys))/(4^l)
                         -3*(genDedekindFast(k-l-1,l+1,a,2*c,bernPolys)-p^(k-l-2)*genDedekindFast(k-l-1,l+1,p*a,2*c,bernPolys))/(2^l)
                         +2*(genDedekindFast(k-l-1,l+1,a,4*c,bernPolys)-p^(k-l-2)*genDedekindFast(k-l-1,l+1,p*a,4*c,bernPolys)))
        result -= S(gmpy2.comb(n,l) * (a/(N*c))^(n-l) * (-1)^l) * sumofDede
    return result

def logPowerSeries(i,p,M):
    """Returns the power series of log_p(M) around residue i, with precision M+int(log(M)).
    w is the teichmuller representative of i. If y=i(modp), log_p(y)=log_p(y/w)."""
    w = S.teichmuller(i)
    f = 0
    for j in range(1,M+math.ceil(log(M))):
        f -= (1/j)*(-1)^j*(y/w-1)^j     
    return f

def gSubi35(i,p,M):  
    """Returns the function g_i(y) as in equation 35"""
    g = 1
    prec = M+math.ceil(log(M))
    for j in range(1,p):
        if j != i:
            g *= (y-j)^prec
    return g 

def hSubi35(i,p,M): 
    """Returns the power series of h_i(y)=log_p(y)/g_i(y) in eqn 35."""
    prec = M+math.ceil(log(M))
    g_i = gSubi35(i,p,M)(y+i)	#This step ensures that the power series expansion is expanded around i. 
    powerg_i = 1/(g_i+O(y^prec))
    powerg_i = powerg_i.polynomial()
    return powerg_i(y-i) * logPowerSeries(i,p,M)

def fEqnlogApprox35(p,M):
    """Returns the f(y) that approximates log_p(y) in eqn 35."""
    f = 0
    for i in range(1,p):
        f += gSubi35(i,p,M)*hSubi35(i,p,M).polynomial()
    return f

def integral_i_helper(f,p,a,N,c,M,bernPolys):
    """Caclulates integral i) in Algo 4.2 to precision M for a given element of the basis"""
    total = 0
    exponents = f.exponents()
    coeffs = f.coefficients()
    if len(coeffs) != len(exponents) : 
        raise Error("Coefficients and Exponents do not match.") 
    totalLen = len(exponents)
    for i in range(totalLen):
        total += coeffs[i]*S(eqn14(p,a,N,c,0,exponents[i]+2,bernPolys)) 
    return total 

def integral_i_dict(f,p,N,M):
    """Finds the integral i) for each element in the basis: [inf]-[1/N],...,[inf]-[(N-1)/N]. 
    Helper method calculates each individual integral. This method returns these as a dictionary"""
    integrals = {}
    bernPolys = bernPolyList(p,M)
    for i in range(1,N):
        integrals[i] = integral_i_helper(f,p,i,N,1,M,bernPolys)
        print("Calculated integral 1:%d/%d" %(i,N-1))
    return integrals

def eqn41(p,a,N,c,n,bernPolys):
    """Returns the integral at Equation 41."""
    result = 0
    for l in range(n+1):
        if N != 4:
            """As we assume N is prime, we only consider d=1 and N. 
            Here n_N = 1 and n_1 = -N"""
            sumofDede = S((p^(n-l)*genDedekindFast(n-l+1,l+1,p*a,c,bernPolys)-p^n*genDedekindFast(n-l+1,l+1,a,c,bernPolys))/(N^l)
                        -N*(p^(n-l)*genDedekindFast(n-l+1,l+1,p*a,N*c,bernPolys)-p^n*genDedekindFast(n-l+1,l+1,a,N*c,bernPolys)))
            """If N = 2 or 3, we need to multiply by 4 and 3 respectively to get integral values."""
            if N == 2:
                sumofDede *= 4
            elif N == 3:
                sumofDede *= 3
        else:
            """If N = 4 and unit is 2[1]-3[2]+[4]:"""
            sumofDede = 2*S((p^(n-l)*genDedekindFast(n-l+1,l+1,p*a,c,bernPolys)-p^n*genDedekindFast(n-l+1,l+1,a,c,bernPolys))/(4^l)
                        -3*(p^(n-l)*genDedekindFast(n-l+1,l+1,p*a,2*c,bernPolys)-p^n*genDedekindFast(n-l+1,l+1,a,2*c,bernPolys))/(2^l)
                         +2*(p^(n-l)*genDedekindFast(n-l+1,l+1,p*a,4*c,bernPolys)-p^n*genDedekindFast(n-l+1,l+1,a,4*c,bernPolys)))
        result -= S(gmpy2.comb(n,l) * (a/(N*c))^(n-l) * (-1)^l) * sumofDede
    return result     
    
def integral_iii_helper(f,p,a,N,c,M,bernPolys):
    """Returns integral iii), with inputs a/Nc, and f which is the polynomial approximation of log_p(x)"""
    total = 0
    exponents = f.exponents()
    coeffs = f.coefficients()
    totalLen = len(exponents)
    for i in range(totalLen):
        total += coeffs[i]*S(eqn41(p,a,N,c,exponents[i],bernPolys)) 
    return total 

def integral_iii_dict(f,p,N,M):
    """Finds the integral iii) for each element in the basis: [inf]-[1/N],...,[inf]-[(N-1)/N]."""
    integrals = {}
    bernPolys = bernPolyList(p,M)
    for i in range(1,N):
        integrals[i] = integral_iii_helper(f,p,i,N,1,M,bernPolys)
        print("Calculated integral 3:%d/%d" %(i,N-1))
    return integrals 


