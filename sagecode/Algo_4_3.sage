"""Paper referenced: https://services.math.duke.edu/~dasgupta/papers/comp.pdf
This section uses the integrals calculated before to run one final algorithm to compute the elliptic units"""
import sympy,math, gmpy2

load("Algo_4_2a.sage")
load("Algo_4_2b.sage")
load("Basis.sage")

S = Qp(p, prec = 2*M + 1, type = 'capped-rel', print_mode = 'series')	#Defines the padic ring of rationals, with specified precision.
R.<y> = PowerSeriesRing(S,sparse=True)	#Allows us to write polynomials in y with coefficients in S.

"""Helper method used to calculate ord_p. See algorithm 4.3, step 1"""
def order_pHelper(T,p,N,a,c):
    """D is just D_{1,1}. alpha denotes the unit"""
    if N != 4:
        """N is prime and unit is [N]-N[1]. 
        When the sum of nd is nonzero (when N is prime), there must be an error correction term in the sum
        The (N-1)/4 correction comes from the error term when sum nd != 0."""
        order = N*genDedekind(1,1,a,N*c)-genDedekind(1,1,a,c) + (N-1)/4
    else:
        """If N = 4 and unit is 2[1]-3[2]+[4]:"""
        order = (-2*genDedekind(1,1,a,4*c)+3*genDedekind(1,1,a,2*c)-genDedekind(1,1,a,c))*2
    return order 

"""See algorithm 4.3, step 1. Calculates ord_p"""
def order_p(T,p,N,a,c,gammas):
    """D is just D_{1,1}, the generalized Dedekind sum. alpha denotes unit is [N]-N[1]. 
    We reduce the calculation by writing [inf]-[a/Nc] in terms of the basis elements, 
    which are calculated by order_pHelper."""
    
    """We call the helper method on each element in the basis and take the sum"""
    totalOrder = 0
    for gamma in gammas[:-1]:
        base = gamma.x
        totalOrder += gamma.coeff*order_pHelper(T,p,N,base,1) 
        """When sum of nd is not 0, we must account for sign change
        Eg. a/c and -a/-c are different when sum nd is not 0. See Eqn 43 of DD2006"""
        if sgn(a)>0 and N != 4 and sgn(gamma.a*base+N*gamma.b) == -1:
            totalOrder += (N-1)/2
        elif sgn(a)<0 and N != 4 and sgn(gamma.a*base+N*gamma.b) == 1:
            totalOrder += (N-1)/2
    if sgn(a) < 0 and N!=4:
        totalOrder+= (N-1)/2
    
    """If N = 2 or 3, we need to multiply by 4 and 3 respectively to get integral values."""
    if N == 2:
        totalOrder *= 4
    elif N == 3:
        totalOrder *= 3
    print("Order of p: %d" %totalOrder)
    return totalOrder

"""Calculates the right hand side of equation 44, which is used to calculate log_beta
See algorithm 4.3, step 2"""
def rhsOfEqn44(a,c,u,v,p,N,s=1):
    """Prop 3.2 but s=1 default. 
    Returns u([infty] - [a/Nc]) for U_(u,v,1)"""
    total = 0
    
    """Reduce a/c"""
    if c < 0:
        a, c = -a, -c
    div = gcd(a,c)
    a, c = a/div, c/div
    
    if N != 4:
        """If N is prime and unit is [N]-N[1]"""
        for l in range(N*c):
            total -= bernPolyTilde(1, a*(l+v/p^s)/(N*c)-u/p^s)*(bernPolyTilde(1,1/c*(l+v/p^s))-N*bernPolyTilde(1,(l+v/p^s)/(N*c)))
    else:
        """If N = 4 and unit is 2[1]-3[2]+[4]:"""
        for l in range(4*c):
            total -= 2*bernPolyTilde(1, a*(l+v/p^s)/(4*c)-u/p^s)*(bernPolyTilde(1,1/c*(l+v/p^s))
                                                               -3*bernPolyTilde(1,(l+v/p^s)/(2*c))
                                                               +2*bernPolyTilde(1,(l+v/p^s)/(4*c)))
    """If N = 2 or 3, we need to multiply by 4 and 3 respectively to get integral values."""
    if N == 2:
        total *= 4
    elif N == 3:
        total *= 3
    return total

"""Used in the log_beta helper method."""
def log_betaTable(table,x):
    """Finds discrete logarithm of x with base beta"""
    for index in range(len(table)):
        if (table[index]-x).valuation()>0:
            return index
    raise "error"

"""Helper method used to calculate log_beta. See Algorithm 4.3, step 2"""
def log_betaHelper(T,a,c,p,beta,tau,N,table,Gamma = None):
    """If Gamma is given, then this is a sum augmented by Gamma."""
    try: 
        A,B,C,D = Gamma.a, Gamma.b, Gamma.c, Gamma.d 
        gammaInv_tau = T(D*tau-B)/T(-C*tau+A)
        total = 0
        for u in range(p):
            for v in range(p):
                if u == 0 and v == 0:
                    pass
                else:
                    total += log_betaTable(table,T(A*u+B*v-(C*u+D*v)*tau))*rhsOfEqn44(a,c,u,v,p,N)
    except:
        """Returns usual log_beta for a/Nc"""
        total = 0
        for u in range(p):
            for v in range(p):
                if u == 0 and v == 0:
                    pass
                else:
                    total += log_betaTable(table,T(u-v*tau))*rhsOfEqn44(a,c,u,v,p,N)
    """Calculate remainder (mod p^2-1) and return"""
    return total

def log_beta(T,a,c,p,beta,tau,N,gammas):
    """Construct discrete logarithm table of beta"""
    table = [T(1)]
    for i in range(1,p^2-1): 
        table.append(table[i-1]*beta)
    
    """Call the helper method on each element of the basis, take the sum """
    total = 0
    for gamma in gammas[:-1]:
        ga, gb, gc, gd = gamma.a, gamma.b, gamma.c, gamma.d 
        base = gamma.x
        total += gamma.coeff*(log_betaHelper(T,base,1,p,beta,tau,N,table,Gamma = gamma))

    print("Order of beta:", total%(p^2-1))
    return total%(p^2-1)

def sumOfMeasures(i,p,a,c,N,s=1):
    """Finds the sum of measures u_bar_m(i+p^sZ) after the pushforward. By default s=1."""
    total = 0
    Modp = IntegerModRing(p^s)
    for y in range(p^s):
        if y%p != 0:
            for x in range(p^s):
                if (x*(Modp(y)^-1)-i)%(p^s) == 0:
                    total += rhsOfEqn44(a,c,x,y,p,N,s)
    return total

"""Helper method used to calculate log_p. See Algorithm 4.3, step 3"""
def log_pHelper(T,taup,i,integral_i,integral_ii,integral_iii,integral_iv):
    """i is the basis m_i, which is [inf] - [i/N] in our case. 
    Returns integral in RHS of 45. 
    Note, equation (45) above should have x-y*tau replace x-y*gamma^-1tau"""
    
    """First and third integrals as calculated in Algo_4_2"""
    first = T(integral_i[i])
    third = T(integral_iii[i])
    
    """Second integral, using equation 36"""
    second = T(0) 
    for j in range(p):
        """first part of eqn 36, can use integral_ii[(i,j,0)] or sumOfMeasures(j,p,i,1,N)"""
        second += log(taup-j)*T(integral_ii[(i,j,0)])
        """second part of eqn 36"""
        for n in range(1,M):
            output = T(integral_ii[(i,j,n)])
            second -= T(output)/((taup-j)^n*n)
    
    """Fourth integral"""
    forth = T(0)
    for n in range(1,M):
        output = T(integral_iv[(i,n)])
        forth -= (taup^n*output)/n
        
    return first + second + third + forth

def log_p(T,a,c,p,N,tau,integral_i,integral_ii,integral_iii,integral_iv,gammas):
    """a_gamma,i can be assume to be +-1 throughout, since we will iterate over all gammas"""
    total = 0
    
    """Call the helper method on each element of the basis, take the sum """
    for gamma in gammas[:-1]:
        ga, gb, gc, gd = gamma.a, gamma.b, gamma.c, gamma.d 
        base = gamma.x
        gammaInv_tau = (gd*tau-gb)/(-gc*tau+ga)
        total += gamma.coeff*log_pHelper(T,gammaInv_tau,base,integral_i,integral_ii,integral_iii,integral_iv)
    print("log_p value: " + str(total + O(p^5)))
    return total

"""Final step of algorithm 4.3, calculates u(alpha, tau)"""
def u_alpha_tau(T,a,c,p,N,tau,beta,integral_i,integral_ii,integral_iii,integral_iv):
    """Write [inf]-[a/Nc] in terms of basis elements"""
    gammas = getGammas(a,c*N)
    """Calculate the three components"""
    exp = log_p(T,a,c,p,N,tau,integral_i,integral_ii,integral_iii,integral_iv,gammas)
    exp = (exp).exp()
    ord_p = order_p(T,p,N,a,c,gammas)
    log_b = log_beta(T,a,c,p,beta,tau,N,gammas)
    return T(p^ord_p) * beta^log_b * exp


