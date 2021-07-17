def loadIntegrals(fileName, N):
    """
    Load integrals that are already calculated.
    Returns: Dictionary of integrals for i), ii), iii), iv) in Algorithm 4.2
    """
    print("Loading integrals from preexisting file.")
    with open(fileName, "r") as file: 
        integral_1_string = file.readline()
        integral_2_string = file.readline()
        integral_3_string = file.readline()
        integral_4_string = file.readline()
    
    integral_1_string = integral_1_string.replace(",", ":")
    integral_1_string = integral_1_string.replace("}", ":")
    integral_1_string = integral_1_string.split(":")
    integral_1 = {}

    for i in range(1,N):
        numbers = integral_1_string[2*i-1].split("+")
        if sage_eval(numbers[0]) == 0:
            integral_1[i] = 0
        else:
            precision = numbers.pop(-1).strip()[2:-1]
            integral_1[i] = O(S(sage_eval(precision)))
            for number in numbers:
                integral_1[i] += S(sage_eval(number))
#     print(integral_1)

    integral_2_string = integral_2_string.replace(",", ":")
    integral_2_string = integral_2_string.replace("}", ":")
    integral_2_string = integral_2_string.split(":")
    integral_2 = {}

    counter = 1
    for i in range(1,N):
        for j in range(p):
            for n in range(M):
                numbers = integral_2_string[4*counter-1].split("+")
                if sage_eval(numbers[0]) == 0:
                    integral_2[(i,j,n)] = 0
                else:
                    precision = numbers.pop(-1).strip()[2:-1]
                    integral_2[(i,j,n)] = O(S(sage_eval(precision)))
                    for number in numbers:
                        integral_2[(i,j,n)] += S(sage_eval(number))
                counter += 1
#     print(integral_2)

    integral_3_string = integral_3_string.replace(",", ":")
    integral_3_string = integral_3_string.replace("}", ":")
    integral_3_string = integral_3_string.split(":")
    integral_3 = {}

    for i in range(1,N):
        numbers = integral_3_string[2*i-1].split("+")
        if sage_eval(numbers[0]) == 0:
            integral_3[i] = 0
        else:
            precision = numbers.pop(-1).strip()[2:-1]
            integral_3[i] = O(S(sage_eval(precision)))
            for number in numbers:
                integral_3[i] += S(sage_eval(number))
#     print(integral_3)

    integral_4_string = integral_4_string.replace(",", ":")
    integral_4_string = integral_4_string.replace("}", ":")
    integral_4_string = integral_4_string.split(":")
    integral_4 = {}

    counter = 1
    for i in range(1,N):
        for n in range(M):
            numbers = integral_4_string[3*counter-1].split("+")
            if sage_eval(numbers[0]) == 0:
                integral_4[(i,n)] = 0
            else:
                precision = numbers.pop(-1).strip()[2:-1]
                integral_4[(i,n)] = O(S(sage_eval(precision)))
                for number in numbers:
                    integral_4[(i,n)] += S(sage_eval(number))
            counter += 1
#     print(integral_4)
    print("All integrals loaded.")

    return integral_1, integral_2, integral_3, integral_4

def getIntegrals(p,N,M):
    """
    Obtain the integrals given 
    p: prime for the p-adic calculations.
    N: Defines the unit used for computation. In this code we implemented the units
        N=4 - [4]-3[2]+2[1]
        N=prime - [N]-N[1]
    M: p-adic precision. Ideally this value is at least M=100.
    Either load the integrals from a preexisting file, or calculates it using the algorithms. 
    Returns: Dictionaries of integrals for i), ii), iii), iv) in Algorithm 4.2
    """
    fileName = "p-%d,N-%d,M-%d.txt"%(p,N,M)
    
    try:
        """Attempt to load integrals from preexisting file."""
        integral_1, integral_2, integral_3, integral_4 = loadIntegrals(fileName, N)
        load("Basis.sage")
        load("Algo_4_2a.sage")
        load("Algo_4_2b.sage")
        load("Algo_4_3.sage")
    except:
        """Calculate the integrals from new file. Integrals will be saved under fileName."""
        print("Did not find valid file containing prexisting integrals. Calculating integrals from algorithms.")
        file = open(fileName, "w")
        """In the following we calculate all the auxilliary integrals. """
        load("Basis.sage")
        load("Algo_4_2a.sage")
        f = fEqnlogApprox35(p,M)
        integral_1 = integral_i_dict(f,p,N,M)
        print("Integral 1 completed.")
        integral_3 = integral_iii_dict(f,p,N,M)
        print("Integral 3 completed.")
        load("Algo_4_2b.sage")
        integral_2 = integral_ii_dict(p,N,M)
        print("Integral 2 completed.")
        integral_4 = integral_iv_dict(p,N,M)
        print("Integral 4 completed.")
        load("Algo_4_3.sage")
        
        file.write(str(integral_1) + "\n")
        file.write(str(integral_2) + "\n")
        file.write(str(integral_3) + "\n")
        file.write(str(integral_4) + "\n")
        file.close()
        
        print("All integrals calculated.")
    # print(integral_1)
    # print(integral_2)
    # print(integral_3)
    # print(integral_4)
    
    return integral_1, integral_2, integral_3, integral_4

def getPoly(D, N, integral_1, integral_2, integral_3, integral_4):
    """
    From the discriminant D of the real quadratic field and the unit specified by N
        N=4 - [4]-3[2]+2[1],
        N=prime - [N]-N[1],
    We compute the polynomial of the extension for this field. 
    Returns: Polynomial of the extension, 0 if there is an error in rational reconstruction of coefficients.
    """
    print("Discriminant: %d\n" %D)
    R.<x> = QQ[]
    K = pari.bnfinit(x^2-D)
    print("Narrow Class Group of quadratic field K:" + str(pari.bnfnarrow(K)[1]) + "\n")
    
    """Construct T, which is Q_p extended by sqrt{D}."""
    T.<d> = S.ext(x^2-D)
    reps = []
    beta = T.primitive_root_of_unity()
    
    """
    totalDeg tells us what to truncate from the precision.
    bqfgam_tPairs is the list of [BQF,Gamma_tau] pairs.
    """ 
    totalDeg = 0
    bqfgam_tPairs = getGamTau(D,N)
    
    """Calculate the Elliptic units."""
    load("Algo_4_3.sage")
    for pair in bqfgam_tPairs:
        bqf, gam_tau = pair[0], pair[1]
        print("Output for BQF %s:" %str(bqf))
        a = bqf(1,0)
        c = bqf(0,1)
        b = bqf(1,1) - a - c
        tau = T(-b/(2*a)) + T(d/(2*a))
        rep = u_alpha_tau(T,int(pair[1][0])+0,int(pair[1][2]/N)+0,p,N,tau,beta,integral_1,integral_2,integral_3,integral_4)
        totalDeg += abs(rep.valuation())+2
        reps.append(rep)
        reps.append(1/rep)
        
#     print("\nElliptic units u(alpha,taui) and u(alpha,taui)^-1:")
#     for rep in reps:
#         print(rep)
    
    """Ring for rational polynomial output"""
    K.<sqrtD> = NumberField(x^2-D)
    X = K['X'].gen()
    extension = K(0)
    
    try:
        output = T(1)
        for rep in reps:
            output *= (x-rep) 
#         print("\nCharacteristic Polynomial in Zp: %s" %str(output))
        exponents = output.exponents()
        output = output.coefficients()

        """
        Rational reconstruction of polynomial.
        -self.rational_reconstruction() is the default choice, but is inefficient as it 
            requires too many digits for reconstruction.
        -We opt for the following method:
            First check if the last 10 digits of the p-adic number is the same. If so, then we assume that this 
            is the recurring digit for the whole number. This works for the currently implemented units as the 
            denominator is p^k or 2*p^k.
            Otherwise, we use simplify(). This always returns a rational number, but it may not be accurate. 
        """
        for i in range(len(exponents)):
            print("\nCoefficients of %s:" %str(x^exponents[i]))
            
            coefficients = output[i].polynomial().coefficients()
            if len(coefficients) == 1:
                """Give 0 for coefficient of sqrt(D) """
                coefficients.append(0)
            t, t_d = coefficients[0], coefficients[1]
            t = t + O(p^(M-totalDeg))
#             print("p-adic value: %s" %str(t))
            
            """Rational reconstruction requires a much higher precision."""
#             t_rational = t.rational_reconstruction()    
            t_pAdicExpansion = list(t.expansion())
            if len(set(t_pAdicExpansion[-10:-1])) == 1:
                recurringDigit = t_pAdicExpansion[-1]
                while len(t_pAdicExpansion) < 2*M:
                    t_pAdicExpansion.append(recurringDigit)
                power = t.valuation()
                t_new = S(0)
                while len(t_pAdicExpansion) > 0:
                    t_new += t_pAdicExpansion.pop(0)*p^power
                    power += 1
                t_new += O(p^power)
                t_rational = simplify(QQ(t_new))
            else:
                t_rational = simplify(QQ(t))   
    
            if int(t_rational) == t_rational:
                print("Rational reconstruction for rational part of coefficient: %s" %str(t_rational))
            else:
                numerator = t_rational.numerator()
                denominator = t_rational.denominator()
                vP = denominator.ord(p)
                q = denominator/(p^vP)
                if q == 1:
                    print("Rational reconstruction for rational part of coefficient: %s/%s^%s" 
                          %(str(numerator),str(p),str(vP)))
                else: 
                    print("Rational reconstruction for rational part of coefficient: %s/%s*%s^%s "
                          %(str(numerator),str(q),str(p),str(vP)))
            extension += t_rational*X^(exponents[i])
            
            t_d = t_d + O(p^(M-totalDeg))
#             print("p-adic value: %s" %str(t_d))

            """Rational reconstruction requires a much higher precision."""
#             t_drational = t_d.rational_reconstruction() 
            t_dpAdicExpansion = list(t_d.expansion())
            if len(set(t_dpAdicExpansion[-10:-1])) == 1:
                recurringDigit = t_dpAdicExpansion[-1]
                while len(t_dpAdicExpansion) < 2*M:
                    t_dpAdicExpansion.append(recurringDigit)
                power = t_d.valuation()
                t_dnew = S(0)
                while len(t_dpAdicExpansion) > 0:
                    t_dnew += t_dpAdicExpansion.pop(0)*p^power
                    power += 1
                t_dnew += O(p^power)
                t_drational = simplify(QQ(t_dnew))
            else:
                t_drational = simplify(QQ(t_d))  

            if int(t_drational) == t_drational:
                print("Rational reconstruction for sqrtD part of coefficient: %s" %str(t_drational))
            else:
                numerator = t_drational.numerator()
                denominator = t_drational.denominator()
                vP = denominator.ord(p)
                q = denominator/(p^vP)
                if q == 1:
                    print("Rational reconstruction for sqrtD part of coefficient: %s/%s^%s" 
                          %(str(numerator),str(p),str(vP)))
                else: 
                    print("Rational reconstruction for sqrtD part of coefficient: %s/%s*%s^%s "
                          %(str(numerator),str(q),str(p),str(vP)))
            extension += t_drational*sqrtD*X^(exponents[i])
        
        """Print the polynomial, and check various properties of the extension."""
        print("\nField extension L/K by characteristic Polynomial for %d: %s\n" %(D,extension.factor()))
        L.<E> = K.extension(extension)
        print("Relative discriminant of L/K is trivial:", str(L.relative_discriminant()) == "Fractional ideal (1)")
        print("L/K is a relative galois extension:", L.is_galois_relative(), "\n")
#         print("Field automorphisms of L/K:", L.automorphisms())
        return extension

    except:
        print('Error in find polynomial or constructing field extension for: ', D)
        return 0


def getBQF(D,N):
    """
    Given discriminant D and N, return BQF's forms such that N|a, a>0. 
    """
    #BQF classes augmented so that N divides a
    BQF_N = []
    #Original classes of BQFs, where for two BQFs in the same class we take the one with positive leading coefficient.
    BQFs = []
    #Original BQFs.
    BQFsAll = BinaryQF_reduced_representatives(D, primitive_only=True)
    
    print("List of BQFs:")
    
    for BQF in BQFsAll:
        print(BQF)
        if BQF(1,0) > 0:
            BQFs.append(BQF)
    
    print("\nRepresentatives of classes of BQFs, augmented to N:")
    
    for Q in BQFs:
        """
        Using the primitive BQFs and the equivalence (a,b,c) ~ (-a,b,-c), we convert them to BQFs that are Heegner, where
        a>0, N|a, b is congruent to a fixed residue, and gcd(a,b,c) = 1.
        See Proposition 1.4 from Darmon, Henri. "Heegner points, Heegner cycles, and congruences." Elliptic curves and related topics. Vol. 4. 2007.
        """
        
        """We first fix a residue (mod 2N)."""
        residue = 0
        while (residue**2-D)%(4*N) != 0 and residue < 4*N:
            residue += 1
        if residue == 4*N:
            print("Error in finding residue satisfying given conditions.")
            raise 
        a = Q(1,0)
        b = Q(1,1)-Q(1,0)-Q(0,1)
        c = Q(0,1)
        if a%N==0 and (b-residue)%(2*N)==0:
            if Q(1,0) < 0:
                Q = BinaryQF(-a,b,-c)
            BQF_N.append(Q)
        else: 
            """Modify Q so that c is coprime to N"""
            """Warning: Currently implementation is only checked (and works) for N=4 and N prime."""
            if gcd(c,N)!=1:
                if gcd(a,N)==1:
                    Q = Q.matrix_action_right(matrix([[0, 1],[-1, 0]])) 
                    a = Q(1,0)
                    b = Q(1,1)-Q(1,0)-Q(0,1)
                    c = Q(0,1)
                elif gcd(b,N)==1:
                    Q = Q.matrix_action_right(matrix([[0, 1],[-1, 1]])) 
                    a = Q(1,0)
                    b = Q(1,1)-Q(1,0)-Q(0,1)
                    c = Q(0,1)
                else: 
                    raise NotImplemented
            """Convert to form where N|A and B=residue(mod 2N)"""
            if a%N == 0 and (b-residue)%(2*N)==0:
                """Found a valid form"""
                if a < 0:
                    Q = BinaryQF(-a,b,-c)
                BQF_N.append(Q)
            else:
                t = 0
                while (t*c-(residue-b)/2)%N != 0:
                    t += 1
                Q = Q.matrix_action_right(matrix([[1, 0],[t, 1]])) 
                a = Q(1,0)
                b = Q(1,1)-Q(1,0)-Q(0,1)
                c = Q(0,1)
                if a < 0:
                    Q = BinaryQF(-a,b,-c)
                BQF_N.append(Q)
    print("\n")
    return BQF_N

def gamTau(a, b, c): 
    """returns matrix {[A, B], [C, D]} as array [A, B, C, D]
    where ax^2 + bxy + cy^2 is our Binary Quadratic Form
    tau = (-b +sqrt(D))/2a"""
    from sympy.solvers.diophantine.diophantine  import diop_quadratic
    from sympy import symbols
    """Solve for minimum solution to Pell's equation."""
    x, y, = symbols("x, y", integer=True)
    Disc = b**2 - 4*a*c
    sols = diop_quadratic(x**2 - Disc*y**2 - 1, 0)
    first = sols.pop()
    u = abs(first[0])
    if(u + first[1]*sqrt(Disc) < 1 or u - first[1]*sqrt(Disc) < 0):
        v = - first[1]
    else:
        v = first[1]
        
    """Let our matrix Gamma_Tau be GamA, GamB, GamC, GamD
    solve u + v*sqrt(D) = GamC * tau + GamD"""
    GamC = 2*a*v
    GamD = u + b*v
    """solve GamA*tau + GamB = (GamC*tau + GamD) * tau. Easier to use GamC*tau + GamD = u+v*sqrt(D)"""
    GamA = u - b*v
    GamB = (v*Disc - (b**2)*v)/(2*a)

    print("GammaTau for BQF (%d,%d,%d): %s"%(a,b,c,str([GamA, GamB, GamC, GamD])))
    
    return [GamA, GamB, GamC, GamD]

def getGamTau(D,N):
    bqfs = getBQF(D,N)
    gamTaus = []
    for bqf in bqfs: 
        a = bqf(1,0)
        c = bqf(0,1)
        b = bqf(1,1) - a - c
        gamTaus.append([bqf,gamTau(a,b,c)])
    return gamTaus 

def getGoodDs(N,p,maxDisc):
    """Returns fundamental D discriminants <= maxDisc satisfying:
    Q(sqrt{D}) has no element of norm -1, prime p does not split, and N splits. Alternatively N|D also works."""
    goodDs = []
    
    """List of discriminants where the negative Pell's Equation has a solution."""
    if maxDisc <= 10000:
        badDs = [2,5,10,13,17,26,29,37,41,50,53,58,61,65,73,74,82,85,89,97,101,106,109,113,122,125,130,137,145,149,157,170,173,
             181,185,193,197,202,218,226,229,233,241,250,257,265,269,274,277,281,290,293,298,313,314,317,325,337,338,346,349,
             353,362,365,370,373,389,394,397,401,409,421,425,433,442,445,449,457,458,461,481,485,493,509,521,530,533,538,541,
             554,557,565,569,577,586,593,601,610,613,617,626,629,634,641,653,661,673,677,685,697,698,701,709,730,733,746,754,
             757,761,769,773,778,785,794,797,809,818,821,829,842,845,853,857,865,877,881,901,914,922,925,929,937,941,949,953,
             962,965,970,977,985,986,997,1009,1010,1013,1018,1021,1025,1033,1037,1042,1049,1061,1066,1069,1073,1082,1090,1093,
             1097,1105,1109,1114,1117,1129,1130,1138,1145,1153,1157,1165,1181,1189,1193,1201,1213,1217,1226,1229,1237,1241,1249,
             1250,1258,1261,1277,1285,1289,1297,1301,1306,1313,1321,1322,1325,1354,1361,1370,1373,1378,1381,1385,1402,1409,1417,
             1418,1429,1433,1445,1450,1453,1465,1466,1481,1489,1490,1493,1514,1522,1525,1546,1549,1553,1565,1570,1585,1586,1594,
             1597,1601,1609,1613,1618,1621,1625,1637,1642,1649,1657,1658,1669,1682,1685,1693,1697,1706,1709,1714,1721,1730,1733,
             1741,1745,1753,1754,1765,1769,1777,1781,1789,1801,1810,1825,1850,1853,1861,1865,1873,1877,1882,1889,1898,1901,1906,
             1913,1921,1930,1933,1937,1949,1970,1973,1985,1993,1994,1997,2017,2026,2029,2042,2050,2053,2069,2074,2081,2089,2113,
             2117,2122,2125,2129,2137,2138,2141,2146,2153,2161,2165,2173,2186,2197,2210,2213,2218,2221,2234,2237,2249,2257,2258,
             2269,2273,2281,2285,2290,2293,2297,2305,2309,2314,2330,2333,2341,2357,2362,2377,2378,2381,2389,2393,2402,2405,2417,
             2425,2426,2437,2441,2458,2465,2473,2474,2477,2482,2501,2509,2521,2522,2545,2549,2554,2557,2561,2570,2581,2593,2602,
             2605,2609,2617,2621,2626,2633,2642,2657,2665,2677,2689,2690,2693,2705,2713,2725,2729,2738,2741,2746,2749,2753,2762,
             2770,2777,2785,2789,2797,2801,2810,2813,2825,2833,2834,2837,2857,2858,2861,2873,2885,2890,2897,2906,2909,2917,2929,
             2930,2941,2953,2957,2965,2969,2977,2986,3001,3026,3029,3037,3041,3049,3050,3061,3065,3074,3077,3085,3089,3098,3109,
             3121,3125,3130,3133,3137,3145,3161,3169,3170,3181,3194,3202,3209,3217,3221,3226,3229,3233,3242,3250,3253,3257,3265,
             3274,3277,3281,3293,3301,3313,3314,3329,3338,3341,3349,3361,3365,3370,3373,3385,3386,3389,3413,3418,3425,3433,3434,
             3445,3449,3457,3461,3466,3469,3482,3485,3490,3517,3529,3530,3533,3538,3541,3545,3554,3557,3562,3578,3581,3589,3593,
             3601,3613,3617,3625,3637,3649,3653,3665,3673,3677,3697,3701,3706,3709,3722,3725,3730,3733,3754,3761,3769,3770,3785,
             3793,3797,3802,3809,3821,3833,3845,3853,3865,3866,3869,3874,3877,3881,3889,3890,3898,3917,3922,3925,3929,3946,3961,
             3970,3973,3977,3985,3986,3989,3994,4001,4013,4021,4033,4045,4049,4057,4058,4073,4082,4093,4097,4106,4121,4129,4133,
             4138,4141,4153,4157,4177,4181,4201,4210,4217,4226,4229,4234,4241,4250,4253,4261,4265,4273,4274,4282,4285,4289,4297,
             4306,4325,4330,4337,4346,4349,4357,4373,4385,4394,4397,4409,4421,4426,4441,4442,4457,4469,4474,4481,4490,4493,4498,
             4505,4513,4514,4517,4537,4538,4549,4553,4561,4570,4573,4586,4589,4594,4597,4610,4618,4621,4625,4637,4649,4657,4666,
             4673,4682,4685,4706,4709,4714,4721,4729,4733,4745,4754,4762,4765,4777,4778,4789,4793,4801,4813,4817,4825,4861,4874,
             4877,4885,4889,4901,4909,4913,4925,4933,4937,4954,4957,4969,4973,4985,4993,5009,5018,5021,5042,5045,5050,5065,5077,
             5081,5090,5098,5101,5113,5114,5122,5153,5161,5162,5165,5185,5189,5197,5209,5210,5213,5233,5234,5237,5242,5245,5261,
             5266,5273,5281,5297,5305,5309,5317,5321,5330,5333,5353,5354,5365,5381,5386,5389,5393,5410,5413,5417,5426,5429,5437,
             5441,5449,5450,5458,5465,5473,5477,5482,5485,5498,5501,5521,5545,5554,5557,5569,5570,5573,5578,5581,5585,5594,5597,
             5617,5618,5626,5629,5641,5653,5657,5669,5674,5689,5693,5701,5713,5717,5722,5729,5737,5741,5749,5765,5770,5777,5785,
             5801,5813,5818,5821,5825,5834,5837,5849,5857,5858,5861,5869,5881,5882,5897,5906,5914,5930,5933,5941,5953,5954,5965,
             5981,5993,6002,6010,6025,6029,6037,6053,6065,6073,6074,6085,6089,6101,6109,6113,6121,6122,6130,6133,6145,6154,6161,
             6170,6173,6178,6185,6197,6205,6217,6218,6221,6229,6242,6245,6250,6253,6257,6266,6269,6274,6277,6301,6305,6317,6322,
             6329,6337,6353,6361,6362,6373,6385,6389,6397,6401,6409,6418,6421,6425,6437,6442,6445,6449,6458,6466,6469,6473,6481,
             6485,6506,6521,6529,6530,6553,6554,6562,6565,6569,6577,6581,6586,6602,6610,6617,6625,6626,6637,6641,6649,6653,6661,
             6673,6682,6689,6698,6701,6709,6725,6730,6733,6737,6746,6749,6757,6761,6770,6778,6781,6793,6817,6826,6829,6833,6841,
             6845,6850,6857,6865,6866,6869,6890,6893,6917,6922,6925,6929,6938,6949,6953,6961,6970,6977,6989,6994,6997,7001,7010,
             7013,7025,7033,7034,7045,7057,7066,7069,7082,7085,7090,7093,7109,7114,7121,7129,7141,7162,7165,7177,7178,7186,7193,
             7202,7213,7226,7229,7237,7241,7250,7253,7261,7265,7274,7297,7306,7309,7321,7325,7330,7333,7345,7349,7354,7369,7373,
             7393,7397,7402,7417,7418,7421,7433,7442,7450,7457,7465,7466,7474,7477,7481,7489,7501,7517,7522,7529,7537,7538,7541,
             7549,7561,7570,7573,7577,7585,7589,7594,7610,7618,7621,7633,7642,7649,7669,7673,7681,7706,7709,7717,7730,7738,7741,
             7745,7753,7754,7757,7762,7765,7789,7793,7801,7817,7825,7829,7834,7837,7841,7850,7853,7858,7873,7877,7897,7901,7913,
             7922,7925,7930,7933,7937,7946,7949,7969,7970,7978,7985,7993,8005,8009,8017,8021,8026,8042,8053,8065,8066,8069,8077,
             8081,8089,8093,8101,8117,8125,8146,8161,8177,8185,8186,8209,8210,8221,8233,8237,8242,8245,8249,8266,8269,8273,8282,
             8285,8290,8293,8297,8314,8317,8321,8329,8345,8353,8354,8362,8369,8377,8381,8389,8410,8425,8429,8450,8458,8461,8465,
             8485,8497,8501,8506,8513,8521,8522,8530,8537,8545,8570,8573,8581,8593,8597,8609,8629,8641,8642,8650,8653,8665,8669,
             8677,8681,8689,8693,8698,8713,8714,8737,8741,8746,8753,8761,8765,8770,8794,8821,8825,8837,8842,8849,8857,8861,8882,
             8885,8893,8905,8906,8917,8929,8933,8938,8941,8945,8969,8986,8989,9001,9010,9013,9026,9029,9034,9041,9049,9050,9061,
             9074,9089,9098,9106,9109,9125,9133,9137,9146,9157,9161,9169,9173,9178,9181,9193,9194,9197,9209,9217,9221,9241,9242,
             9250,9257,9265,9274,9277,9281,9290,9293,9298,9305,9325,9337,9341,9349,9365,9370,9377,9385,9389,9397,9410,9413,9418,
             9421,9433,9434,9437,9442,9458,9461,9466,9473,9497,9509,9521,9529,9530,9533,9553,9565,9577,9578,9586,9593,9601,9605,
             9613,9626,9629,9634,9649,9661,9665,9673,9677,9685,9689,9697,9698,9721,9722,9725,9733,9749,9754,9769,9770,9773,9778,
             9781,9802,9805,9817,9818,9829,9833,9857,9865,9866,9881,9893,9901,9914,9925,9929,9938,9941,9946,9949,9953,9965,9970,
             9973,9985,9997]  
    else:
        badDs = []
        from sympy.solvers.diophantine.diophantine  import diop_quadratic
        from sympy import symbols
        """Check of solvability of negative Pell's equation."""
        x, y, = symbols("x, y", integer=True)
        for D in range(1,maxDisc+1):
            sols = diop_quadratic(x**2 - D*y**2 + 1, 0)
            if len(sols) > 0:
                badDs.append(D)
 
    for D in range(maxDisc+1):
        if N == 4:
            """We need the discriminant to be 1 Mod 8 and square free to eliminante non-fundamental discriminants."""
            if (D+0).is_squarefree() and D%8 == 1 and kronecker(D,p) == -1:
                if D not in badDs:
                    goodDs.append(D)
        else:
            if (D+0).is_squarefree() and D%4 == 1 and kronecker(D,p) == -1 and (kronecker(D,N) == 1 or D%N == 0):
                if D not in badDs:
                    goodDs.append(D)
            elif D%4 == 0 and (int(D/4)+0).is_squarefree() and kronecker(D,p) == -1 and (kronecker(D,N) == 1 or D%N == 0):
                if (int(D/4) not in badDs) and (int(D/4) not in goodDs):
                    goodDs.append(D)
    print("Output valid discriminants from 1 to %d, for p = %d and N = %d:" %(maxDisc,p,N))    
    print(goodDs) 
    return goodDs


