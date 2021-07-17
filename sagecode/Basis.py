"""N prime. The basis used will be [inf]-[1/N], [inf]-[2/N], ... , [inf]-[(N-1)/N]
N = 4. Basis used [inf]-[1/N]
"""
class Gamma:
    def __init__(self, a, c):      #given numerator and denominator, construct rest of matrix
        """Attributes defined for Gamma:
        a,b,c,d: The matrix coefficients
        num,denom: Numerator and denomnator that generates Gamma by 4.1. Denominator recorded has been divided by N.  
        coeff: +1 or -1, + if c*d >0, - if c*d<0
        x: The numerator in the basis element [inf]-[x/N]
        """
        if c < 0:
            a, c = -a, -c
        self.num = int(a)
        self.denom = int(c)
        
        if self.denom == 0:					#case of c = 0 is the infinity case
            """If our element is null and we want to return the gammas"""
            self.a = 0
            self.b = 0
            self.c = 0
            self.d = 0				#d will never equal 0 unless its infinity, so this is what we can use to check infinity
            self.coeff = 0
            self.x = 0
        elif N != 4 and self.denom == 1 and self.num >= 1 and self.num <= N-1:
            """If the element we pass in is a basis element"""
            self.a = 1
            self.b = 0
            self.c = 0
            self.d = 1
            self.coeff = 1
            self.x = a
        else:
            base = N*c
            d = ""
   
            #If N=4, change to =1, else =x. 
            if N == 4:
                d = (a+0).inverse_mod(base)
                if d > base/2: 
                    d -= base
                b = (a*d - 1)/(N*c)
                self.x = 1
            else:
                d = a.inverse_mod(base) - base
                b = (a*d - 1)/(base)
                self.x = max(int(abs(d/c)), 1)								#x is not the unique choice 
                
            if c*d > 0:												#second matrix in 4.1
                self.a = int(a-N*b)
                self.b = int(b)
                self.c = int(N*c-N*d)
                self.d = int(d)
                self.coeff = 1
            else:													#first matrix
                self.a = int(a)
                self.b = int(b)
                self.c = int(N*c)
                self.d = int(d)
                self.coeff = -1
    
    def __str__(self):     #for print statements
        return (str(self.a) + " " + str(self.b) + " " + str(self.c) + " " + str(self.d))
    
def getGammas(num, denom):
    """returns an array of every gamma used, which gives our basis"""
    a = []							
    gam = Gamma(num, (denom/N))			#constructs the first gamma given the fraction
    a.append(gam)
    return getGammasHelper(gam, a)		#sends to the recursive function
    
def getGammasHelper(prev, gammas):		
    """Helper method to the above. Recursively calls itself as done in algorithm 4.1"""
    if prev.d == 0:					#if infinity, we are done
        return gammas				
    else:						
        if prev.coeff == 1:
            num = prev.a  					#create new fraction
            denom = prev.c
            gam = Gamma(num, (denom/N))					#create new matrix
            gammas.append(gam)
            return getGammasHelper(gam, gammas)			#run it again
        else:
            num = prev.a*prev.x + prev.b*N
            denom = prev.c*prev.x + prev.d*N
            gam = Gamma(num, (denom/N))
            gammas.append(gam)
            return getGammasHelper(gam, gammas)


