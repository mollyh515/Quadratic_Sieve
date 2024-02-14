import math

def is_prime(n):      #checks if n is prime
    if n<=1:
        return False
    for i in range(2,int(math.sqrt(n))+1):
        if n%i == 0:
            return False
    return True

def LJ(a,p):
    a=a%p   #now 0<=m<n
    t=1
    while a!=0:
        
        while a%2==0:
            a=int(a//2)
            if p%8==3 or p%8==5:
                t=-t
        
        a,p=p,a
        if a%4==3 and p%4==3:
            t=-t
        a=a%p
    
    if p==1:
        return t
    else:
        return 0

def QC_357(a,p):
    #faster for p=3,5,7(mod 8)
    if p%8==3 or p%8==7:
        x=pow(a,int((p+1)//4),p)   #computes a^((p+1)/4)(mod p)
        return x
    
    if p%8==5:
        x=pow(a,int((p+3)//8),p)   #computes a^((p+3)/8)(mod p)
        x_squared=pow(x,2,p)  #computes x^2(mod p)
        if x_squared%p==a:
            return x
        else:
            x=(x*pow(2,int((p-1)//4),p))%p  #computes x*2^((p-1)/4)(mod p)
            return x
        
def QC_18(a,p):
    #works for any odd prime
    b=2
    while pow(b,int((p-1)//2),p)!=p-1:   #computes Legendre of (b/p) until (b/p)=-1 since (b/p)=b^((p-1)/2)(mod p) by Euler's Criterion
        b+=1
    
    u=p-1
    s=0
    while u%2==0:
        s+=1
        u=int(u//2)    #must be int or else we will have .0 at the end


    d=pow(a,u,p)
    f=pow(b,u,p)

    m=0
    for i in range(s):   #stops at s-1
        g=pow(d*(f**m),2**(s-1-i),p)
        if g%p==p-1:   #g=p-1(mod p)=-1(mod p)
            m=m+(2**i)
    
    x=(pow(a,int((u+1)//2),p)*pow(f,int(m//2), p)) % p
    return x

def qr_modp(n,p_vals):
    qr_modp=[]
    for p in p_vals:
        if LJ(n,p)==1:        #QR
            qr_modp.append(p)
    return qr_modp


def Gaussian_Elimination(AT):
    #takes in transpose of matrix A and returns null space solution to λA=0
    pivot_cols=[False]*len(AT[0]) #True if pivot column in A^T

    for i in range(len(AT)):  #iterate through rows
        for j in range(len(AT[0])): #iterate through columns
            if AT[i][j]==1:   #pivot
                pivot_cols[j]=True
                for k in range(len(AT)):
                    if k!=i and AT[k][j]==1:
                        for idx in range(len(AT[k])):
                            AT[k][idx]=(AT[k][idx]+AT[i][idx])%2 #Add Row i to Row 2
                break 

    pivot_positions=[row.index(1) if 1 in row else -1 for row in AT] #finds row indices of pivots

    AT=[row for _, row in sorted(zip(pivot_positions,AT)) if _!=-1]+[row for _, row in sorted(zip(pivot_positions,AT)) if _==-1] #sorts the pivots

    null_space=[]

    for j in range(len(AT[0])): #construct the null space of A^T
        if not pivot_cols[j]:   #not pivot column (free variable)
            null_vector=[0]*len(AT[0])
            null_vector[j]=1
            for i in range(len(AT)):
                for k in range(len(AT[0])):
                    if pivot_cols[k] and AT[i][k]==1 and AT[i][j]==1:  #current fv in current row, current pivot in current row
                        null_vector[k]=1  #only one pivot per row
                        break
            null_space.append(null_vector)


    return(null_space)



def QS(n):
    #1 Initialization
    B=math.sqrt(math.e**(math.sqrt(math.log(n)*math.log(math.log(n))))) 
    p_vals=[p for p in range(3, int(B)+1) if is_prime(p)]
    odd_primes=qr_modp(n,p_vals)
    factor_base=[-1, 2]+odd_primes
    K=len(odd_primes)+1
    print(K)
    sol=[]
    for pk in odd_primes:  #excluding  2
        if pk%8==3 or pk%8==5 or pk%8==7:
            tk=QC_357(n,pk)
            sol.append((pk,tk,-tk))
        else:
            tk=QC_18(n,pk)
            sol.append((pk,tk,-tk))

    #2 Sieveing
    B_smooth=[]
    M1=[]
    N=math.isqrt(n)+1
    j=0
    #while len(B_smooth)<K+8:   #n4 and n5
    while len(B_smooth)<K+2:  #n1, n2, and n3
        x=N+j
        #x^2-n
        z=pow(x,2)-n
        w=z
        vector=[x,z]
        exps=[0]*len(factor_base)
        for idx, pk in enumerate(factor_base):  #including 2
            exp=0
            if w>0 and pk==-1:    #otherwise we will keep dividing a positive number by -1
                continue
            if w<0 and pk==-1:   #odd numbers have only one factor of -1
                exps[idx]=1
                w=w//pk
                continue
            elif w%pk==0:           #a factor
                while w%pk==0:
                    exp+=1
                    w=w//pk
            exps[idx]=exp
            if w==1:       #completely factored so B-smooth since we are only checking primes less than or equal to B
                B_smooth.append(z)   #z is w before factoring
                vector+=exps  
                M1.append(vector)    #only want those that are B-smooth in the matrix
                break
        j+=1

    #3 Linear Algebra
    M2=[[j%2 for j in i[2:]] for i in M1] #Constructs (K+2)x(K+1) matrix formed by the exponent vectors of M1 reduced modulo 2
    M2T=[[M2[j][i] for j in range(len(M2))] for i in range(len(M2[0]))] #(λM)^T=(M^T)(λ^T)
    x_vals=[i[0] for i in M1]
    print(x_vals)
    possible_solns=Gaussian_Elimination(M2T)
    if len(possible_solns)==1:
        l=possible_solns
    else:
        for soln in possible_solns:
            y=1
            x=1
            exp_vector1=[0]*len(M2[0])
            for i in range(len(soln)):
                if soln[i]==1:
                    x*=x_vals[i]
            for i in range(len(soln)):
                if soln[i]==1:
                    for j in range(len(exp_vector1)):
                        exp_vector1[j]+=M1[i][j+2]
            exp_vector2=[i//2 for i in exp_vector1]
            for idx,pk in enumerate(factor_base):
                if exp_vector2[idx]==0:
                    continue
                else:
                    y*=pk**exp_vector2[idx]
            y=abs(y)
            x=pow(x,1,n)   #modulo n
            y=pow(y,1,n)   #modulo n
            m1=math.gcd(x+y,n)
            m2=math.gcd(x-y,n)
            if (m1!=1 and m1!=n) and (m2!=1 and m2!=n):
                l=soln
                break
    print('x=',x)
    print('y=',y)  
    print('l=',l)
    print('gcd(x+y,n)=',m1)
    print('gcd(x-y,n)=',m2)

n1=3215031751
n2=9912409831
n3=37038381852397
#n4=341550071728321
#n5=31868712526338419047
print(QS(n1))
print(QS(n2))
print(QS(n3))
#print(QS(n4))
#print(QS(n5))