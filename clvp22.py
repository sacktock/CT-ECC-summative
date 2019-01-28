#function HammingG
#input: a number r
#output: G, the generator matrix of the (2^r-1,2^r-r-1) Hamming code
def hammingGeneratorMatrix(r):
    n = 2**r-1
    
    #construct permutation pi
    pi = []
    for i in range(r):
        pi.append(2**(r-i-1))
    for j in range(1,r):
        for k in range(2**j+1,2**(j+1)):
            pi.append(k)

    #construct rho = pi^(-1)
    rho = []
    for i in range(n):
        rho.append(pi.index(i+1))

    #construct H'
    H = []
    for i in range(r,n):
        H.append(decimalToVector(pi[i],r))

    #construct G'
    GG = [list(i) for i in zip(*H)]
    for i in range(n-r):
        GG.append(decimalToVector(2**(n-r-i-1),n-r))

    #apply rho to get Gtranpose
    G = []
    for i in range(n):
        G.append(GG[rho[i]])

    #transpose    
    G = [list(i) for i in zip(*G)]

    return G


#function decimalToVector
#input: numbers n and r (0 <= n<2**r)
#output: a string v of r bits representing n
def decimalToVector(n,r): 
    v = []
    for s in range(r):
        v.insert(0,n%2)
        n //= 2
    return v

def message(a):
    L = len(a)
    r=2
    
    while 2**r -2*r -1 < L:
        r+=1
        
    m = decimalToVector(L,r)
    m.extend(a)
    k = 2**r - r - 1

    while len(m) < k:
        m.append(0)

    print(r)
    return m

def hammingEncoder(m):
    r =2
    L= len(m)
    c =[]
    while 2**r-r-1 != L:
        r+=1
        if 2**r-r-1 > L:
            return c
    G = hammingGeneratorMatrix(r)
    for j in range(0,2**r-1):
        t=0
        for i in range(0,L):
            t += G[i][j]*m[i]
        c.append(t%2)

    return c


def repetitionEncoder(m,n):
    if type(m) != list:
        return
    if type(n) != int:
        return
    if len(m) != 1:
        return
    return m*n

def repetitionDecoder(v):
    k = sum(v)
    P = (len(v)+1) //2
    if k < P:
        return [0]
    elif k > P:
        return [1]
    else:
        return []
            

    
    
