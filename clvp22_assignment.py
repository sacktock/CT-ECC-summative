import copy

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

##############
#My Functions#
##############

#function parityCheckMatrixGenerator
#input: r
#output: H, a matrix with dimensions 2^r-1,r for checking codewords
def parityCheckMatrixGenerator(r):
    H =[]
    for i in range(1,2**r):
        H.append(decimalToVector(i,r))
    return H

#function vectorToDecimal
#input: v, a binary vector
#output: t, the decimal representation of v
def vectorToDecimal(v):
    t=0
    L = len(v)
    for i in range(0,L):
        t += v[i]*2**(L-i-1)
    return t

#function matrix multiplication
#input: A and B, 2 matrices
#output: matrix C in binary space
def matrixMultiplication(A,B):
    if len(A) != len(B):
        return []
    
    C=[]
    for i in range(0,len(A[0])):
        t=0
        for j in range(0,len(B)):
            t+=A[j][i]*B[j]
        C.append(t%2)

    return C
            
    
#function message
#input: a, a binary vector of any positive length
#output: m, a binary vector of length 2^r - r - 1
def message(a):
    if type(a) != list:
        return []
    
    L = len(a)
    r=2
    
    while 2**r -2*r -1 < L:
        r+=1
        
    m = decimalToVector(L,r)
    m.extend(a)
    k = 2**r - r - 1

    while len(m) < k:
        m.append(0)

    return m

#function hammingEncoder
#input: m, a binary vector of length 2^r - r - 1
#output: c, a codeword binary vector of length 2^r -1
def hammingEncoder(m):
    if type(m) != list:
        return []
    
    r =2
    L= len(m)
    
    while 2**r-r-1 != L:
        r+=1
        if 2**r-r-1 > L:
            return []
        
    G = hammingGeneratorMatrix(r)
    
    c =matrixMultiplication(G,m)
    
    #for j in range(0,2**r-1):
     #   t=0
      #  for i in range(0,L):
       #     t += G[i][j]*m[i]
        #c.append(t%2)
    
    return c

#function hammingDecoder
#input: v, a binary vector of length 2^r - 1
#output: c, a codeword binary vector of length 2^r - 1
def hammingDecoder(v):
    if type(v) != list:
        return []
    
    r=2
    L= len(v)

    while 2**r-1 != L:
        r+=1
        if 2**r-1 > L:
            return []
        
    H = parityCheckMatrixGenerator(r)

    c=copy.deepcopy(v)
    i=0
    
    while True: #method 2 - local search
        z= matrixMultiplication(H,c)
        #for j in range(0,r):
        #    t=0
         #   for k in range(0,L):
          #      t += H[k][j]*c[k]
           # z.append(t%2)
        if sum(z) == 0:
            break
        elif i >= L:
            return []
        else:
            c=copy.deepcopy(v)
            c[i] = (c[i] + 1) % 2
        i+=1

    return c

#function messageFromCodeword
#input: c, a codeword binary vector of length 2^r - 1
#output: m, a binary vector of length 2^r - r - 1
def messageFromCodeword(c):
    if type(c) != list:
        return []
    
    r=2
    L=len(c)
    
    while 2**r-1 != L:
        r+=1
        if 2**r-1 > L:
            return []
    m=copy.copy(c)
    
    for i in range(0,r):
        del m[2**i-1]
        
    return m

#function dataFromMessage
#input: m, a binary vector of length 2^r - r - 1
#output: a, a binary vector of any length - raw data
def dataFromMessage(m):
    if type(m) != list:
        return []
    
    r=2
    L=len(m)
    
    while 2**r -r-1 != L:
        r+=1
        if 2**r -r -1 > L:
            return []
    
    d = vectorToDecimal(m[:r])

    if d+r > len(m):
        a = []
    elif sum(m[r+d:]) == 0:
        a = m[r:r+d]
    else:
        a = []

    return a
        
def repetitionEncoder(m,n):
    if type(m) != list:
        return []
    
    if type(n) != int:
        return []
    
    if len(m) != 1:
        return []
    
    return m*n

def repetitionDecoder(v):
    k = sum(v)
    P = len(v)/2
    
    if k < P:
        return [0]
    elif k > P:
        return [1]
    else:
        return []
            

    
    
