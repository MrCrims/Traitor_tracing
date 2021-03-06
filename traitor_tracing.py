# coding : utf-8
'''
@creat_time = 2022/3/9,21:20
@auther = MrCrimson
Emal : mrcrimson@163.com
'''
import gmpy2
from Crypto.Util import number
import sympy,random,math
import numpy as np
#Q1：有限域上线性空间的基？
gamma = [[]]# gamma matrix store codewords
k = 5
l = 15 # l>=2k+2
q = number.getPrime(64)

def gamma_generation():
    return

def key_generation():
    g = sympy.primitive_root(q) #  g is a generator of Group q

    PBK = []
    all_r = []
    y = 1
    theta_up = 0
    # Gen Public Key:
    # for i in range(1,2k) random choose r in Fq ,compute h = g^r
    for i in range(2*k):
        r = random.randint(0,q-1) # as q is quite big I think it's impossible to meet 2 same r
        alpha = random.randint(0,q-1)
        all_r.append(r)
        theta_up += r*alpha % q
        h = pow(g,r,q)
        PBK.append(h)
        y = y*pow(h,alpha) % q
    # y = π h^α
    PBK.append(y)
    PRK_all = []

    for i in range(l):
        theta_down = 0
        for j in range(2*k):
            theta_down += all_r[j]*gamma[i][j] % q
        theta = theta_up * gmpy2.invert(theta_down,q) % q
        PRK_all.append(theta)
    return PBK,PRK_all

def Encryption(message,PBK):
    y = PBK[-1]
    a = random.random(0,q-1)
    C = []
    C.append(message*pow(y,a,q)%q)
    for h in PBK[:-1]:
        C.append(pow(h,a,q))
    return C

def Decryption(C,PRK,id):
    H = C[1:]
    U = 1
    for i in range(2*k):
        U = U * pow(H[i],gamma[id][i],q)
    message = C[0] * gmpy2.invert(pow(U,PRK,q),q) % q
    return message

def generate_gamma():
    # generate A matrix
    A = []
    for i in range(l-2*k):
        row = []
        for j in range(1,l+1):
            row.append(pow(j,i,q))
        A.append(row)
    # use schmi-gram process to get a set of basis
    

    print(np.linalg.matrix_rank(A))
    return

# Minimal Access Black box tracing against arbitrary pirates

def mini_access_black_box_tracing(Decoder,T_suspect,PBK,all_random,H):
    # Decoder is a pirate decoder
    # T_suspect is a set of PRK suspected for creating Decoder
    # all_random is a set including all random bits used in key generation
    def compute_P(lamda_,q,T_suspect,H):
        def make_CT(PBK,T_suspect,H):
            CT = []
            for i in range(len(T_suspect)):
                CT_i = []
                S = random.randint(q)
                for d in T_suspect[:i]:
                    W = 1
                    for j in range(2*k):
                        W = W*pow(H[j],d[j],q) % q
                    CT_i.append([S,W,H])
                CT.append(CT_i)
            return CT
        CT = make_CT(PBK,T_suspect,H)
        P = []
        for i in range(k):
            CW_i = CT[i]
            count = 0
            for j in range(1,lamda_):
                C = CW_i[random.randint(len((CW_i)))]
                if Decoder.decode(C):
                    count += 1
            p = count/lamda_
            P.append(p)
        return P
    epsilon = 0.1
    lamda_ = 64*k*math.log(k/epsilon)
    P = compute_P(lamda_,q,T_suspect,H)
    if P[k]<0.75:
        print("No guilty!")
        return False
    for j in range(len(P)):
        if math.fabs(P[j+1]-P[j])>1/(2*k):
            print("No."+str(j)+"is guilty")
            return j
if __name__ == '__main__':
    generate_gamma()