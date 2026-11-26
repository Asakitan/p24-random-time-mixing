"""Verify the exact P24 torus simplex (standard library, Python >= 3.9)."""
import json
from math import comb
from pathlib import Path

F = (-1, -5, 25, 5, -50, -1, 35, 0, -10, 0, 1)  # ascending powers
P = 65537  # prime

def get_modes():
    # Check the substitution z^10 F(z+z^-1)=Phi_50(z).
    phi = [0]*21
    for k, a in enumerate(F):
        for j in range(k+1):
            phi[10+k-2*j] += a*comb(k,j)
    assert phi == [(-1)**(k//5) if k%5==0 else 0 for k in range(21)]
    h = [(2,)+(0,)*9, (0,1)+(0,)*8]
    for r in range(2,25):
        a = (0,)+h[-1]
        h.append(tuple(a[i]-h[-2][i]-a[10]*F[i] for i in range(10)))
    return h[1:]

def determinant_mod(A):
    a=[[e%P for e in row] for row in A]
    result=1
    for k in range(len(a)):
        i=next((i for i in range(k,len(a)) if a[i][k]),None)
        if i is None:
            return 0
        if k!=i:
            a[k],a[i]=a[i],a[k]
            result=-result
        pivot=a[k][k]
        result=(result*pivot)%P
        inv=pow(pivot,-1,P)
        for i in range(k+1,len(a)):
            t=a[i][k]*inv%P
            if t:
                for j in range(k+1,len(a)):
                    a[i][j]=(a[i][j]-t*a[k][j])%P
                a[i][k]=0
    return result%P

def verify():
    dec = json.loads(Path(__file__).with_name("certificate.json").read_text())
    G, u, D = dec["phase"], dec["numerator"], dec["denominator"]
    assert len(G)==len(u)==145 and all(len(g)==10 and all(isinstance(v,int) and 0<=v<6 for v in g) for g in G)
    assert all(isinstance(z,int) and z>0 for z in u) and sum(u)==D
    v=get_modes()
    targets={}
    for i in range(24):
        for j in range(i):
            w=tuple(a-b for a,b in zip(v[i],v[j]))
            w=min(w,tuple(-a for a in w))
            d=-int((i-j)%2==0)
            assert w not in targets or targets[w]==d
            targets[w]=d
    assert len(targets)==144
    A=[[1]*145]
    for w in sorted(targets):
        row=[(2,1,-1,-2,-1,1)[sum(a*b for a,b in zip(w,g))%6] for g in G]
        assert 12*sum(a*b for a,b in zip(row,u)) == targets[w]*D
        A.append(row)
    det=determinant_mod(A)
    assert det!=0, "simplex is singular modulo chosen prime"
    print(f"P24: 276 mode pairs, 144 gaps, 145 positive weights; det mod {P} = {det}.")

if __name__=="__main__":
    verify()
