import numpy as np

def GetAlpha(N=5, ReturnAll=False):
    Data = {}
    AlphaData = {}

    F = open("P30-%d.txt"%(N))

    alpha = np.array([0,0.05,0.1,0.15,0.2,0.4,0.6,0.8,1.0])
    for L in F:
        X = L.split()
        if len(X)<3: continue

        ID = X[0]
        Ealpha = np.array([float(x) for x in X[1:]])

        if ReturnAll: AlphaData[ID]=(alpha,Ealpha)

        pLin = np.polyfit(alpha, Ealpha, 1)
        aLin = np.roots(pLin)[0]

        if aLin>=0. and aLin<=1.:
            ii = np.abs(Ealpha).argmin()
            i0 = max(ii-1,0)
            i1 = i0+3
            if i1>len(alpha):
                i1 = len(alpha)
                i0 = i1-3
            p = np.polyfit(alpha[i0:i1],Ealpha[i0:i1],2)
        elif aLin<0.:
            p = np.polyfit(alpha[:5], Ealpha[:5], 2)
        else:
            p = np.polyfit(alpha[-3:], Ealpha[-3:], 1)
        q = np.roots(p)

        if np.max(np.imag(q))>0.:
            q = [-1000.,]

        if np.abs(np.sum(np.abs(np.imag(q))**2))>0.:
            p = np.polyfit(alpha[-2:], Ealpha[-2:], 1)

        if len(q)==1 or (np.abs(q[0])<np.abs(q[1])):
            a0 = q[0]
        else:
            a0 = q[1]

        Data[ID] = a0

    if ReturnAll:
        return Data, AlphaData
    return Data

