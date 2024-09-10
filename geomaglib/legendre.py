import math
import numpy as np
def Flattened_Chaos_Legendre1(nmax, theta, epsilon = 1e-2):
    """
    Takes nmax and colatitude (degrees) as inputs
    Outputs a 2 dimensional numpy array which contains the associated
    legendre polynomials (Pnm) and the respective derivatives (dPnm).
    In these respective arrays the values are arranged in order major order:
    i.e.
    P10,P20,P30...Pnmax0, P11,P21,P31...Pnmax1

    """

    costh = np.cos(np.radians(theta))
    sinth = np.sqrt(1-costh*costh)

    # Pnm = np.zeros((nmax+1, nmax+2) + costh.shape)
    # print(np.shape(costh))
    # dPnm = np.zeros((int((nmax+1)*(nmax+2)/2) ,len(costh)))
    Pnm = []
    dPnm = []
    Pnm.append(np.ones(len(costh)))  # is copied into trailing dimensions
    dPnm.append(np.zeros(len(costh)))


    rootn = np.sqrt(np.arange(2 * nmax**2 + 1))

    # Recursion relations after Langel "The Main Field" (1987),
    # eq. (27) and Table 2 (p. 256)
    for m in range(nmax):
        if(m == 1):
            Pnm.append(sinth)
            dPnm.append(costh)
        # Buffer normalization factor for [n,n] to [n,n-1]
        c2 = rootn[m + m + 1]
        #Pnm_tmp holds Pnm[m,m] with normalization factor
        Pnm_tmp = c2 * Pnm[-1]
        #Pnm[m,m-1] = Pnm[m-1,m-1]*cos*c2
        Pnm.append(costh * Pnm_tmp)
        #Hold this for the derivative of diagonal later
        dPnm_diag_tmp = Pnm[-1]
        #Derivative of previous diagonal * cos - sin * previous diagonal
        dPnm.append(dPnm[-1]*costh*c2 - sinth * Pnm_tmp)

        for n in range(m+2, nmax+1):
            d = n * n - m * m
            e = n + n - 1


            Pnm.append((e * costh *  Pnm[-1] - rootn[d-e] * Pnm[-2])
                         / rootn[d])
            #Must be after Pnm is set to be setting derivative with same n
            #dPnm[n,m] = (n*cos*Pnm[n,m] - d^1/2 * Pnm[n-1,m])/sin**2

            dPnm.append((n*costh*Pnm[-1] - rootn[d]* Pnm[-2])/sinth)

        if m > 0:#Diagonal append Pnm[m,m] = sin^m(theta)
            Pnm.append(sinth*Pnm_tmp / rootn[m+m+2])
            dPnm.append((dPnm_diag_tmp* rootn[m+1] * np.sqrt(.5)))#/(2*(m+1)+1))

    #Pnm.pop(0)
    #dPnm.pop(0)
    return[Pnm, dPnm]