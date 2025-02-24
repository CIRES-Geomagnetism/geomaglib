import math
from typing import Union

import numpy as np
import warnings

def Flattened_Chaos_Legendre1(nmax: int, theta: Union[float, list[float]], epsilon: float = 1e-6, silence_warnings: bool = False) -> list[list[float]]:
    """
    Inputs:
    nmax (int): maximum degree of legendre polynomial (all orders for each degree are calculated)
    theta (float or np.array(float)): colatitude values where the legendre polynomials are calculated
    epsilon (float): By how much pole values are shifted away from poles to avoid divide by zero (Should be > 1e-6)
    silence_warnings (bool): True -> do not display warning that theta values < epsilon have been shifted. False -> display

    Outputs:
    Takes nmax and colatitude (degrees) as inputs
    Outputs a 2 dimensional numpy array which contains the associated
    legendre polynomials (Pnm) and the respective derivatives (dPnm).
    In these respective arrays the values are arranged in order major order:
    i.e.
    P10,P20,P30...Pnmax0, P11,P21,P31...Pnmax1
    """

    if np.isscalar(theta):
        theta = np.array([theta])
    if (np.isclose(0, min(theta), epsilon)) or (np.isclose(max(theta), 180, epsilon)):
        if not silence_warnings:
            warnings.warn(f'Input coordinates include the poles. They have been shifted by {epsilon}')
        mask_0 = np.isclose(theta, 0, atol=epsilon)
        mask_180 = np.isclose(theta, 180, atol=epsilon)
        theta[mask_0 | mask_180] = epsilon

    if min(theta) < 0.0 or max(theta) > 180.0:
        raise ValueError('Colatitude outside bounds [0, 180].')

    costh = np.cos(np.radians(theta))
    sinth = np.sqrt(1 - costh * costh)

    Pnm = []
    dPnm = []
    Pnm.append(np.ones(len(costh)))  # is copied into trailing dimensions
    dPnm.append(np.zeros(len(costh)))

    rootn = np.sqrt(np.arange(2 * nmax ** 2 + 1))

    # Recursion relations after Langel "The Main Field" (1987),
    # eq. (27) and Table 2 (p. 256)
    for m in range(nmax):
        if m == 1:
            Pnm.append(sinth)
            dPnm.append(costh)
        # Buffer normalization factor for [n,n] to [n,n-1]
        c2 = rootn[m + m + 1]
        # Pnm_tmp holds Pnm[m,m] with normalization factor
        Pnm_tmp = c2 * Pnm[-1]
        # Pnm[m,m-1] = Pnm[m-1,m-1] * cos * c2
        Pnm.append(costh * Pnm_tmp)
        # Hold this for the derivative of diagonal later
        dPnm_diag_tmp = Pnm[-1]
        # Derivative of previous diagonal * cos - sin * previous diagonal
        dPnm.append(dPnm[-1] * costh * c2 - sinth * Pnm_tmp)

        for n in range(m + 2, nmax + 1):
            d = n * n - m * m
            e = n + n - 1

            Pnm.append((e * costh * Pnm[-1] - rootn[d - e] * Pnm[-2]) / rootn[d])
            dPnm.append((n * costh * Pnm[-1] - rootn[d] * Pnm[-2]) / sinth)

        if m > 0:  # Diagonal append Pnm[m,m] = sin^m(theta)
            Pnm.append(sinth * Pnm_tmp / rootn[m + m + 2])
            dPnm.append((dPnm_diag_tmp * rootn[m + 1] * np.sqrt(0.5)))
    if nmax == 1:
        Pnm.append(sinth)
        dPnm.append(costh)
    return [Pnm, dPnm]

def get_index(n: int, m:int, nmax: int) -> int:
    """
    Get the index of Flattened_Chaos_Legendre1
    Args:
        n: degree number
        m: order number
        nmax: maximum degree or order

    Returns:
        The index of legendre polynomial value in Flattened_Chaos_Legendre1

    The order of legendre polynomial in Flattened_Chaos_Legendre1()
        {degree}{order}
        00 (degree=0, order=0)
        10
        20
        30
        ...
        {nmax}0
        11
        21
        ...
        {nmax}1
        22
        32
        ...

    """
    if m == 0:
        return n
    else:
        return int(m * ((nmax + 1) + (nmax + 1 - (m - 1))) // 2 + (n - m))

