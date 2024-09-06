import numpy as np
import math
from collections import defaultdict
import warnings


#Leg_SHA_optimized_ind
def Leg_SHA_optimized_ind(rho, theta, phi, nmax, mmax, get='all', epsilon = 1e-2, silence_warnings = False):
    """
    Written by Collin Kadlecek 7/2024 Collin.kadlecek@noaa.gov
    Computes matrices to connect the radial, colatitude and azimuthal field
    components to the magnetic potential field in terms of spherical harmonic
    coefficients (Schmidt quasi-normalized). In degree order 
    i.e. g10,g11, h11, g20, g21, h21, g22, h22, ... ,g_nmax_nmax h_nmax_nmax

    Parameters
    ----------

    radius : ndarray, shape (...)
        Array containing the radius in earth radii.
    theta : ndarray, shape (...)
        Array containing the colatitude in degrees
        :math:`[0^\\circ,180^\\circ]`.
    phi : ndarray, shape (...)
        Array containing the longitude in degrees.
    nmax : int, positive
        Maximum degree of the sphercial harmonic expansion.
    mmax : int, positive
        Maximum order of the spherical harmonic expansion (defaults to
        ``nmax``). For ``mmax = 0``, for example, only the zonal terms are
        returned.
    get : {'r', 'all'}, optional
        Return just the radial component of the field, or all components

    Returns
    -------
    A_radius, A_theta, A_phi : ndarray, shape (3, nmax*(nmax+2), len(theta))
        Matrices for radial, colatitude and azimuthal field components. The
        second dimension ``M`` varies depending on ``nmax``, ``nmin`` and
        ``mmax``.

    Algorithm description
    --------
    The implimentation of this function follows the structure of the associated
    legendre polynomial generation function (legendre_poly) created by Nils Olsen from the 
    Chaos magpy package. This function recursively generates the associated legendre
    polynomials. Leg_SHA follows this recursive generations structure. However, it 
    also generates each polynomials derivative for each n,m index. Further, only
    the necessary legendre polynomials are held in memory. Upon each polynomials
    generation, it's associated matrix coefficients are calculated. 

    There are two recursion relations used for the polynomial generation. These
    are equations 26 and 27 from The Main GeoMagnetic Field Langel Chapter 4. 26 is 
    used to generate diagonal (n == m) indeces, and 27 for all other. 

    A different source is used to generate the derivatives. There are two recursion 
    relations which follow the structures of the polynomial generation:
    NOTE:
    dP[l,k] refers to the derivative of legendre polynomial with degree = l and order = k
    (where degree >= order)
    costh = cos(theta) sinth = sin(theta)

    For off diagonal terms (m = n-1) which is a special case of 27:
    dP[n,n-1] = dP[n-1,n-1]*costh*sqrt(n+n+1)-sinth*P[n-1,n-1]*sqrt(n+n+1)
    
    For diagonal terms (m = n) following 26:
    dP[n,n] = P[n,n-1]*sqrt(n)/sqrt(2)

    For all other terms:
    dP[n,m] = (n*costh*P[n,m] - sqrt(n**2-m**2)*P[n-1,m])/sinth
    

    Warnings
    --------
    The function can also return the design matrix for the field components at
    the geographic poles, i.e. where ``theta == 0.`` or ``theta == 180.``.
    However, users should be careful when doing this because the vector basis
    for spherical geocentric coordinates,
    :math:`{{\\mathbf{e}_r, \\mathbf{e}_\\theta, \\mathbf{e}_\\phi}}`,
    depends on longitude, which is not well defined at the poles. That is,
    at the poles, any value for the longitude maps to the same location in
    euclidean coordinates but gives a different vector basis in spherical
    geocentric coordinates. Nonetheless, by choosing a specific value for the
    longitude at the poles, users can define the vector basis, which then
    establishes the meaning of the spherical geocentric components. The vector
    basis for the horizontal components is defined as

    .. math::

        \\mathbf{e}_\\theta &= \\cos\\theta\\cos\\phi\\mathbf{e}_x -
            \\cos\\theta\\sin\\phi\\mathbf{e}_y - \\sin\\theta\\mathbf{e}_z\\\\
        \\mathbf{e}_\\phi &= -\\sin\\phi\\mathbf{e}_x +
            \\cos\\phi\\mathbf{e}_y

    Hence, at the geographic north pole as given by ``theta = 0.`` and
    ``phi = 0.`` (chosen by the user), the returned design matrix ``A_theta``
    will be for components along the direction
    :math:`\\mathbf{e}_\\theta = \\mathbf{e}_x` and ``A_phi`` for components
    along :math:`\\mathbf{e}_\\phi = \\mathbf{e}_y`. However,
    if ``phi = 180.`` is chosen, ``A_theta`` will be for components along
    :math:`\\mathbf{e}_\\theta = -\\mathbf{e}_x` and ``A_phi``
    along :math:`\\mathbf{e}_\\phi = -\\mathbf{e}_y`.

    Examples
    --------
    Create the design matrices given 4 locations on the Earth's surface:

    >>> r = 1 = 6371.2km
    >>> theta = np.array([90., 100., 110., 120.])
    >>> phi = 0.
    >>> A_radius, A_theta, A_phi = design_gauss(r, theta, phi, nmax=1)
    >>> A_radius
    array([[ 1.22464680e-16, 2.00000000e+00, 0.00000000e+00],
           [-3.47296355e-01, 1.96961551e+00, 0.00000000e+00],
           [-6.84040287e-01, 1.87938524e+00, 0.00000000e+00],
           [-1.00000000e+00, 1.73205081e+00, 0.00000000e+00]])

    Say, only the axial dipole coefficient is non-zero, what is the magnetic
    field in terms of spherical geocentric components?

    >>> m = np.array([-30000, 0., 0.])  # model coefficients in nT
    >>> Br = A_radius @ m; Br
    array([-3.67394040e-12, 1.04188907e+04, 2.05212086e+04, 3.00000000e+04])
    >>> Bt = A_theta @ m; Bt
    array([-30000. , -29544.23259037, -28190.77862358, -25980.76211353])
    >>> Bp = A_phi @ m; Bp
    array([0., 0., 0., 0.])

    A more complete example is given below:

    .. code-block:: python

        import numpy as np
        from chaosmagpy.model_utils import design_gauss, synth_values

        nmax = 5  # desired truncation degree, i.e. 35 model coefficients
        coeffs = np.arange(35)  # example model coefficients

        # example locations
        N = 10
        radius = 1 = 6371.2 km # Earth's surface
        phi = np.linspace(-180., 180., num=N)  # azimuth in degrees
        theta = np.linspace(1., 179., num=N)  # colatitude in degrees

        # compute design matrices to compute predictions of the
        # field components from the model coefficients, each is of shape N x 35
        A_radius, A_theta, A_phi = design_gauss(r, theta, phi, nmax=nmax)

        # magnetic components computed from the model
        Br_pred = A_radius @ coeffs
        Bt_pred = A_theta @ coeffs
        Bp_pred = A_phi @ coeffs

        # compute the magnetic components directly from the coefficients
        Br, Bt, Bp = synth_values(coeffs, radius, theta, phi)
        np.linalg.norm(Br - Br_pred)  # approx 0.
        np.linalg.norm(Bt - Bt_pred)  # approx 0.
        np.linalg.norm(Bp - Bp_pred)  # approx 0.

    """
    del_indexes = []
    def calc_sha(k,n,m,rn1,Pnm,dPnm, cos_phi_mem,sin_phi_mem, sin_theta):
        k = k - 1
        dPnm = -dPnm
        k = 2*k
        nr = (n+1) * rn1
        tr = rn1*inverse_sinth
        A_r[k] = nr*Pnm*cos_phi_mem
        A_theta[k] = rn1*dPnm*cos_phi_mem
        A_phi[k] = -tr*Pnm*(-m*sin_phi_mem)
        if(m!=0):
            k+=1
            #h terms as in h11, h21, h22 etc...
            A_r[k] = nr*Pnm*sin_phi_mem
            A_theta[k] = rn1*dPnm*sin_phi_mem
            A_phi[k] = -tr*Pnm*m*cos_phi_mem
        else:
            del_indexes.append(k+1)
            A_r[k +1] = np.nan + np.zeros(len(costh))
            A_theta[k +1] = np.nan + np.zeros(len(costh))
            A_phi[k +1] = np.nan + np.zeros(len(costh))

    if get=='r':
        calcTP=False
    else:
        calcTP=True
    N_data=np.size(theta,0)

    theta_min = np.amin(theta)
    theta_max = np.amax(theta)

    # check if poles are included
    if (np.isclose(0,theta_min,epsilon)) or (np.isclose(theta_max,180,epsilon)):
        if(not silence_warnings):
            warnings.warn(f'Input coordinates include the poles. They have been shifted by {epsilon}')
        mask_0 = np.isclose(theta, 0, atol=epsilon)
        mask_180 = np.isclose(theta, 180, atol=epsilon)
        theta[mask_0 | mask_180] = epsilon
        poles = True
    if (theta_min <= 0.0) or (theta_max >= 180.0):
        raise ValueError('Colatitude outside bounds [0, 180].')
    if (np.isclose(theta_max, 3.14, atol= .01)):
        warnings.warn('Please ensure that your inputs are in earth radii (6371.2km) and degrees for both longitude and colatitude')

    if np.isscalar(theta) is not True and theta.ndim > 1:
        dim2=np.size(theta,1)
    else:
        dim2 = 1
    if dim2 > 1:
        theta=theta.transpose()
        phi=phi.transpose()

    costh = np.cos(np.deg2rad(theta))
    sinth = np.sqrt(1-costh*costh)


    inverse_sinth = 1/sinth

    if np.isscalar(rho) is True:
        rho=rho*np.ones(N_data)
    N_nm=mmax*(mmax + 2) + mmax
    #This looks strange, but np.zeros takes a tuple as input, thus the extra paren is necessary
    A_r=np.zeros((N_nm,N_data))
    if(calcTP):
        A_theta=np.zeros((N_nm,N_data))
        A_phi=np.zeros((N_nm,N_data))

    rootn = np.sqrt(np.arange(2 * nmax**2 + 1))
    # Compute the power array
    powers = -(np.arange(nmax) + 3)

    # Raise each element in rho to each power
    r2n = np.transpose(rho[:, np.newaxis] ** powers)
    # Recursion relations after Langel "The Main Field" (1987),
    # eq. (27) and Table 2 (p. 256)
    phi_rad = np.deg2rad(phi)
    Pnm = []
    dPnm = []


    Pnm.append(np.ones(len(costh)))  # is copied into trailing dimensions
    dPnm.append(np.zeros(len(costh)))
    Pnm.append(np.ones(len(costh)))  # is copied into trailing dimensions
    dPnm.append(np.zeros(len(costh)))
    k=0 #Cycles through all the coefficients for degree n, order m model.


    for m in range(nmax):

        cos_phi = np.cos(phi_rad*m)
        sin_phi = np.sin(phi_rad*m)
        if(m == 1):
            Pnm.append(sinth)
            dPnm.append(costh)
            k = 2
            calc_sha(k,1,1,r2n[0],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)


        # Buffer normalization factor for [n-1,n-1] to [n,n-1]
        c2 = rootn[m + m + 1]
        #Pnm_tmp holds Pnm[m,m] with normalization factor
        Pnm_tmp = c2 * Pnm[-1]
        Pnm.append(costh * Pnm_tmp)

        #Hold this for the derivative of diagonal later
        dPnm_diag_tmp = Pnm[-1]
        #Derivative of previous diagonal * cos - sin * previous diagonal
        dPnm.append(dPnm[-1]*costh*c2 - sinth * Pnm_tmp)
        kstart = int((m+1)*(m+2)/2+m)
        calc_sha(kstart,m+1,m,r2n[m],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)

        Pnm.pop(0)
        dPnm.pop(0)
        k = kstart
        for n in range(m+2, nmax+1):
            d = n * n - m * m
            e = n + n - 1
            k = k + n

            Pnm.append((e * costh *  Pnm[-1] - rootn[d-e] * Pnm[-2])
                         / rootn[d])
            #Must be after Pnm is set to be setting derivative with same n
            #dPnm[n,m] = (n*cos*Pnm[n,m] - d^1/2 * Pnm[n-1,m])/sin**2

            dPnm.append((n*costh*Pnm[-1] - rootn[d]* Pnm[-2])*inverse_sinth)
            calc_sha(k,n,m,r2n[n-1],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)
            Pnm.pop(0)
            dPnm.pop(0)

        if m > 0:#Diagonal append Pnm[m,m] = sin^m(theta)
            cos_phi = np.cos(phi_rad*(m+1))
            sin_phi = np.sin(phi_rad*(m+1))
            Pnm.append(sinth*Pnm_tmp / rootn[m+m+2])
            dPnm.append((dPnm_diag_tmp* rootn[m+1] * np.sqrt(.5)))#/(2*(m+1)+1))
            k = kstart + 1
            calc_sha(k,m+1,m+1,r2n[m],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)
            Pnm.pop(0)
            dPnm.pop(0)


    mask = ~np.isnan(A_r).any(axis=1)
    A_r = A_r[mask]
    A_theta = A_theta[mask]
    A_phi = A_phi[mask]

    if calcTP:
        return A_r, A_theta, A_phi

    else:
        return A_r


    # ### Leg_SHA_optimized
def Leg_SHA_optimized(rho, theta, phi, nmax, mmax, get='all',epsilon = 1e-2, silence_warnings = False):
    """
    Written by Collin Kadlecek 7/2024, Collin.kadlecek@noaa.gov
    Computes matrices to connect the radial, colatitude and azimuthal field
    components to the magnetic potential field in terms of spherical harmonic
    coefficients (Schmidt quasi-normalized). In order major order
    i.e. g10,g20,g30...gnmax0, g11,h11,g21,h21...

    Parameters
    ----------

    radius : ndarray, shape (...)
        Array containing the radius in earth radii.
    theta : ndarray, shape (...)
        Array containing the colatitude in degrees
        :math:`[0^\\circ,180^\\circ]`.
    phi : ndarray, shape (...)
        Array containing the longitude in degrees.
    nmax : int, positive
        Maximum degree of the sphercial harmonic expansion.
    mmax : int, positive
        Maximum order of the spherical harmonic expansion (defaults to
        ``nmax``). For ``mmax = 0``, for example, only the zonal terms are
        returned.
    get : {'r', 'all'}, optional
        Return just the radial component of the field, or all components

    Returns
    -------
    A_radius, A_theta, A_phi : ndarray, shape (3, nmax*(nmax+2), len(theta))
        Matrices for radial, colatitude and azimuthal field components. The
        second dimension ``M`` varies depending on ``nmax``, ``nmin`` and
        ``mmax``.

    Algorithm description
    --------
    The implimentation of this function follows the structure of the associated
    legendre polynomial generation function (legendre_poly) created by Nils Olsen from the
    Chaos magpy package. This function recursively generates the associated legendre
    polynomials. Leg_SHA follows this recursive generations structure. However, it
    also generates each polynomials derivative for each n,m index. Further, only
    the necessary legendre polynomials are held in memory. Upon each polynomials
    generation, it's associated matrix coefficients are calculated.

    There are two recursion relations used for the polynomial generation. These
    are equations 26 and 27 from The Main GeoMagnetic Field Langel Chapter 4. 26 is
    used to generate diagonal (n == m) indeces, and 27 for all other.

    The derivative is generated with a recurrence relation from relation from Holmes and
    Featherstone. There are two recursion relations which follow the structures of the
    polynomial generation:
    NOTE:
    dP[l,k] refers to the derivative of legendre polynomial with degree = l and order = k
    (where degree >= order)
    costh = cos(theta) sinth = sin(theta)

    For off diagonal terms (m = n-1) which is a special case of 27:
    dP[n,n-1] = dP[n-1,n-1]*costh*sqrt(n+n+1)-sinth*P[n-1,n-1]*sqrt(n+n+1)

    For diagonal terms (m = n) following 26:
    dP[n,n] = P[n,n-1]*sqrt(n)/sqrt(2)

    For all other terms:
    dP[n,m] = (n*costh*P[n,m] - sqrt(n**2-m**2)*P[n-1,m])/sinth


    Warnings
    --------
    The function can also return the design matrix for the field components at
    the geographic poles, i.e. where ``theta == 0.`` or ``theta == 180.``.
    However, users should be careful when doing this because the vector basis
    for spherical geocentric coordinates,
    :math:`{{\\mathbf{e}_r, \\mathbf{e}_\\theta, \\mathbf{e}_\\phi}}`,
    depends on longitude, which is not well defined at the poles. That is,
    at the poles, any value for the longitude maps to the same location in
    euclidean coordinates but gives a different vector basis in spherical
    geocentric coordinates. Nonetheless, by choosing a specific value for the
    longitude at the poles, users can define the vector basis, which then
    establishes the meaning of the spherical geocentric components. The vector
    basis for the horizontal components is defined as

    .. math::

        \\mathbf{e}_\\theta &= \\cos\\theta\\cos\\phi\\mathbf{e}_x -
            \\cos\\theta\\sin\\phi\\mathbf{e}_y - \\sin\\theta\\mathbf{e}_z\\\\
        \\mathbf{e}_\\phi &= -\\sin\\phi\\mathbf{e}_x +
            \\cos\\phi\\mathbf{e}_y

    Hence, at the geographic north pole as given by ``theta = 0.`` and
    ``phi = 0.`` (chosen by the user), the returned design matrix ``A_theta``
    will be for components along the direction
    :math:`\\mathbf{e}_\\theta = \\mathbf{e}_x` and ``A_phi`` for components
    along :math:`\\mathbf{e}_\\phi = \\mathbf{e}_y`. However,
    if ``phi = 180.`` is chosen, ``A_theta`` will be for components along
    :math:`\\mathbf{e}_\\theta = -\\mathbf{e}_x` and ``A_phi``
    along :math:`\\mathbf{e}_\\phi = -\\mathbf{e}_y`.

    Examples
    --------
    Create the design matrices given 4 locations on the Earth's surface:

    >>> r = 1 = 6371.2km
    >>> theta = np.array([90., 100., 110., 120.])
    >>> phi = 0.
    >>> A_radius, A_theta, A_phi = design_gauss(r, theta, phi, nmax=1)
    >>> A_radius
    array([[ 1.22464680e-16, 2.00000000e+00, 0.00000000e+00],
           [-3.47296355e-01, 1.96961551e+00, 0.00000000e+00],
           [-6.84040287e-01, 1.87938524e+00, 0.00000000e+00],
           [-1.00000000e+00, 1.73205081e+00, 0.00000000e+00]])

    Say, only the axial dipole coefficient is non-zero, what is the magnetic
    field in terms of spherical geocentric components?

    >>> m = np.array([-30000, 0., 0.])  # model coefficients in nT
    >>> Br = A_radius @ m; Br
    array([-3.67394040e-12, 1.04188907e+04, 2.05212086e+04, 3.00000000e+04])
    >>> Bt = A_theta @ m; Bt
    array([-30000. , -29544.23259037, -28190.77862358, -25980.76211353])
    >>> Bp = A_phi @ m; Bp
    array([0., 0., 0., 0.])

    A more complete example is given below:

    .. code-block:: python

        import numpy as np
        from chaosmagpy.model_utils import design_gauss, synth_values

        nmax = 5  # desired truncation degree, i.e. 35 model coefficients
        coeffs = np.arange(35)  # example model coefficients

        # example locations
        N = 10
        radius = 1 = 6371.2 km # Earth's surface
        phi = np.linspace(-180., 180., num=N)  # azimuth in degrees
        theta = np.linspace(1., 179., num=N)  # colatitude in degrees

        # compute design matrices to compute predictions of the
        # field components from the model coefficients, each is of shape N x 35
        A_radius, A_theta, A_phi = design_gauss(r, theta, phi, nmax=nmax)

        # magnetic components computed from the model
        Br_pred = A_radius @ coeffs
        Bt_pred = A_theta @ coeffs
        Bp_pred = A_phi @ coeffs

        # compute the magnetic components directly from the coefficients
        Br, Bt, Bp = synth_values(coeffs, radius, theta, phi)
        np.linalg.norm(Br - Br_pred)  # approx 0.
        np.linalg.norm(Bt - Bt_pred)  # approx 0.
        np.linalg.norm(Bp - Bp_pred)  # approx 0.

    """

    #Locally define a function for calculating SHA
    def calc_sha(k,n,m,rn1,Pnm,dPnm, cos_phi_mem,sin_phi_mem, sin_theta):
        dPnm = -dPnm
        nr = (n+1) * rn1
        tr = rn1*inverse_sinth
        #g terms as in g11, g21, g22 etc...
        A_r[k] = nr*Pnm*cos_phi_mem
        A_theta[k] = rn1*dPnm*cos_phi_mem
        A_phi[k] = tr*Pnm*(m*sin_phi_mem)
        k+=1

        if(m!=0):

            #h terms as in h11, h21, h22 etc...
            A_r[k] = nr*Pnm*sin_phi_mem
            A_theta[k] = rn1*dPnm*sin_phi_mem
            A_phi[k] = -tr*Pnm*m*cos_phi_mem
            k = k + 1
        # else:
        #     A_r[k] = np.nan + np.zeros(len(costh))
        #     A_theta[k] = np.nan + np.zeros(len(costh))
        #     A_phi[k] = np.nan + np.zeros(len(costh))
        return k


    theta_min = np.amin(theta)
    theta_max = np.amax(theta)

    # check if poles are included
    if (np.isclose(0,theta_min,atol = epsilon)) or (np.isclose(theta_max,180,atol = epsilon)):
        if(not silence_warnings):
            strtemp = f'Input coordinates include the poles. They have been shifted by {epsilon}'
            warnings.warn(strtemp, stacklevel=2)
        mask_0 = np.isclose(theta, 0, atol=epsilon)
        mask_180 = np.isclose(theta, 180, atol=epsilon)
        theta[mask_0 | mask_180] = epsilon
        poles = True
    if (theta_min < 0.0) or (theta_max > 180.0):
        raise ValueError('Colatitude outside bounds [0, 180].')

    if(not silence_warnings):
        if (np.isclose(theta_max, 3.14, atol= .01)):
            warnings.warn('Please ensure that your inputs are in earth radii (6371.2km)\n and degrees for both longitude and colatitude')

    if get=='r':
        calcTP=False
    else:
        calcTP=True
    N_data=np.size(theta,0)
    if (np.size(phi,0) != N_data or np.size(rho,0) != N_data):
        raise ValueError(r"Leg_SHA_optimized: rho, theta, and phi have different dimensions")
    if np.isscalar(theta) is not True and theta.ndim > 1:
        dim2=np.size(theta,1)
    else:
        dim2 = 1
    if dim2 > 1:
        theta=theta.transpose()
        phi=phi.transpose()

    costh = np.cos(np.deg2rad(theta))
    sinth = np.sqrt(1-costh*costh)
    inverse_sinth = 1/sinth
    if np.isscalar(rho) is True:
        rho=rho*np.ones(N_data)
    # N_nm=mmax*(mmax + 2) + (nmax + mmax)*(2*mmax + 1)
    N_nm=mmax*(mmax + 2)# (nmax - mmax)*(2*mmax + 1)
    # calculate sub-matrices for s=0, p=0
    #This looks strange, but np.zeros takes a tuple as input, thus the extra paren is necessary
    A_r=np.zeros((N_nm,N_data))
    if(calcTP):
        A_theta=np.zeros((N_nm,N_data))
        A_phi=np.zeros((N_nm,N_data))

    rootn = np.sqrt(np.arange(2 * nmax**2 + 1))
    # Compute the power array
    powers = -(np.arange(nmax) + 3)

    # Raise each element in rho to each power
    r2n = np.transpose(rho[:, np.newaxis] ** powers)
    phi_rad = np.deg2rad(phi)
    Pnm = []
    dPnm = []

    #This is done twice to begin the Pnm queue
    Pnm.append(np.ones(len(costh)))  # is copied into trailing dimensions
    dPnm.append(np.zeros(len(costh)))
    Pnm.append(np.ones(len(costh)))
    dPnm.append(np.zeros(len(costh)))
    k=0 #Cycles through all the coefficients for degree n, order m model.

    # k = calc_sha(k,0,0,r2n[0],1,0,np.cos(phi_rad),np.sin(phi_rad),sinth)



    cos_phi = np.ones(len(phi))
    sin_phi = cos_phi-1
    for m in range(nmax):


        if(m == 1):

            Pnm.append(sinth)
            dPnm.append(costh)
            k = calc_sha(k,1,1,r2n[0],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)


        # Buffer normalization factor for [n-1,n-1] to [n,n-1]
        c2 = rootn[m + m + 1]
        #Pnm_tmp holds Pnm[m,m] with normalization factor
        Pnm_tmp = c2 * Pnm[-1]
        #Pnm[m,m-1] = Pnm[m-1,m-1]*cos*c2
        Pnm.append(costh * Pnm_tmp)

        #Hold this for the derivative of diagonal later
        dPnm_diag_tmp = Pnm[-1]
        #Derivative of previous diagonal * cos - sin * previous diagonal
        dPnm.append(dPnm[-1]*costh*c2 - sinth * Pnm_tmp)
        k = calc_sha(k,m+1,m,r2n[m],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)
        Pnm.pop(0)
        dPnm.pop(0)

        for n in range(m+2, nmax+1):
            d = n * n - m * m
            e = n + n - 1


            Pnm.append((e * costh *  Pnm[-1] - rootn[d-e] * Pnm[-2])
                         / rootn[d])
            #dPnm[n,m] = (n*cos*Pnm[n,m] - d^1/2 * Pnm[n-1,m])/sin**2

            dPnm.append((n*costh*Pnm[-1] - rootn[d]* Pnm[-2])*inverse_sinth)
            k = calc_sha(k,n,m,r2n[n-1],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)
            Pnm.pop(0)
            dPnm.pop(0)

        cos_phi = np.cos(phi_rad*(m+1))
        sin_phi = np.sin(phi_rad*(m+1))
        if m > 0:#Diagonal append Pnm[m,m] = sin^m(theta)

            Pnm.append(sinth*Pnm_tmp / rootn[m+m+2])
            dPnm.append((dPnm_diag_tmp* rootn[m+1] * np.sqrt(.5)))

            k = calc_sha(k,m+1,m+1,r2n[m],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)

            Pnm.pop(0)
            dPnm.pop(0)

    # dP(n,m) = Pnm(m,n+1) is the derivative of P(n,m) vrt. theta
    if calcTP:
        return A_r, A_theta, A_phi
    else:
        return A_r

def Leg_SHA_optimized_with_nmin(rho, theta, phi, nmin, nmax, get='all',epsilon = 1e-2, silence_warnings = False):
    """
    Written by Collin Kadlecek 7/2024, Collin.kadlecek@noaa.gov
    Computes matrices to connect the radial, colatitude and azimuthal field
    components to the magnetic potential field in terms of spherical harmonic
    coefficients (Schmidt quasi-normalized). In order major order
    i.e. g10,g20,g30...gnmax0, g11,h11,g21,h21...

    Parameters
    ----------

    radius : ndarray, shape (...)
        Array containing the radius in earth radii.
    theta : ndarray, shape (...)
        Array containing the colatitude in degrees
        :math:`[0^\\circ,180^\\circ]`.
    phi : ndarray, shape (...)
        Array containing the longitude in degrees.
    nmax : int, positive
        Maximum degree of the sphercial harmonic expansion.
    mmax : int, positive
        Maximum order of the spherical harmonic expansion (defaults to
        ``nmax``). For ``mmax = 0``, for example, only the zonal terms are
        returned.
    get : {'r', 'all'}, optional
        Return just the radial component of the field, or all components

    Returns
    -------
    A_radius, A_theta, A_phi : ndarray, shape (3, nmax*(nmax+2), len(theta))
        Matrices for radial, colatitude and azimuthal field components. The
        second dimension ``M`` varies depending on ``nmax``, ``nmin`` and
        ``mmax``.

    Algorithm description
    --------
    The implimentation of this function follows the structure of the associated
    legendre polynomial generation function (legendre_poly) created by Nils Olsen from the
    Chaos magpy package. This function recursively generates the associated legendre
    polynomials. Leg_SHA follows this recursive generations structure. However, it
    also generates each polynomials derivative for each n,m index. Further, only
    the necessary legendre polynomials are held in memory. Upon each polynomials
    generation, it's associated matrix coefficients are calculated.

    There are two recursion relations used for the polynomial generation. These
    are equations 26 and 27 from The Main GeoMagnetic Field Langel Chapter 4. 26 is
    used to generate diagonal (n == m) indeces, and 27 for all other.

    The derivative is generated with a recurrence relation from relation from Holmes and
    Featherstone. There are two recursion relations which follow the structures of the
    polynomial generation:
    NOTE:
    dP[l,k] refers to the derivative of legendre polynomial with degree = l and order = k
    (where degree >= order)
    costh = cos(theta) sinth = sin(theta)

    For off diagonal terms (m = n-1) which is a special case of 27:
    dP[n,n-1] = dP[n-1,n-1]*costh*sqrt(n+n+1)-sinth*P[n-1,n-1]*sqrt(n+n+1)

    For diagonal terms (m = n) following 26:
    dP[n,n] = P[n,n-1]*sqrt(n)/sqrt(2)

    For all other terms:
    dP[n,m] = (n*costh*P[n,m] - sqrt(n**2-m**2)*P[n-1,m])/sinth


    Warnings
    --------
    The function can also return the design matrix for the field components at
    the geographic poles, i.e. where ``theta == 0.`` or ``theta == 180.``.
    However, users should be careful when doing this because the vector basis
    for spherical geocentric coordinates,
    :math:`{{\\mathbf{e}_r, \\mathbf{e}_\\theta, \\mathbf{e}_\\phi}}`,
    depends on longitude, which is not well defined at the poles. That is,
    at the poles, any value for the longitude maps to the same location in
    euclidean coordinates but gives a different vector basis in spherical
    geocentric coordinates. Nonetheless, by choosing a specific value for the
    longitude at the poles, users can define the vector basis, which then
    establishes the meaning of the spherical geocentric components. The vector
    basis for the horizontal components is defined as

    .. math::

        \\mathbf{e}_\\theta &= \\cos\\theta\\cos\\phi\\mathbf{e}_x -
            \\cos\\theta\\sin\\phi\\mathbf{e}_y - \\sin\\theta\\mathbf{e}_z\\\\
        \\mathbf{e}_\\phi &= -\\sin\\phi\\mathbf{e}_x +
            \\cos\\phi\\mathbf{e}_y

    Hence, at the geographic north pole as given by ``theta = 0.`` and
    ``phi = 0.`` (chosen by the user), the returned design matrix ``A_theta``
    will be for components along the direction
    :math:`\\mathbf{e}_\\theta = \\mathbf{e}_x` and ``A_phi`` for components
    along :math:`\\mathbf{e}_\\phi = \\mathbf{e}_y`. However,
    if ``phi = 180.`` is chosen, ``A_theta`` will be for components along
    :math:`\\mathbf{e}_\\theta = -\\mathbf{e}_x` and ``A_phi``
    along :math:`\\mathbf{e}_\\phi = -\\mathbf{e}_y`.

    Examples
    --------
    Create the design matrices given 4 locations on the Earth's surface:

    >>> r = 1 = 6371.2km
    >>> theta = np.array([90., 100., 110., 120.])
    >>> phi = 0.
    >>> A_radius, A_theta, A_phi = design_gauss(r, theta, phi, nmax=1)
    >>> A_radius
    array([[ 1.22464680e-16, 2.00000000e+00, 0.00000000e+00],
           [-3.47296355e-01, 1.96961551e+00, 0.00000000e+00],
           [-6.84040287e-01, 1.87938524e+00, 0.00000000e+00],
           [-1.00000000e+00, 1.73205081e+00, 0.00000000e+00]])

    Say, only the axial dipole coefficient is non-zero, what is the magnetic
    field in terms of spherical geocentric components?

    >>> m = np.array([-30000, 0., 0.])  # model coefficients in nT
    >>> Br = A_radius @ m; Br
    array([-3.67394040e-12, 1.04188907e+04, 2.05212086e+04, 3.00000000e+04])
    >>> Bt = A_theta @ m; Bt
    array([-30000. , -29544.23259037, -28190.77862358, -25980.76211353])
    >>> Bp = A_phi @ m; Bp
    array([0., 0., 0., 0.])

    A more complete example is given below:

    .. code-block:: python

        import numpy as np
        from chaosmagpy.model_utils import design_gauss, synth_values

        nmax = 5  # desired truncation degree, i.e. 35 model coefficients
        coeffs = np.arange(35)  # example model coefficients

        # example locations
        N = 10
        radius = 1 = 6371.2 km # Earth's surface
        phi = np.linspace(-180., 180., num=N)  # azimuth in degrees
        theta = np.linspace(1., 179., num=N)  # colatitude in degrees

        # compute design matrices to compute predictions of the
        # field components from the model coefficients, each is of shape N x 35
        A_radius, A_theta, A_phi = design_gauss(r, theta, phi, nmax=nmax)

        # magnetic components computed from the model
        Br_pred = A_radius @ coeffs
        Bt_pred = A_theta @ coeffs
        Bp_pred = A_phi @ coeffs

        # compute the magnetic components directly from the coefficients
        Br, Bt, Bp = synth_values(coeffs, radius, theta, phi)
        np.linalg.norm(Br - Br_pred)  # approx 0.
        np.linalg.norm(Bt - Bt_pred)  # approx 0.
        np.linalg.norm(Bp - Bp_pred)  # approx 0.

    """
    if(nmin == 0):
        nmin = -1
    #Locally define a function for calculating SHA
    def calc_sha(k,n,m,rn1,Pnm,dPnm, cos_phi_mem,sin_phi_mem, sin_theta):
        if(n < nmin):
            return k
        dPnm = -dPnm
        nr = (n+1) * rn1
        tr = rn1*inverse_sinth
        #g terms as in g11, g21, g22 etc...
        A_r[k] = nr*Pnm*cos_phi_mem
        A_theta[k] = rn1*dPnm*cos_phi_mem
        A_phi[k] = tr*Pnm*(m*sin_phi_mem)
        k+=1

        if(m!=0):

            #h terms as in h11, h21, h22 etc...
            A_r[k] = nr*Pnm*sin_phi_mem
            A_theta[k] = rn1*dPnm*sin_phi_mem
            A_phi[k] = -tr*Pnm*m*cos_phi_mem
            k = k + 1
        # else:
        #     A_r[k] = np.nan + np.zeros(len(costh))
        #     A_theta[k] = np.nan + np.zeros(len(costh))
        #     A_phi[k] = np.nan + np.zeros(len(costh))
        return k


    theta_min = np.amin(theta)
    theta_max = np.amax(theta)

    # check if poles are included
    if (np.isclose(0,theta_min,atol = epsilon)) or (np.isclose(theta_max,180,atol = epsilon)):
        if(not silence_warnings):
            strtemp = f'Input coordinates include the poles. They have been shifted by {epsilon}'
            warnings.warn(strtemp, stacklevel=2)
        mask_0 = np.isclose(theta, 0, atol=epsilon)
        mask_180 = np.isclose(theta, 180, atol=epsilon)
        theta[mask_0 | mask_180] = epsilon
        poles = True
    if (theta_min < 0.0) or (theta_max > 180.0):
        raise ValueError('Colatitude outside bounds [0, 180].')

    if(not silence_warnings):
        if (np.isclose(theta_max, 3.14, atol= .01)):
            warnings.warn('Please ensure that your inputs are in earth radii (6371.2km) and degrees for both longitude and colatitude')

    if get=='r':
        calcTP=False
    else:
        calcTP=True
    N_data=np.size(theta,0)
    if (np.size(phi,0) != N_data or np.size(rho,0) != N_data):
        raise ValueError(r"Leg_SHA_optimized: rho, theta, and phi have different dimensions")
    if np.isscalar(theta) is not True and theta.ndim > 1:
        dim2=np.size(theta,1)
    else:
        dim2 = 1
    if dim2 > 1:
        theta=theta.transpose()
        phi=phi.transpose()

    costh = np.cos(np.deg2rad(theta))
    sinth = np.sqrt(1-costh*costh)
    inverse_sinth = 1/sinth
    if np.isscalar(rho) is True:
        rho=rho*np.ones(N_data)
    # N_nm=mmax*(mmax + 2) + (nmax + mmax)*(2*mmax + 1)
    N_nm=nmax*(nmax + 2) - (nmin+1)*(nmin-1)
    # calculate sub-matrices for s=0, p=0
    #This looks strange, but np.zeros takes a tuple as input, thus the extra paren is necessary
    A_r=np.zeros((N_nm,N_data))
    if(calcTP):
        A_theta=np.zeros((N_nm,N_data))
        A_phi=np.zeros((N_nm,N_data))

    rootn = np.sqrt(np.arange(2 * nmax**2 + 1))
    # Compute the power array
    powers = -(np.arange(nmax) + 3)

    # Raise each element in rho to each power
    r2n = np.transpose(rho[:, np.newaxis] ** powers)
    phi_rad = np.deg2rad(phi)
    Pnm = []
    dPnm = []

    #This is done twice to begin the Pnm queue
    Pnm.append(np.ones(len(costh)))  # is copied into trailing dimensions
    dPnm.append(np.zeros(len(costh)))
    Pnm.append(np.ones(len(costh)))
    dPnm.append(np.zeros(len(costh)))
    k=0 #Cycles through all the coefficients for degree n, order m model.

    # k = calc_sha(k,0,0,r2n[0],1,0,np.cos(phi_rad),np.sin(phi_rad),sinth)



    cos_phi = np.ones(len(phi))
    sin_phi = cos_phi-1
    for m in range(nmax):


        if(m == 1):

            Pnm.append(sinth)
            dPnm.append(costh)
            k = calc_sha(k,1,1,r2n[0],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)


        # Buffer normalization factor for [n-1,n-1] to [n,n-1]
        c2 = rootn[m + m + 1]
        #Pnm_tmp holds Pnm[m,m] with normalization factor
        Pnm_tmp = c2 * Pnm[-1]
        #Pnm[m,m-1] = Pnm[m-1,m-1]*cos*c2
        Pnm.append(costh * Pnm_tmp)

        #Hold this for the derivative of diagonal later
        dPnm_diag_tmp = Pnm[-1]
        #Derivative of previous diagonal * cos - sin * previous diagonal
        dPnm.append(dPnm[-1]*costh*c2 - sinth * Pnm_tmp)
        k = calc_sha(k,m+1,m,r2n[m],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)
        Pnm.pop(0)
        dPnm.pop(0)

        for n in range(m+2, nmax+1):
            d = n * n - m * m
            e = n + n - 1


            Pnm.append((e * costh *  Pnm[-1] - rootn[d-e] * Pnm[-2])
                         / rootn[d])
            #dPnm[n,m] = (n*cos*Pnm[n,m] - d^1/2 * Pnm[n-1,m])/sin**2

            dPnm.append((n*costh*Pnm[-1] - rootn[d]* Pnm[-2])*inverse_sinth)
            k = calc_sha(k,n,m,r2n[n-1],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)
            Pnm.pop(0)
            dPnm.pop(0)

        cos_phi = np.cos(phi_rad*(m+1))
        sin_phi = np.sin(phi_rad*(m+1))
        if m > 0:#Diagonal append Pnm[m,m] = sin^m(theta)

            Pnm.append(sinth*Pnm_tmp / rootn[m+m+2])
            dPnm.append((dPnm_diag_tmp* rootn[m+1] * np.sqrt(.5)))

            k = calc_sha(k,m+1,m+1,r2n[m],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)

            Pnm.pop(0)
            dPnm.pop(0)

    # dP(n,m) = Pnm(m,n+1) is the derivative of P(n,m) vrt. theta
    if calcTP:
        return A_r, A_theta, A_phi
    else:
        return A_r


def Leg_SHA_optimized_ind_nmin(rho, theta, phi, nmin, nmax, get='all', epsilon = 1e-2, silence_warnings = False):
    """
    Written by Collin Kadlecek 7/2024 Collin.kadlecek@noaa.gov
    Computes matrices to connect the radial, colatitude and azimuthal field
    components to the magnetic potential field in terms of spherical harmonic
    coefficients (Schmidt quasi-normalized). In degree order 
    i.e. g10,g20,g30...gnmax0, g11,h11,g21,h21...

    Parameters
    ----------

    radius : ndarray, shape (...)
        Array containing the radius in earth radii.
    theta : ndarray, shape (...)
        Array containing the colatitude in degrees
        :math:`[0^\\circ,180^\\circ]`.
    phi : ndarray, shape (...)
        Array containing the longitude in degrees.
    nmax : int, positive
        Maximum degree of the sphercial harmonic expansion.
    mmax : int, positive
        Maximum order of the spherical harmonic expansion (defaults to
        ``nmax``). For ``mmax = 0``, for example, only the zonal terms are
        returned.
    get : {'r', 'all'}, optional
        Return just the radial component of the field, or all components

    Returns
    -------
    A_radius, A_theta, A_phi : ndarray, shape (3, nmax*(nmax+2), len(theta))
        Matrices for radial, colatitude and azimuthal field components. The
        second dimension ``M`` varies depending on ``nmax``, ``nmin`` and
        ``mmax``.

    Algorithm description
    --------
    The implimentation of this function follows the structure of the associated
    legendre polynomial generation function (legendre_poly) created by Nils Olsen from the 
    Chaos magpy package. This function recursively generates the associated legendre
    polynomials. Leg_SHA follows this recursive generations structure. However, it 
    also generates each polynomials derivative for each n,m index. Further, only
    the necessary legendre polynomials are held in memory. Upon each polynomials
    generation, it's associated matrix coefficients are calculated. 

    There are two recursion relations used for the polynomial generation. These
    are equations 26 and 27 from The Main GeoMagnetic Field Langel Chapter 4. 26 is 
    used to generate diagonal (n == m) indeces, and 27 for all other. 

    A different source is used to generate the derivatives. There are two recursion 
    relations which follow the structures of the polynomial generation:
    NOTE:
    dP[l,k] refers to the derivative of legendre polynomial with degree = l and order = k
    (where degree >= order)
    costh = cos(theta) sinth = sin(theta)

    For off diagonal terms (m = n-1) which is a special case of 27:
    dP[n,n-1] = dP[n-1,n-1]*costh*sqrt(n+n+1)-sinth*P[n-1,n-1]*sqrt(n+n+1)
    
    For diagonal terms (m = n) following 26:
    dP[n,n] = P[n,n-1]*sqrt(n)/sqrt(2)

    For all other terms:
    dP[n,m] = (n*costh*P[n,m] - sqrt(n**2-m**2)*P[n-1,m])/sinth
    

    Warnings
    --------
    The function can also return the design matrix for the field components at
    the geographic poles, i.e. where ``theta == 0.`` or ``theta == 180.``.
    However, users should be careful when doing this because the vector basis
    for spherical geocentric coordinates,
    :math:`{{\\mathbf{e}_r, \\mathbf{e}_\\theta, \\mathbf{e}_\\phi}}`,
    depends on longitude, which is not well defined at the poles. That is,
    at the poles, any value for the longitude maps to the same location in
    euclidean coordinates but gives a different vector basis in spherical
    geocentric coordinates. Nonetheless, by choosing a specific value for the
    longitude at the poles, users can define the vector basis, which then
    establishes the meaning of the spherical geocentric components. The vector
    basis for the horizontal components is defined as

    .. math::

        \\mathbf{e}_\\theta &= \\cos\\theta\\cos\\phi\\mathbf{e}_x -
            \\cos\\theta\\sin\\phi\\mathbf{e}_y - \\sin\\theta\\mathbf{e}_z\\\\
        \\mathbf{e}_\\phi &= -\\sin\\phi\\mathbf{e}_x +
            \\cos\\phi\\mathbf{e}_y

    Hence, at the geographic north pole as given by ``theta = 0.`` and
    ``phi = 0.`` (chosen by the user), the returned design matrix ``A_theta``
    will be for components along the direction
    :math:`\\mathbf{e}_\\theta = \\mathbf{e}_x` and ``A_phi`` for components
    along :math:`\\mathbf{e}_\\phi = \\mathbf{e}_y`. However,
    if ``phi = 180.`` is chosen, ``A_theta`` will be for components along
    :math:`\\mathbf{e}_\\theta = -\\mathbf{e}_x` and ``A_phi``
    along :math:`\\mathbf{e}_\\phi = -\\mathbf{e}_y`.

    Examples
    --------
    Create the design matrices given 4 locations on the Earth's surface:

    >>> r = 1 = 6371.2km
    >>> theta = np.array([90., 100., 110., 120.])
    >>> phi = 0.
    >>> A_radius, A_theta, A_phi = design_gauss(r, theta, phi, nmax=1)
    >>> A_radius
    array([[ 1.22464680e-16, 2.00000000e+00, 0.00000000e+00],
           [-3.47296355e-01, 1.96961551e+00, 0.00000000e+00],
           [-6.84040287e-01, 1.87938524e+00, 0.00000000e+00],
           [-1.00000000e+00, 1.73205081e+00, 0.00000000e+00]])

    Say, only the axial dipole coefficient is non-zero, what is the magnetic
    field in terms of spherical geocentric components?

    >>> m = np.array([-30000, 0., 0.])  # model coefficients in nT
    >>> Br = A_radius @ m; Br
    array([-3.67394040e-12, 1.04188907e+04, 2.05212086e+04, 3.00000000e+04])
    >>> Bt = A_theta @ m; Bt
    array([-30000. , -29544.23259037, -28190.77862358, -25980.76211353])
    >>> Bp = A_phi @ m; Bp
    array([0., 0., 0., 0.])

    A more complete example is given below:

    .. code-block:: python

        import numpy as np
        from chaosmagpy.model_utils import design_gauss, synth_values

        nmax = 5  # desired truncation degree, i.e. 35 model coefficients
        coeffs = np.arange(35)  # example model coefficients

        # example locations
        N = 10
        radius = 1 = 6371.2 km # Earth's surface
        phi = np.linspace(-180., 180., num=N)  # azimuth in degrees
        theta = np.linspace(1., 179., num=N)  # colatitude in degrees

        # compute design matrices to compute predictions of the
        # field components from the model coefficients, each is of shape N x 35
        A_radius, A_theta, A_phi = design_gauss(r, theta, phi, nmax=nmax)

        # magnetic components computed from the model
        Br_pred = A_radius @ coeffs
        Bt_pred = A_theta @ coeffs
        Bp_pred = A_phi @ coeffs

        # compute the magnetic components directly from the coefficients
        Br, Bt, Bp = synth_values(coeffs, radius, theta, phi)
        np.linalg.norm(Br - Br_pred)  # approx 0.
        np.linalg.norm(Bt - Bt_pred)  # approx 0.
        np.linalg.norm(Bp - Bp_pred)  # approx 0.

    """
    del_indexes = []
    if(nmin == 0):
        nmin = nmin-1
    def calc_sha(k,n,m,rn1,Pnm,dPnm, cos_phi_mem,sin_phi_mem, sin_theta):
        k = k - 1
        k = 2*k
        if(n < nmin):
            del_indexes.append(k+1)
            A_r[k] = np.nan + np.zeros(len(costh))
            A_theta[k] = np.nan + np.zeros(len(costh))
            A_phi[k] = np.nan + np.zeros(len(costh))

            A_r[k +1] = np.nan + np.zeros(len(costh))
            A_theta[k +1] = np.nan + np.zeros(len(costh))
            A_phi[k +1] = np.nan + np.zeros(len(costh))
            return
        dPnm = -dPnm
        nr = (n+1) * rn1
        tr = rn1*inverse_sinth
        A_r[k] = nr*Pnm*cos_phi_mem
        A_theta[k] = rn1*dPnm*cos_phi_mem
        A_phi[k] = -tr*Pnm*(-m*sin_phi_mem)
        if(m!=0):
            k+=1
            #h terms as in h11, h21, h22 etc...
            A_r[k] = nr*Pnm*sin_phi_mem
            A_theta[k] = rn1*dPnm*sin_phi_mem
            A_phi[k] = -tr*Pnm*m*cos_phi_mem
        else:
            del_indexes.append(k+1)
            A_r[k +1] = np.nan + np.zeros(len(costh))
            A_theta[k +1] = np.nan + np.zeros(len(costh))
            A_phi[k +1] = np.nan + np.zeros(len(costh))

    if get=='r':
        calcTP=False
    else:
        calcTP=True
    N_data=np.size(theta,0)

    theta_min = np.amin(theta)
    theta_max = np.amax(theta)

    # check if poles are included
    if (np.isclose(0,theta_min,epsilon)) or (np.isclose(theta_max,180,epsilon)):
        if(not silence_warnings):
            warnings.warn(f'Input coordinates include the poles. They have been shifted by {epsilon}')
        mask_0 = np.isclose(theta, 0, atol=epsilon)
        mask_180 = np.isclose(theta, 180, atol=epsilon)
        theta[mask_0 | mask_180] = epsilon
        poles = True
    if (theta_min <= 0.0) or (theta_max >= 180.0):
        raise ValueError('Colatitude outside bounds [0, 180].')
    if (np.isclose(theta_max, 3.14, atol= .01)):
        warnings.warn('Please ensure that your inputs are in earth radii (6371.2km) and degrees for both longitude and colatitude')

    if np.isscalar(theta) is not True and theta.ndim > 1:
        dim2=np.size(theta,1)
    else:
        dim2 = 1
    if dim2 > 1:
        theta=theta.transpose()
        phi=phi.transpose()

    costh = np.cos(np.deg2rad(theta))
    sinth = np.sqrt(1-costh*costh)


    inverse_sinth = 1/sinth

    if np.isscalar(rho) is True:
        rho=rho*np.ones(N_data)
    N_nm=nmax*(nmax + 2) + nmax
    #This looks strange, but np.zeros takes a tuple as input, thus the extra paren is necessary
    A_r=np.zeros((N_nm,N_data))
    if(calcTP):
        A_theta=np.zeros((N_nm,N_data))
        A_phi=np.zeros((N_nm,N_data))

    rootn = np.sqrt(np.arange(2 * nmax**2 + 1))
    # Compute the power array
    powers = -(np.arange(nmax) + 3)

    # Raise each element in rho to each power
    r2n = np.transpose(rho[:, np.newaxis] ** powers)
    # Recursion relations after Langel "The Main Field" (1987),
    # eq. (27) and Table 2 (p. 256)
    phi_rad = np.deg2rad(phi)
    Pnm = []
    dPnm = []


    Pnm.append(np.ones(len(costh)))  # is copied into trailing dimensions
    dPnm.append(np.zeros(len(costh)))
    Pnm.append(np.ones(len(costh)))  # is copied into trailing dimensions
    dPnm.append(np.zeros(len(costh)))
    k=0 #Cycles through all the coefficients for degree n, order m model.


    for m in range(nmax):

        cos_phi = np.cos(phi_rad*m)
        sin_phi = np.sin(phi_rad*m)
        if(m == 1):
            Pnm.append(sinth)
            dPnm.append(costh)
            k = 2
            calc_sha(k,1,1,r2n[0],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)


        # Buffer normalization factor for [n-1,n-1] to [n,n-1]
        c2 = rootn[m + m + 1]
        #Pnm_tmp holds Pnm[m,m] with normalization factor
        Pnm_tmp = c2 * Pnm[-1]
        Pnm.append(costh * Pnm_tmp)

        #Hold this for the derivative of diagonal later
        dPnm_diag_tmp = Pnm[-1]
        #Derivative of previous diagonal * cos - sin * previous diagonal
        dPnm.append(dPnm[-1]*costh*c2 - sinth * Pnm_tmp)
        kstart = int((m+1)*(m+2)/2+m)
        calc_sha(kstart,m+1,m,r2n[m],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)

        Pnm.pop(0)
        dPnm.pop(0)
        k = kstart
        for n in range(m+2, nmax+1):
            d = n * n - m * m
            e = n + n - 1
            k = k + n

            Pnm.append((e * costh *  Pnm[-1] - rootn[d-e] * Pnm[-2])
                         / rootn[d])
            #Must be after Pnm is set to be setting derivative with same n
            #dPnm[n,m] = (n*cos*Pnm[n,m] - d^1/2 * Pnm[n-1,m])/sin**2

            dPnm.append((n*costh*Pnm[-1] - rootn[d]* Pnm[-2])*inverse_sinth)
            calc_sha(k,n,m,r2n[n-1],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)
            Pnm.pop(0)
            dPnm.pop(0)

        if m > 0:#Diagonal append Pnm[m,m] = sin^m(theta)
            cos_phi = np.cos(phi_rad*(m+1))
            sin_phi = np.sin(phi_rad*(m+1))
            Pnm.append(sinth*Pnm_tmp / rootn[m+m+2])
            dPnm.append((dPnm_diag_tmp* rootn[m+1] * np.sqrt(.5)))#/(2*(m+1)+1))
            k = kstart + 1
            calc_sha(k,m+1,m+1,r2n[m],Pnm[-1],dPnm[-1],cos_phi,sin_phi,sinth)
            Pnm.pop(0)
            dPnm.pop(0)


    mask = ~np.isnan(A_r).any(axis=1)
    A_r = A_r[mask]
    A_theta = A_theta[mask]
    A_phi = A_phi[mask]

    if calcTP:
        return A_r, A_theta, A_phi

    else:
        return A_r

def permute_SHA(nmin, nmax):
    if(nmin != 0):
        N_nm = nmax*(nmax + 2)- (nmin+1)*(nmin-1)
    if(nmin == 0):
        N_nm = nmax*(nmax + 2)
    # print('total range = ', N_nm)
    permute = np.zeros((N_nm,N_nm))
    i = 0
    for m in range(0,nmax):

        kstart = int((m+1)*(m+2)/2+m)
        if(m+1 > nmin):
            permute[i,kstart] = 1
        i += 1
        k = kstart
        for n in range(m+2,nmax+1):
            k = k + n
            permute[i,k] = 1
            i += 1
        permute[i,kstart + 1] = 1
        i += 1
    return(permute)


def legendre_manoj(lat, nMax):
    """
    Computes  all of the Schmidt-semi normalized associated Legendre
    functions up to degree nMax. If nMax <= 16, function MAG_PcupLow is used.
    Otherwise MAG_PcupHigh is called.
    INPUT  Latitude (degrees)
    nMax        	integer 	 ( Maxumum degree of spherical harmonic secular model)
                    LegendreFunction Pointer to data structure with the following elements
                                                    double *Pcup;  (  pointer to store Legendre Function  )
                                                    double *dPcup; ( pointer to store  Derivative of Legendre function )
    LegendreFunction  Calculated Legendre variables in the data structure
    """

    sin_phi = np.sin(np.radians(lat))

    #    if nMax <= 16 or (1-abs(sin_phi)) < 10**-10:
    [P, dP] = MAG_PcupLow(sin_phi, nMax)
    # hoping that overflow doesn't occur
    #    else:
    #        raise ValueError("Implement Pcup High for this", nMax)

    return [P, dP]


def MAG_PcupLow(x, nMax):
    """ This function evaluates all of the Schmidt-semi normalized associated Legendre
        functions up to degree nMax.
        Calling Parameters:
                INPUT
                        nMax:	 Maximum spherical harmonic degree to compute.
                        x:		cos(colatitude) or sin(latitude).
                OUTPUT
                        Pcup:	A vector of all associated Legendgre polynomials evaluated at
                                        x up to nMax.
                   dPcup: Derivative of Pcup(x) with respect to latitude
        Notes: Overflow may occur if nMax > 20 , especially for high-latitudes.
        Use MAG_PcupHigh for large nMax.
    Written by Manoj Nair, June, 2009 . Manoj.C.Nair@Noaa.Gov.
    Note: In geomagnetism, the derivatives of ALF are usually found with
    respect to the colatitudes. Here the derivatives are found with respect
    to the latitude. The difference is a sign reversal for the derivative of
    the Associated Legendre Functions.
    """

    NumTerms = ((nMax + 1) * (nMax + 2) // 2)
    try:
        Pcup = np.zeros((NumTerms, x.size))
        dPcup = np.zeros((NumTerms, x.size))
    except AttributeError:

        Pcup = np.zeros((NumTerms, 1))
        dPcup = np.zeros((NumTerms, 1))
    Pcup[0] = 1.0
    dPcup[0] = 0.0
    # sin (geocentric latitude) - sin_phi
    z = np.sqrt((1.0 - x) * (1.0 + x))
    schmidtQuasiNorm = np.zeros(NumTerms)

    # schmidtQuasiNorm = (double *) malloc((NumTerms + 1) * sizeof ( double))

    # First, Compute the Gauss-normalized associated Legendre functions
    for n in range(1, nMax + 1):
        for m in range(n + 1):
            index = n * (n + 1) // 2 + m
            if n == m:
                index1 = (n - 1) * n // 2 + m - 1;
                Pcup[index] = z * Pcup[index1];
                dPcup[index] = z * dPcup[index1] + x * Pcup[index1];
            elif n == 1 and m == 0:
                index1 = (n - 1) * n // 2 + m;
                Pcup[index] = x * Pcup[index1];
                dPcup[index] = x * dPcup[index1] - z * Pcup[index1];
            elif n > 1 and n != m:
                index1 = (n - 2) * (n - 1) // 2 + m;
                index2 = (n - 1) * n // 2 + m;
                if m > n - 2:
                    Pcup[index] = x * Pcup[index2]
                    dPcup[index] = x * dPcup[index2] - z * Pcup[index2]
                else:
                    k = float(((n - 1) * (n - 1)) - (m * m)) / float((2 * n - 1) * (2 * n - 3))
                    Pcup[index] = x * Pcup[index2] - k * Pcup[index1];
                    dPcup[index] = x * dPcup[index2] - z * Pcup[index2] - k * dPcup[index1];

    # Compute the ratio between the the Schmidt quasi-normalized associated Legendre
    # functions and the Gauss-normalized version.

    schmidtQuasiNorm[0] = 1.0;
    for n in range(1, nMax + 1):
        index = n * (n + 1) // 2
        index1 = (n - 1) * n // 2
        # for m = 0
        schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * float((2 * n - 1)) / float(n)
        for m in range(1, n + 1):
            index = (n * (n + 1) // 2 + m)
            index1 = (n * (n + 1) // 2 + m - 1)
            schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * np.sqrt(
                float((n - m + 1) * (2 if m == 1 else 1)) / float(n + m))

    # Converts the  Gauss-normalized associated Legendre
    #          functions to the Schmidt quasi-normalized version using pre-computed
    #          relation stored in the variable schmidtQuasiNorm

    for n in range(1, nMax + 1):
        for m in range(n + 1):
            index = (n * (n + 1) // 2 + m);
            Pcup[index] = Pcup[index] * schmidtQuasiNorm[index];
            dPcup[index] = -dPcup[index] * schmidtQuasiNorm[index];
        # The sign is changed since the new WMM routines use derivative with respect to latitude
        # instead of co-latitude

    return [Pcup, dPcup]

def PcupHigh(lat: float, nMax: int):
    """This function evaluates all of the Schmidt-semi normalized associated Legendre
        functions up to degree nMax. The functions are initially scaled by
        10^280 sin^m in order to minimize the effects of underflow at large m
        near the poles (see Holmes and Featherstone 2002, J. Geodesy, 76, 279-299).
        Note that this function performs the same operation as MAG_PcupLow.
        However this function also can be used for high degree (large nMax) models.

        Calling Parameters:
                INPUT
                        nMax:	 Maximum spherical harmonic degree to compute.
                        x:		sin(latitude)

                OUTPUT
                        Pcup:	A vector of all associated Legendgre polynomials evaluated at
                                        x up to nMax. The lenght must by greater or equal to (nMax+1)*(nMax+2)/2.
                  dPcup:   Derivative of Pcup(x) with respect to latitude

                CALLS : none
        Notes:



      Adopted from the FORTRAN code written by Mark Wieczorek September 25, 2005.

      Manoj Nair, Nov, 2009 Manoj.C.Nair@Noaa.Gov

      Change from the previous version
      The prevous version computes the derivatives as
      dP(n,m)(x)/dx, where x = sin(latitude) (or cos(colatitude) ).
      However, the WMM Geomagnetic routines requires dP(n,m)(x)/dlatitude.
      Hence the derivatives are multiplied by sin(latitude).
      Removed the options for CS phase and normalizations.

      Note: In geomagnetism, the derivatives of ALF are usually found with
      respect to the colatitudes. Here the derivatives are found with respect
      to the latitude. The difference is a sign reversal for the derivative of
      the Associated Legendre Functions.

      The derivatives can't be computed for latitude = |90| degrees."""

    #double pm2, pm1, pmm, plm, rescalem, z, scalef;
    #double *f1, *f2, *PreSqr;
    #int k, kstart, m, n, NumTerms;
    x = math.sin(np.deg2rad(lat))
    NumTerms = int((nMax + 1) * (nMax + 2) / 2)

    z = math.sqrt((1.0 - x) * (1.0 + x))


    if z == 0:
        raise ValueError("z becomes error")

    f1 = np.zeros(NumTerms + 1)
    PreSqr = np.zeros(NumTerms + 1)
    f2 = np.zeros(NumTerms + 1)
    Pcup = np.zeros(NumTerms + 1)
    dPcup = np.zeros(NumTerms + 1)

    scalef = 1.0e-308

    for n in range(2 * nMax + 2):
        PreSqr[n] = math.sqrt(n * 1.0)

    k = 2
    for n in range(2, nMax + 1):
        k = k + 1
        f1[k] = float(2 * n - 1) / n
        f2[k] = float(n - 1) / n

        for m in range(1, n - 1):
            k = k + 1
            f1[k] = float(2 * n - 1) / PreSqr[n + m] / PreSqr[n - m]
            f2[k] = PreSqr[n - m - 1] * PreSqr[n + m - 1] / PreSqr[n + m] / PreSqr[n - m]

        k = k + 2

    pm2 = 1.0
    Pcup[0] = 1.0
    dPcup[0] = 0.0

    pm1 = x
    Pcup[1] = pm1
    dPcup[1] = z
    k = 1

    for n in range(2, nMax + 1):
        k = k + n
        plm = f1[k] * x * pm1 - f2[k] * pm2
        Pcup[k] = plm
        dPcup[k] = n * (pm1 - x * plm) / z
        pm2 = pm1
        pm1 = plm

    pmm = PreSqr[2] * scalef
    rescalem = 1.0 / scalef
    kstart = 0

    for m in range(1, nMax):
        rescalem = rescalem * z

        #Calculate Pcup(m, m)
        kstart = int(kstart + m + 1)
        pmm = pmm * PreSqr[2 * m + 1] / PreSqr[2 * m]
        Pcup[kstart] = pmm * rescalem / PreSqr[2 * m + 1]
        dPcup[kstart] = -float(m * x * Pcup[kstart] / z)
        pm2 = pmm / PreSqr[2 * m + 1]
        #Calculate Pcup(m + 1, m)
        k = kstart + m + 1
        pm1 = x * PreSqr[2 * m + 1] * pm2
        Pcup[k] = pm1 * rescalem
        dPcup[k] = float((pm2 * rescalem) * PreSqr[2 * m + 1] - x * (m + 1) * Pcup[k]) / z
        #Calculate Pcup(n, m)

        for n in range(m + 2, nMax + 1):
            k = k + n
            plm = x * f1[k] * pm1 - f2[k] * pm2
            Pcup[k] = plm * rescalem
            dPcup[k] = (PreSqr[n + m] * PreSqr[n - m] * (pm1 * rescalem) - float(n) * x * Pcup[k]) / z
            pm2 = pm1
            pm1 = plm

    #Calculate Pcup(nMax, nMax)

    rescalem = rescalem * z
    kstart = kstart + nMax + 1
    pmm = pmm / PreSqr[2 * nMax]
    Pcup[kstart] = pmm * rescalem
    dPcup[kstart] = -float(nMax * x * Pcup[kstart]) / z

    return Pcup, dPcup
