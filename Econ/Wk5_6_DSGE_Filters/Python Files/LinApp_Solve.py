"""
Python module for using the method outlined by Uhlig (1997) to solve a
log-linearized RBC model for policy functions.

Original adaptation of MATLAB code done by Spencer Lyon in May 2012

Additional work has been done by Chase Coleman
"""
import scipy as sp
import numpy as np
from numpy import hstack, vstack, zeros, dot, eye, kron
from scipy import linalg as la
from numpy import linalg as npla


def _nullSpaceBasis(A):
    """
    This funciton will find the basis of the null space of the matrix A.

    Parameters
    ----------
    A : array_like, dtype=float
        The matrix you want the basis for

    Returns
    -------
    vecs : array_like, dtype=float
        A numpy matrix containing the vectors as row vectors.

    Notes
    -----
    If A is an empty matrix, an empty matrix is returned.

    """
    if A:
        U, s, Vh = la.svd(A)
        vecs = np.array([])
        toAppend = A.shape[1] - s.size
        s = np.append(s, zeros((1, toAppend)))
        for i in range(0, s.size):
            if s[i] == 0:
                vecs = Vh[-toAppend:, :]
        if vecs.size == 0:
            vecs = zeros((1, A.shape[1]))
        return np.mat(vecs)
    else:
        return zeros((0, 0))

Uhlig_Solve(AA=None, BB=None, CC=None, DD=None, FF=None, GG=None, HH=None,
              JJ=None, KK=None, LL=None, MM=None, NN=None):
    """
    This function mimics the behavior of Harald Uhlig's solve.m and
    calc_qrs.m files in Uhlig's toolkit.

    In order to use this function, the user must have log-linearized the
    model they are dealing with to be in the following form (assume that
    y corresponds to the model's "jump variables", z represents the
    exogenous state variables and x is for endogenous state variables.
    nx, ny, nz correspond to the number of variables in each category.)
    The inputs to this function are the matrices found in the following
    equations.

    .. math::

        Ax_t + Bx_t-1 + Cy_t + Dz_t = 0

        E\{Fx_{t+1} + Gx_t + Hx_{t-1} + Jy_{t+1} +
           Ky_t + Lz_{t+1} Mz_t \} = 0

    The purpose of this function is to find the recursive equilibrium
    law of motion defined by the following equations.

    .. math::

        X_t = PX_{t-1} + Qz_t

        Y_t = RY_{t-1} + Sz_t

    Following outline given in Uhhlig (1997), we solve for :math:`P` and
    :math:`Q` using the following set of equations:

    .. math::

        FP^2 + Gg + H =0

        FQN + (FP+G)Q + (LN + M)=0

    Once :math:`P` and :math:`Q` are known, one ca solve for :math:`R`
    and :math:`S` using the following equations:

    .. math::

        R = -C^{-1}(AP + B)

        S = - C^{-1}(AQ + D)

    Parameters
    ----------
    AA : array_like, dtype=float, shape=(ny, nx)
        The matrix represented above by :math:`A`. It is the matrix of
        derivatives of the Y equations with repsect to :math:`X_t`
    BB : array_like, dtype=float, shape=(ny, nx)
        The matrix represented above by :math:`B`. It is the matrix of
        derivatives of the Y equations with repsect to
        :math:`X_{t-1}`.
    CC : array_like, dtype=float, shape=(ny, ny)
        The matrix represented above by :math:`C`. It is the matrix of
        derivatives of the Y equations with repsect to :math:`Y_t`
    DD : array_like, dtype=float, shape=(ny, nz)
        The matrix represented above by :math:`C`. It is the matrix of
        derivatives of the Y equations with repsect to :math:`Z_t`
    FF : array_like, dtype=float, shape=(nx, nx)
        The matrix represetned above by :math:`F`. It the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`X_{t+1}`
    GG : array_like, dtype=float, shape=(nx, nx)
        The matrix represetned above by :math:`G`. It the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`X_t`
    HH : array_like, dtype=float, shape=(nx, nx)
        The matrix represetned above by :math:`H`. It the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`X_{t-1}`
    JJ : array_like, dtype=float, shape=(nx, ny)
        The matrix represetned above by :math:`J`. It the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`Y_{t+1}`
    KK : array_like, dtype=float, shape=(nx, ny)
        The matrix represetned above by :math:`K`. It the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`Y_t`
    LL : array_like, dtype=float, shape=(nx, nz)
        The matrix represetned above by :math:`L`. It the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`Z_{t+1}`
    MM : array_like, dtype=float, shape=(nx, nz)
        The matrix represetned above by :math:`M`. It the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`Z_t`
    NN : array_like, dtype=float, shape=(nz, nz)
        The autocorrelation matrix for the exogenous state
        vector z.

    Returns
    -------
    P : array_like, dtype=float, shape=(nx, nx)
        The matrix :math:`P` in the law of motion for endogenous state
        variables described above.
    Q : array_like, dtype=float, shape=(nx, nz)
        The matrix :math:`P` in the law of motion for endogenous state
        variables described above.
    R : array_like, dtype=float, shape=(ny, nx)
        The matrix :math:`P` in the law of motion for endogenous state
        variables described above.
    S : array_like, dtype=float, shape=(ny, nz)
        The matrix :math:`P` in the law of motion for endogenous state
        variables described above.

    References
    ----------
    .. [1] Uhlig, H. (1999): "A toolkit for analyzing nonlinear dynamic
       stochastic models easily," in Computational Methods for the Study
       of Dynamic Economies, ed. by R. Marimon, pp. 30-61. Oxford
       University Press.

    """
    #The original coding we did used the np.matrix form for our matrices so we
    #make sure to set our inputs to numpy matrices.
    AA = np.matrix(AA)
    BB = np.matrix(BB)
    CC = np.matrix(CC)
    DD = np.matrix(DD)
    FF = np.matrix(FF)
    GG = np.matrix(GG)
    HH = np.matrix(HH)
    JJ = np.matrix(JJ)
    KK = np.matrix(KK)
    LL = np.matrix(LL)
    MM = np.matrix(MM)
    NN = np.matrix(NN)
    #Tolerance level to use
    TOL = .000001

    # Here we use matrices to get pertinent dimensions.
    nx = FF.shape[1]
    l_equ = CC.shape[0]
    ny = CC.shape[1]
    nz = LL.shape[1]

    k_exog = min(NN.shape)

    # The following if and else blocks form the
    # Psi, Gamma, Theta Xi, Delta mats
    if l_equ == 0:
        if CC.any():
            # This blcok makes sure you don't throw an error with an empty CC.
            CC_plus = la.pinv(CC)
            CC_0 = _nullSpaceBasis(CC.T)
        else:
            CC_plus = np.mat([])
            CC_0 = np.mat([])
        Psi_mat = FF
        Gamma_mat = -GG
        Theta_mat = -HH
        Xi_mat = np.mat(vstack((hstack((Gamma_mat, Theta_mat)),
                        hstack((eye(nx), zeros((nx, nx)))))))
        Delta_mat = np.mat(vstack((hstack((Psi_mat, zeros((nx, nx)))),
                           hstack((zeros((nx, nx)), eye(nx))))))

    else:
        CC_plus = la.pinv(CC)
        CC_0 = _nullSpaceBasis(CC.T)
        Psi_mat = vstack((zeros((l_equ - ny, nx)), FF \
                            - dot(dot(JJ, CC_plus), AA)))
        if CC_0.size == 0:
            # This block makes sure you don't throw an error with an empty CC.
            Gamma_mat = vstack((dot(CC_0, AA), dot(dot(JJ, CC_plus), BB) \
                        - GG + dot(dot(KK, CC_plus), AA)))
            Theta_mat = vstack((dot(CC_0, AA), dot(dot(KK, CC_plus), BB) \
                                    - HH))
        else:
            Gamma_mat = dot(dot(JJ, CC_plus), BB) - GG +\
                               dot(dot(KK, CC_plus), AA)
            Theta_mat = dot(dot(KK, CC_plus), BB) - HH
        Xi_mat = vstack((hstack((Gamma_mat, Theta_mat)), \
                            hstack((eye(nx), zeros((nx, nx))))))
        Delta_mat = vstack((hstack((Psi_mat, np.mat(zeros((nx, nx))))),\
                                hstack((zeros((nx, nx)), eye(nx)))))

    # Now we need the generalized eigenvalues/vectors for Xi with respect to
    # Delta. That is eVals and eVecs below.

    eVals, eVecs = la.eig(Xi_mat, Delta_mat)
    if npla.matrix_rank(eVecs) < nx:
        print('Error: Xi is not diagonalizable, stopping')

    # From here to line 158 we Diagonalize Xi, form Lambda/Omega and find P.
    else:
        Xi_sortabs = np.sort(abs(eVals))
        Xi_sortindex = np.argsort(abs(eVals))
        Xi_sortedVec = np.array([eVecs[:, i] for i in Xi_sortindex]).T
        Xi_sortval = Xi_sortabs
        Xi_select = np.arange(0, nx)
        if np.imag(Xi_sortedVec[nx - 1]).any():
            if (abs(Xi_sortval[nx - 1] - sp.conj(Xi_sortval[nx])) < TOL):
                drop_index = 1
                cond_1 = (abs(np.imag(Xi_sortval[drop_index])) > TOL)
                cond_2 = drop_index < nx
                while cond_1 and cond_2:
                    drop_index += 1
                if drop_index >= nx:
                    print('There is an error. Too many complex eigenvalues.'
                          +' Quitting')
                else:
                    print('droping the lowest real eigenvalue. Beware of' +
                          ' sunspots')
                    Xi_select = np.array([np.arange(0, drop_index - 1),\
                                          np.arange(drop_index + 1, nx)])

        # Here Uhlig computes stuff if user chose "Manual roots" I skip it.
        if max(abs(Xi_sortval[Xi_select])) > 1 + TOL:
            print('It looks like we have unstable roots. This might not work')
        if abs(max(abs(Xi_sortval[Xi_select])) - 1) < TOL:
            print('Check the model to make sure you have a unique steady' +
                  ' state we are having problems with convergence.')
        Lambda_mat = np.diag(Xi_sortval[Xi_select])
        Omega_mat = Xi_sortedVec[nx:2 * nx, Xi_select]
        #Omega_mat = sp.reshape(Omega_mat,\
        #        (math.sqrt(Omega_mat.size),math.sqrt(Omega_mat.size)))
        if npla.matrix_rank(Omega_mat) < nx:
            print("Omega matrix is not invertible, Can't solve for P")
        else:
            PP = dot(dot(Omega_mat, Lambda_mat), la.inv(Omega_mat))
            PP_imag = np.imag(PP)
            PP = np.real(PP)
            if (sum(sum(abs(PP_imag))) / sum(sum(abs(PP))) > .000001).any():
                print("A lot of P is complex. We will continue with the" +
                      " real part and hope we don't lose too much information")

    # The code from here to the end was from he Uhlig file cacl_qrs.m.
    # I think for python it fits better here than in a separate file.

    # The if and else below make RR and VV depending on our model's setup.
    if l_equ == 0:
        RR = zeros((0, nx))
        VV = hstack((kron(NN.T, FF) + kron(eye(k_exog), \
            (dot(FF, PP) + GG)), kron(NN.T, JJ) + kron(eye(k_exog), KK)))

    else:
        RR = - dot(CC_plus, (dot(AA, PP) + BB))
        VV = sp.vstack((hstack((kron(eye(k_exog), AA), \
                        kron(eye(k_exog), CC))), hstack((kron(NN.T, FF) +\
                        kron(eye(k_exog), dot(FF, PP) + dot(JJ, RR) + GG),\
                        kron(NN.T, JJ) + kron(eye(k_exog), KK)))))

    # Now we use LL, NN, RR, VV to get the QQ, RR, SS matrices.
    if (npla.matrix_rank(VV) < k_exog * (nx + ny)):
        print("Sorry but V is not invertible. Can't solve for Q and S")
    else:
        LL = sp.mat(LL)
        NN = sp.mat(NN)
        LLNN_plus_MM = dot(LL, NN) + MM

        if DD.any():
            impvec = vstack([DD.T, np.reshape(LLNN_plus_MM,
                                              (nx * k_exog, 1), 'F')])
        else:
            impvec = np.reshape(LLNN_plus_MM, (nx * k_exog, 1), 'F')

        QQSS_vec = np.matrix(la.solve(-VV, impvec))

        if (max(abs(QQSS_vec)) == sp.inf).any():
            print("We have issues with Q and S. Entries are undefined." +
                  " Probably because V is no inverible.")

        QQ = np.reshape(np.matrix(QQSS_vec[0:nx * k_exog, 0]),
                        (nx, k_exog), 'F')

        SS = np.reshape(QQSS_vec[(nx * k_exog):((nx + ny) * k_exog), 0],\
                        (ny, k_exog), 'F')

        #Build WW - we don't use this, but Uhlig defines it so we do too.
        WW = sp.vstack((
        hstack((eye(nx), zeros((nx, k_exog)))),
        hstack((dot(RR, la.pinv(PP)), (SS - dot(dot(RR, la.pinv(PP)), QQ)))),
        hstack((zeros((k_exog, nx)), eye(k_exog)))))

    return PP, QQ, RR, SS
