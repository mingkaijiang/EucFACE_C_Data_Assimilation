#!/usr/bin/env python
"""
DALEC ensemble Kalman filter (EnKF) ...

Reference:
---------
* Williams et al. (2005) An improved analysis of forest carbon dynamics using
  data assimilation. Global Change Biology 11, 89–105.
* Evensen (2003) The Ensemble Kalman Filter: theoretical formulation and
  practical implementation. Ocean Dynamics, 53, 343-367

"""
import sys
import matplotlib.pyplot as plt
import pandas as pd

__author__  = "Martin De Kauwe"
__version__ = "1.0 (07.03.2018)"
__email__   = "mdekauwe@gmail.com"

import numpy as np
import sys

def main(fname):

    obs = None
    met = pd.read_csv(fname)

    # initialise structures
    p = GenericClass()
    c = GenericClass()

    setup_initial_conditions(p, c)

    (A, err_var, err_type, ens_var, q) = setup_holding_matrices_vectors(c)

    A = initialise_ensemble(p, c, A)

    (err_var, err_type, q) = initialise_error_stuff(c, err_var, err_type, q)

    for i in range(len(met)):
        (A, ens_var, q) = forecast(c, p, met, i, A, err_var,
                                       err_type, ens_var, q)

        # Recalcualte model forecast where observations are avaliable
        #if c.nrobs > 0:
        #    analysis(A, c, obs)

        dump_output(c, A)

class GenericClass:
    pass

def setup_holding_matrices_vectors(c):

    # Setup matrix holding ensemble members
    A = np.zeros((c.ndims, c.nrens))

    err_var = np.zeros(c.ndims)
    err_type = np.zeros(c.ndims)
    ens_var = np.zeros(c.ndims)

    # model noise vector
    q = np.zeros((c.ndims, c.nrens))

    return A, err_var, err_type, ens_var, q

def setup_initial_conditions(p, c):

    # dalec model values
    p.t1 = 4.41E-06
    p.t2 = 0.473267
    p.t3 = 0.314951
    p.t4 = 0.434401
    p.t5 = 0.00266518
    p.t6 = 2.06E-06
    p.t7 = 2.48E-03
    p.t8 = 2.28E-02
    p.t9 = 2.65E-06
    p.cf0 = 57.7049
    p.cw0 = 769.863
    p.cr0 = 101.955
    p.cl0 = 40.4494
    p.cs0 = 9896.7

    # acm parameterisation
    p.a0 = 2.155;
    p.a1 = 0.0142;
    p.a2 = 217.9;
    p.a3 = 0.980;
    p.a4 = 0.155;
    p.a5 = 2.653;
    p.a6 = 4.309;
    p.a7 = 0.060;
    p.a8 = 1.062;
    p.a9 = 0.0006;

    # location - oregon
    p.lat = 44.4
    p.sla = 111.

    # timestep
    p.delta_t = 1.0

    # specified time decorrelation length s.
    p.tau = 1.0

    # The factor a should be related to the time step used, eqn 32
    # (i.e. this is zero)
    p.alpha = 1.0 - (p.delta_t / p.tau)

    # setup factor to ensure variance growth over time becomes independent of
    # alpha and delta_timestep (as long as the dynamical model is linear).
    p.rho = setup_stochastic_model_error(p)

    c.nrobs = 0 # Number of observations
    c.ndims = 16  # Dimension of model state
    c.nrens = 200 # Number of ensemble members
    c.max_params = 15
    c.seed = 0

    c.POS_RA = 0
    c.POS_AF = 1
    c.POS_AW = 2
    c.POS_AR = 3
    c.POS_LF = 4
    c.POS_LW = 5
    c.POS_LR = 6
    c.POS_CF = 7
    c.POS_CW = 8
    c.POS_CR = 9
    c.POS_RH1 = 10
    c.POS_RH2 = 11
    c.POS_D = 12
    c.POS_CL = 13
    c.POS_CS = 14
    c.POS_GPP = 15

def initialise_ensemble(p, c, A):

    A[c.POS_RA,:] = 1.0 + np.random.normal(0.0, 0.1 * 1.0, c.nrens)
    A[c.POS_AF,:] = 0.3 + np.random.normal(0.0, 0.1 * 0.3, c.nrens)
    A[c.POS_AW,:] = 0.3 + np.random.normal(0.0, 0.1 * 0.3, c.nrens)
    A[c.POS_AR,:] = 0.3 + np.random.normal(0.0, 0.1 * 0.3, c.nrens)
    A[c.POS_LF,:] = 0.3 + np.random.normal(0.0, 0.1 * 0.3, c.nrens)
    A[c.POS_LW,:] = 0.3 + np.random.normal(0.0, 0.1 * 0.3, c.nrens)
    A[c.POS_LR,:] = 0.3 + np.random.normal(0.0, 0.1 * 0.3, c.nrens)
    A[c.POS_RH1,:] = 0.3 + np.random.normal(0.0, 0.1 * 0.3, c.nrens)
    A[c.POS_RH2,:] = 0.3 + np.random.normal(0.0, 0.1 * 0.3, c.nrens)
    A[c.POS_D,:] = 0.3 + np.random.normal(0.0, 0.1 * 0.3, c.nrens)
    A[c.POS_GPP,:] = 1.0 + np.random.normal(0.0, 0.1 * 1.0, c.nrens)
    A[c.POS_CF,:] = p.cf0 + np.random.normal(0.0, 0.1 * p.cf0, c.nrens)
    A[c.POS_CW,:] = p.cw0 + np.random.normal(0.0, 0.1 * p.cw0, c.nrens)
    A[c.POS_CR,:] = p.cr0 + np.random.normal(0.0, 0.1 * p.cr0, c.nrens)
    A[c.POS_CL,:] = p.cl0 + np.random.normal(0.0, 0.1 * p.cl0, c.nrens)
    A[c.POS_CS,:] = p.cs0 + np.random.normal(0.0, 0.1 * p.cs0, c.nrens)

    return A

def initialise_error_stuff(c, err_var, err_type, q):

    # default error variances for the state vector elements
    err_var[c.POS_RA] = 0.2
    err_var[c.POS_AF] = 0.2
    err_var[c.POS_AW] = 0.2
    err_var[c.POS_AR] = 0.2
    err_var[c.POS_LF] = 0.5
    err_var[c.POS_LW] = 0.5
    err_var[c.POS_LR] = 0.5
    err_var[c.POS_CF] = 0.2
    err_var[c.POS_CW] = 0.2
    err_var[c.POS_CR] = 0.2
    err_var[c.POS_RH1] = 0.2
    err_var[c.POS_RH2] = 0.2
    err_var[c.POS_D] = 0.2
    err_var[c.POS_CL] = 0.2
    err_var[c.POS_CS] = 0.2
    err_var[c.POS_GPP] = 0.2

    err_type[c.POS_RA] = 1
    err_type[c.POS_AF] = 1
    err_type[c.POS_AW] = 1
    err_type[c.POS_AR] = 1
    err_type[c.POS_LF] = 0
    err_type[c.POS_LW] = 0
    err_type[c.POS_LR] = 0
    err_type[c.POS_CF] = 1
    err_type[c.POS_CW] = 1
    err_type[c.POS_CR] = 1
    err_type[c.POS_RH1] = 1
    err_type[c.POS_RH2] = 1
    err_type[c.POS_D] = 1
    err_type[c.POS_CL] = 1
    err_type[c.POS_CS] = 1
    err_type[c.POS_GPP] = 1

    # initial noise matrix
    q = np.random.normal(0.0, 1.0, c.ndims * c.nrens).reshape(c.ndims, c.nrens)

    return err_var, err_type, q

def setup_stochastic_model_error(p):
    """
    Set up stochastic model error according to eqn 42 - Evenson, 2003. This
    ensures that the variance growth over time becomes independent of alpha and
    delta_t (as long as the dynamical model is linear).
    """

    # number of timesteps per time unit
    n = 1.0

    num = (1.0 - p.alpha)**2
    den = n - 2.0 * p.alpha * n * p.alpha**2 + (2.0 * p.alpha)**(n + 1.0)

    return np.sqrt(1.0 / p.delta_t * num / den)

def forecast(c, p, met, i, A, err_var, err_type, ens_var, q):

    A_tmp = np.zeros((c.ndims, c.nrens))
    A_mean = np.zeros(c.ndims)

    # generate model prediction
    lai = np.maximum(0.1, A[c.POS_CF,:]  / p.sla)
    gpp = acm(met, p, lai, i)

    A_tmp[c.POS_GPP,:] = gpp
    A_tmp[c.POS_RA,:] = gpp * p.t2
    A_tmp[c.POS_AF,:] = gpp * p.t3 * (1. - p.t2)
    A_tmp[c.POS_AR,:] = gpp * p.t4 * (1. - p.t2)
    A_tmp[c.POS_AW,:] = gpp * (1. - p.t3 - p.t4) * (1. - p.t2)
    A_tmp[c.POS_LF,:] = A[c.POS_CF,:] * p.t5
    A_tmp[c.POS_LW,:] = A[c.POS_CW,:] * p.t6
    A_tmp[c.POS_LR,:] = A[c.POS_CR,:] * p.t7
    A_tmp[c.POS_RH1,:] = np.exp(0.0693 * met.temp[i]) * A[c.POS_CL,:] * p.t8
    A_tmp[c.POS_RH2,:] = np.exp(0.0693 * met.temp[i]) * A[c.POS_CS,:] * p.t9
    A_tmp[c.POS_D,:] = np.exp(0.0693 * met.temp[i]) * A[c.POS_CL,:] * p.t1
    A_tmp[c.POS_CF,:] = A[c.POS_CF,:] + A[c.POS_AF,:] - A[c.POS_LF,:]
    A_tmp[c.POS_CW,:] = A[c.POS_CW,:] + A[c.POS_AW,:] - A[c.POS_LW,:]
    A_tmp[c.POS_CR,:] = A[c.POS_CR,:] + A[c.POS_AR,:] - A[c.POS_LR,:]
    A_tmp[c.POS_CL,:] = A[c.POS_CL,:] + A[c.POS_LF,:] + A[c.POS_LR,:] - \
                        A[c.POS_RH1,:] - A[c.POS_D,:]
    A_tmp[c.POS_CS,:] = A[c.POS_CS,:] + A[c.POS_D,:] + \
                        A[c.POS_LW,:] - A[c.POS_RH2,:]

    A = A_tmp.copy()

    # Calculate ensemble average
    for i in range(c.ndims):
        A_mean[i] = np.sum(A[i,:]) / float(c.nrens)

    # add the noise model into the ensemble
    for i in range(c.ndims):
        for j in range(c.nrens):
            A[i,j] += np.sqrt(p.delta_t) * p.rho * \
                        np.sqrt(ens_var[i]) * q[i,j]

    # simulate the time evolution of model errors
    for i in range(c.ndims):
        for j in range(c.nrens):
            q_previous_time_step = q[i,j]
            q[i,j] = p.alpha * q_previous_time_step + \
                            np.sqrt(1.0 - p.alpha**2) * \
                            np.random.normal(0.0, 1.0)

    # calculate the error variance
    ens_var = assign_model_errors(c, ens_var, err_var, err_type, A_mean)

    return A, ens_var, q

def assign_model_errors(c, ens_var, err_var, err_type, A_mean):

    for i in range(c.ndims):
        if err_type[i] == 0:
            ens_var[i] = err_var[i]
        else:
            ens_var[i] = err_var[i] * np.fabs(A_mean[i])

    return ens_var

def acm(met, p, lai, i):

    trange = 0.5 * (met.maxt[i] - met.mint[i])
    gs = np.fabs(met.psid[i])**p.a9 / (p.a5 * met.rtot[i] + trange)
    pp = lai * met.nit[i] / gs * p.a0 * np.exp(p.a7 * met.maxt[i])
    qq = p.a2 - p.a3
    ci = 0.5 * (met.ca[i] + qq - pp + np.sqrt((met.ca[i] + qq - pp)**2.0 -\
         4.0 * (met.ca[i] * qq - pp * p.a2)))
    e0 = p.a6 * lai**2 / (lai**2 + p.a8)
    dec = -23.4 * np.cos((360.0 * (met.doy[i]+ 10.0) / 365.0) * \
            np.pi / 180.0) * np.pi / 180.0
    m = np.tan(p.lat * np.pi / 180.0) * np.tan(dec)
    if m >= 1.0:
        dayl = 24.0
    elif m <= -1.0:
        dayl = 0.0
    else:
        dayl = 24.0 * np.arccos(-m) / np.pi
    cps = e0 * met.rad[i] * gs * (met.ca[i] - ci) / \
            (e0 * met.rad[i] + gs * (met.ca[i] - ci))
    gpp = cps * (p.a1 * dayl + p.a4)

    return gpp

def analysis(A, c, obs):

    """
    The standard analysis eqn: A = A + Pe H^T(H Pe H^T + Re)^-1 (D - H A) can
    be reformed using D' = D - HA, Pe = A'(A')^T, Re = YY^T
    (where Y symbolises Gamma) such that it
    becomes A = A + A' A'^T H^T(HA' A'^T H^T + YY^T)^-1 D'
    """

    sig_sum = 0.0
    sig_sum1 = 0.0

    # Minimum of nrobs and nrens
    nrmin = np.minimum(c.nrobs+1, c.nrens)

    I = np.zeros((c.nrens,c.nrens))
    U = np.zeros((c.nrobs,nrmin))
    D = np.zeros((c.nrobs,c.nrens))
    S = np.zeros((c.nrobs,nrmin))
    E = np.zeros((c.nrobs,nrmin))
    H = np.zeros((c.nrobs,c.ndims))
    HA = np.zeros((c.nrobs,c.nrens))
    ES = np.zeros((c.nrobs,c.nrens))
    X1 = np.zeros((nrmin,c.nrens))
    X2 = np.zeros((nrmin,c.nrens))
    X3 = np.zeros((c.nrobs,c.nrens))
    X4 = np.zeros((c.nrens,c.nrens))
    Reps = np.zeros((c.ndims,c.nrobs))
    A_tmp = np.zeros((c.ndims,c.nrens))
    A_dash = np.zeros((c.ndims,c.nrens))

    sig = np.zeros(nrmin)
    D_mean = np.zeros(c.nrobs)
    E_mean = np.zeros(c.nrobs)
    S_mean = np.zeros(c.nrobs)
    HA_mean = np.zeros(c.nrobs)
    A_mean = np.zeros(c.nrobs)

def dump_output(c, A):

    #x = np.sum(A[c.POS_GPP,:])
    #x2 = np.sum(A[c.POS_GPP,:]**2)

    x = np.sum(A[c.POS_CF,:])
    x2 = np.sum(A[c.POS_CF,:]**2)

    ensemble_member_avg = x / float(c.nrens)
    ensemble_member_stdev_error = np.sqrt((x2 - \
                                    (x**2) / float(c.nrens) ) /\
                                    float(c.nrens))

    print(ensemble_member_avg, ensemble_member_stdev_error)

if __name__ == "__main__":

    fname = "data/dalec_drivers.OREGON.no_obs.csv"
    main(fname)
