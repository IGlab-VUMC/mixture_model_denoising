## Functions for running negative binomial mixture model for Wasdin et. al. 2023

# Loading in necessary packages
import numpy as np
import pandas as pd
import scipy.optimize
import scipy.stats as st

import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Seq import Seq
import Levenshtein as Lev
import warnings

def separate_vrc01(df):
    '''Input the LIBRA-seq output dataframe. Returns a non-VRC01 and VRC01 dataframe (in that order).'''
    vrc01_cdr3 = 'ACTAGGGGAAAAAACTGTGATTACAATTGGGACTTCGAACAC'
    lseq_cdr3s = pd.DataFrame(df['CDR3_IMGT.H'])
    vrc01_list = []
    for i, row in lseq_cdr3s.iterrows():
        if len(row.values[0]) == len(vrc01_cdr3):

            # Using overlap of 0.99 as criteria for VRC01
            hd = (Lev.distance(row.values[0], vrc01_cdr3))
            if hd/len(row.values[0]) < 0.01:
                vrc01_list.append(i)
    vrc01_df = df.loc[vrc01_list]

    non_vrc01_list = [i for i in df.index if i not in vrc01_list]
    lseq_non_vrc01 = df.loc[non_vrc01_list]

    return lseq_non_vrc01, vrc01_df

def convert_params_statsmodels(model_noise_nb):
    '''Convert mixed distribution parameters to be n/p for use with Scipy nbinom'''
    mu_noise = np.exp(model_noise_nb.params[0])
    n_noise = 1/model_noise_nb.params[1]
    p_noise = n_noise/(n_noise + mu_noise)
    beta_noise = (1/p_noise)-1 

    return mu_noise, n_noise, p_noise, beta_noise


def convert_params_mle_to_scipy(param_array):
    '''Convert mixed distribution parameters to be n/p for use with Scipy nbinom'''
    n_signal = param_array[0]
    p_signal = 1/(1 + param_array[1])
    n_noise = param_array[2]
    p_noise= 1/(1 + param_array[3])
    w = param_array[4]

    return n_signal, p_signal, n_noise, p_noise, w

# The following functions were taken and modified from:
# http://bebi103.caltech.edu.s3-website-us-east-1.amazonaws.com/2019a/content/lessons/lesson_09/mixture_models.html

def log_like_mix(alpha1, b1, alpha2, b2, w, n):
    """Log-likeihood of binary Negative Binomial mixture model."""
    # Fix nonidentifieability be enforcing values of w
    if w < 0 or w > 1:
        return -np.inf

    # Physical bounds on parameters
    if alpha1 < 0 or alpha2 < 0 or b1 < 0 or b2 < 0:
        return -np.inf

    logx1 = st.nbinom.logpmf(n, alpha1, 1/(1+b1))
    logx2 = st.nbinom.logpmf(n, alpha2, 1/(1+b2))

    # Multipliers for log-sum-exp
    lse_coeffs = np.tile([w, 1-w], [len(n), 1]).transpose()

    # log-likelihood for each measurement
    log_likes = scipy.special.logsumexp(np.vstack([logx1, logx2]), axis=0, b=lse_coeffs)

    return np.sum(log_likes)
def log_like_iid_nbinom(params, n):
    """Log likelihood for i.i.d. NBinom measurements, parametrized
    by alpha, b=1/beta."""
    alpha, b = params

    if alpha <= 0 or b <= 0:
        return -np.inf

    return np.sum(st.nbinom.logpmf(n, alpha, 1/(1+b)))


def mle_iid_nbinom(n):
    """Perform maximum likelihood estimates for parameters for i.i.d.
    NBinom measurements, parametrized by alpha, b=1/beta"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        res = scipy.optimize.minimize(
            fun=lambda params, n: -log_like_iid_nbinom(params, n),
            x0=np.array([3, 3]),
            args=(n,),
            method='Powell'
        )

    if res.success:
        return res.x
    else:
        raise RuntimeError('Convergence failed with message', res.message)

def initial_guess_mix(n, w_guess):
    """Generate initial guess for mixture model."""
    n_low = n[n < np.percentile(n, 100*w_guess)]
    n_high = n[n >= np.percentile(n, 100*w_guess)]

    alpha1, b1 = mle_iid_nbinom(n_low)
    alpha2, b2 = mle_iid_nbinom(n_high)

    return alpha1, b1, alpha2, b2


# custom_guess = (n_noise, p_noise, p_signal, n_signal1)
def mle_mix(n, w_guess):
    """Obtain MLE estimate for parameters for binary mixture
    of Negative Binomials."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        res = scipy.optimize.minimize(
            fun=lambda params, n: -log_like_mix(*params, n),
            x0=[*initial_guess_mix(n, w_guess), w_guess],
            # x0=[*custom_guess, w_guess],

            args=(n,),
            method='Powell',
            tol=1e-8,
        )

    if res.success:
        return res.x
    else:
        raise RuntimeError('Convergence failed with message', res.message)
    
def custom_mle_mix(n, w_guess, initials):
    """Obtain MLE estimate for parameters for binary mixture
    of Negative Binomials."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        res = scipy.optimize.minimize(
            fun=lambda params, n: -log_like_mix(*params, n),
            # x0=[*initial_guess_mix(n, w_guess), w_guess],
            x0=[*initials, w_guess],


            args=(n,),
            method='Powell',
            tol=1e-4,
        )

    if res.success:
        return res
    else:
        return res
        # raise RuntimeError('Convergence failed with message', res.message)



# Back to my functions

def calculate_component_probabilities(x_signal, umi, popt):

    probdf = pd.DataFrame(x_signal).copy()
    n_signal, p_signal, n_noise, p_noise, weight = convert_params(popt)

    probdf['pmf'] = (st.nbinom.pmf(probdf[umi], n_signal, p_signal))
    probdf['pmf_noise'] = (st.nbinom.pmf(probdf[umi], n_noise, p_noise))

    probdf['pmf'] = (st.nbinom.pmf(probdf[umi], n_signal, p_signal))
    probdf['pmf_noise'] = (st.nbinom.pmf(probdf[umi], n_noise, p_noise))

    pxa = (weight)*(probdf['pmf'])
    pxb = (1-weight)*(probdf['pmf_noise'])
    Px = pxa + pxb
    Pa = pxa/Px
    Pb = pxb/Px
    probdf['probsA'] = Pa
    probdf['probsB'] = Pb

    if st.nbinom.median(n_signal, p_signal) < st.nbinom.median(n_noise, p_noise):
        probdf['probsB'] = Pa
        probdf['probsA'] = Pb
    return probdf