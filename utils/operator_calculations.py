import numpy as np
import numpy.ma as ma
from scipy.sparse.linalg import eigs
from scipy.sparse import diags,coo_matrix
import deeptime.markov.tools.estimation as dt_estimation
import deeptime.markov.tools.analysis as dt_analysis


def segment_maskedArray(tseries,min_size=50):
    '''
    Segments  time series in case it has missing data
    '''
    if len(tseries.shape)>1:
        mask = ~np.any(tseries.mask,axis=1)
    else:
        mask = ~tseries.mask
    segments = np.where(np.abs(np.diff(np.concatenate([[False], mask, [False]]))))[0].reshape(-1, 2)
    return segments


def get_count_matrix(labels,lag,nstates):
    observable_seqs = ma.compress_rows(ma.vstack([labels[:-lag],labels[lag:]]).T)

    row = observable_seqs[:,0]
    col = observable_seqs[:,1]

    data = np.ones(row.size)
    C = coo_matrix((data, (row, col)), shape=(nstates, nstates))
    # export to output format
    count_matrix = C.tocsr()

    return count_matrix


def largest_connected_submatrix(C):
    return dt_estimation.largest_connected_submatrix(C)


def transition_matrix(labels,lag,return_connected=False):
    nstates = np.max(labels)+1
    count_matrix = get_count_matrix(labels,lag,nstates)
    connected_count_matrix = largest_connected_submatrix(count_matrix)
    P = dt_estimation.transition_matrix(connected_count_matrix)
    if return_connected:
        lcs = dt_estimation.largest_connected_set(count_matrix)
        return lcs,P
    else:
        return P

def get_connected_labels(labels,lcs):
    final_labels = ma.zeros(labels.shape,dtype=int)
    for key in np.argsort(lcs):
        final_labels[labels==lcs[key]]=key+1
    final_labels[final_labels==0] = ma.masked
    final_labels-=1
    return final_labels

def sorted_spectrum(R,k=5,which='LR'):
    eigvals,eigvecs = eigs(R,k=k,which=which)
    sorted_indices = np.argsort(eigvals.real)[::-1]
    return eigvals[sorted_indices],eigvecs[:,sorted_indices]


def get_reversible_transition_matrix(P):
    probs = stationary_distribution(P)
    P_hat = diags(1/probs)*P.transpose()*diags(probs)
    R=(P+P_hat)/2
    return R


def stationary_distribution(P):
    probs = dt_analysis.stationary_distribution(P)
    return probs


def get_entropy(labels):
    P = transition_matrix(labels,1)
    probs = stationary_distribution(P)
    logP = P.copy()
    logP.data = np.log(logP.data)
    return (-diags(probs).dot(P.multiply(logP))).sum()


def simulate(P,state0,iters):
    '''
    Monte Carlo simulation of the markov chain characterized by the matrix P
    state0: initial system
    iters: number of iterations of the simulation
    '''
    states = np.zeros(iters,dtype=int)
    states[0]=state0
    state=state0
    for k in range(1,iters):
        new_state = np.random.choice(np.arange(P.shape[1]),p=list(np.hstack(P[state,:].toarray())))
        state=new_state
        states[k]=state
    return states


def state_lifetime(states,tau):
    '''
    Get distribution of lifetimes of each state in states
    tau is the sampling time of the states
    '''
    durations=[]
    for state in np.sort(np.unique(states.compressed())):
        gaps = states==state
        gaps_boundaries = np.where(np.abs(np.diff(np.concatenate([[False], gaps, [False]]))))[0].reshape(-1, 2)
        durations.append(np.hstack(np.diff(gaps_boundaries))*tau)
    return durations
