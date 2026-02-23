import numpy as np
import numpy.ma as ma
import operator_calculations as op_calc
from scipy.signal import find_peaks
from sklearn.cluster import KMeans
from joblib import Parallel, delayed
import skfuzzy as fuzz


def calculate_coherence(X,c,inv_measure, P):
    labels = np.zeros(len(X),dtype=int)
    labels[X<=c] = 1
    rho_sets = [(inv_measure[labels==idx]@(P[labels==idx,:][:,labels==idx])).sum()/inv_measure[labels==idx].sum()
                  for idx in range(2)]
    return rho_sets


def calculate_measures(X,c,inv_measure,P):
    labels = np.zeros(len(X),dtype=int)
    labels[X<=c] = 1
    measures = [(inv_measure[labels==idx]).sum() for idx in range(2)]
    return measures


def optimal_partition(phi2,inv_measure,P,return_rho = True):
    '''
    Find coherent sets so as to maximize the measure of coherence, rho
    phi2: first nontrivial eigenfunction of R
    inv_measure: steady-state distribution
    P: Markov transition matrix
    '''
    X = phi2
    c_range = np.sort(phi2)[1:-1]
    rho_c = np.zeros(len(c_range))
    rho_sets = Parallel(n_jobs=-1)(delayed(calculate_coherence)(X,c,inv_measure,P) for c in c_range)
    measures = Parallel(n_jobs=-1)(delayed(calculate_measures)(X,c,inv_measure,P) for c in c_range)
    rho_sets = np.vstack(rho_sets)
    measures = np.hstack(measures)
    rho_c = np.min(rho_sets,axis=1)
    peaks, heights = find_peaks(rho_c, height=0.5)
    if len(peaks)==0:
        print('No prominent coherent set')
        return None
    else:
        idx = peaks[np.argmax(heights['peak_heights'])]

        c_opt = c_range[idx]
        kmeans_labels = np.zeros(len(X),dtype=int)
        kmeans_labels[X<=c_opt] = 1

        if return_rho:
            return c_range,rho_sets,measures,idx,kmeans_labels
        else:
            return kmeans_labels


def get_distorted_eigfs(P_ensemble,num_eigfs_max,tau=1,dt=1,commute=True):
    """
    Rescale the eigenvectors of P such that the split location falls at zero
    """
    inv_measure = op_calc.stationary_distribution(P_ensemble)
    R = op_calc.get_reversible_transition_matrix(P_ensemble)
    eigvals,eigvecs = op_calc.sorted_spectrum(R,k=num_eigfs_max+1)
    sorted_indices = np.argsort(eigvals.real)[::-1]
    eigvals = eigvals[sorted_indices][1:].real
    eigvals[np.abs(eigvals-1)<1e-12] = np.nan
    eigvals[eigvals<1e-12] = np.nan

    t_imp = -(tau*dt)/np.log(np.abs(eigvals))

    eigfunctions = eigvecs.real/np.linalg.norm(eigvecs.real,axis=0)
    distorted_eigfs = np.zeros((eigfunctions.shape[0], eigfunctions.shape[1]-1))
    for i in range(1,eigfunctions.shape[1]):
        print(i)
        phi = eigfunctions[:,i]
        c_range,rho_sets,measures,split_idx,coh_labels  = optimal_partition(phi,inv_measure,P_ensemble,return_rho=True)

        split_idx = split_idx+1
        sort_range = np.sort(phi)
        if np.abs(sort_range[split_idx+1] - sort_range[split_idx]) > np.abs(sort_range[split_idx-1] - sort_range[split_idx]):
            next_phi_loc = split_idx+1
            midpoint =  0.5*(sort_range[split_idx] +  sort_range[next_phi_loc])
            den_min = np.abs((phi.min() - midpoint))
            den_max = np.abs((phi.max() - midpoint))

            del_x_min = (sort_range[split_idx] - midpoint)/den_min
            del_x_max = (sort_range[next_phi_loc]-midpoint)/den_max

            neg_range = np.linspace(-1,del_x_min, len(sort_range[0:split_idx+1]))
            pos_range = np.linspace(del_x_max,1,len(sort_range[next_phi_loc:]))
        else:
            next_phi_loc = split_idx-1
            midpoint =  0.5*(sort_range[split_idx] +  sort_range[next_phi_loc])
            den_min = np.abs((phi.min() - midpoint))
            den_max = np.abs((phi.max() - midpoint))

            del_x_min = (sort_range[next_phi_loc] - midpoint)/den_min
            del_x_max = (sort_range[split_idx]-midpoint)/den_max

            neg_range = np.linspace(-1,del_x_min, len(sort_range[0:next_phi_loc+1]))
            pos_range = np.linspace(del_x_max,1,len(sort_range[split_idx:]))

        distort_r = np.hstack([neg_range,pos_range])
        distort = np.zeros(phi.shape)

        distort = np.zeros(phi.shape)
        pos = [np.where(phi == a)[0][0] for a in np.sort(phi)]

        for j in range(phi.shape[0]):
            distort[pos[j]] = distort_r[j]
        if commute==True:
            distorted_eigfs[:,i-1] = distort*np.sqrt(t_imp[i-1]/2)
        else:
            distorted_eigfs[:,i-1] = distort

    return distorted_eigfs, inv_measure


def kmeans_propagation(distorted_eigfs,n_eigf,state_to_split,hard_assign):
    """
    Perform Kmeans based coarse-graining on the required number of eigenvectors to discover substates
    """
    X = distorted_eigfs[:,:n_eigf]
    idx = np.where(hard_assign == state_to_split)[0]
    km = KMeans(n_clusters=2,random_state=123,n_init=1000).fit(X[idx,:])
    hard_assign[idx] = km.labels_ +np.max(hard_assign)+10

    final_hard_assign = np.zeros(hard_assign.shape,dtype=int)
    for new_idx,label in enumerate(np.sort(np.unique(hard_assign))):
        final_hard_assign[hard_assign==label]=new_idx

    return final_hard_assign


def fuzzy_propagation(ev,n_eigf,indices,g_to_split,hard_assign,rb):
    """
    Perform fuzzy c-means based coarse-graining on the required number of eigenvectors to discover substates
    and return the posterior probabilities.
    """
    X = ev[indices[g_to_split],:n_eigf]
    c, cr,_,_,_,_,_= fuzz.cmeans(X.T, 2,m=2,error=0.005,maxiter=5000)

    cr_preds = np.zeros((len(np.hstack(indices)),2))
    for i in np.unique(hard_assign):
        if i == g_to_split:
            cr_preds[indices[i],:] = cr.T
        else:
            X = ev[indices[i],:n_eigf]
            cr_pred,_,_,_,_,_ = fuzz.cmeans_predict(X.T, c,2,error=0.005,maxiter=5000)
            cr_preds[indices[i],:] = cr_pred.T

    new_rb = np.multiply(cr_preds, np.vstack([rb[:,g_to_split],rb[:,g_to_split]]).T)
    not_g_to_split = np.delete(np.arange(len(indices)),g_to_split)
    hard_assign[indices[g_to_split]] = np.argmax(cr,axis=0) + np.max(hard_assign) + 1

    final_rb = np.zeros((rb.shape[0],len(np.unique(hard_assign))))
    final_rb[:,:len(np.unique(hard_assign))-2] = rb[:,not_g_to_split]
    final_rb[:,len(np.unique(hard_assign))-2:] = new_rb

    final_hard_assign = np.zeros(hard_assign.shape,dtype=int)

    for new_idx,label in enumerate(np.sort(np.unique(hard_assign))):
        final_hard_assign[hard_assign==label]=new_idx

    return final_hard_assign,final_rb


def relative_fuzzy_partition(distorted_eigfs,inv_measure,cgs=10,response=False):

    """
    A hierarchical subdivision fo dynamics on the distorted eigenvectors

    Parameters:
    distorted_eigfs: distorted eigenfunctions of the transfer operator approximation
    inv_measure: invariant measure of the whole transfer operator
    num_eigfs_max: Maximum number of eigenfunctions to use (?)
    n_groups: The number of groups to split into. If None, then split will run until scale separation stops

    Returns:
    labels_tree: The cluster labels after a hard assign
    responsibilities: The posterior P(C|X)
    """
    X = distorted_eigfs[:,:1]
    print('CG = 2',flush=True)

    c, cr,_,_,_,_,_= fuzz.cmeans(distorted_eigfs[:,0][:,np.newaxis].T, 2,m=2,error=0.0005,maxiter=1000)
    hard_assign = np.argmax(cr,axis=0)

    labels_tree=np.zeros((cgs-1,len(hard_assign)),dtype=int)

    respons = np.zeros((distorted_eigfs.shape[0],2))
    respons = cr.T
    if response ==  True:
        response_tree = []
        response_tree.append(cr.T)

    labels_tree[0,:] = hard_assign

    measures_iter = []
    measures_iter.append([(inv_measure[hard_assign==g]).sum() for g in np.unique(hard_assign)])

    for cg in range(3, cgs+1):
        print('CG = {}'.format(cg),flush=True)
        measures = [(inv_measure[hard_assign==g]).sum() for g in np.unique(hard_assign)]
        indices_groups = [np.where(hard_assign==g)[0] for g in np.unique(hard_assign)]
        measures_sort_idx = np.argsort(measures)[::-1]

        for m in measures_sort_idx:
            idx = np.where(hard_assign == m)[0]
            if len(idx) > 1:
                g_to_split = m
                break

        n_eigf = int(np.ceil(np.log2(cg)))
        hard_assign,respons = fuzzy_propagation(distorted_eigfs,n_eigf,indices_groups,g_to_split,hard_assign,respons)
        labels_tree[cg-2,:] = np.copy(hard_assign)
        if response == True:
            response_tree.append(respons)

        measures_iter.append([(inv_measure[hard_assign==g]).sum() for g in np.unique(hard_assign)])
    if response == True:
        return labels_tree,response_tree
    else:
        return labels_tree
