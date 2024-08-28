import math
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.metrics import precision_recall_fscore_support, precision_recall_curve, auc
import matplotlib.pyplot as plt
# from random import seed
# from random import randint
import random
import pickle


def get_biggest_complex_distances(pep_results, mesure_type, keep_monomers=False):
    mesures = dict()
    for pep_idx in pep_results.keys():
        best_s = None
        best_n_chains = 0
        #select the structure with most chians
        for s in pep_results[pep_idx]:
            n_chains = pep_results[pep_idx][s]['chains']
            if n_chains > best_n_chains:
                best_n_chains = n_chains
                best_s = s
        
        if best_n_chains==1:
            if keep_monomers:
                mesures[pep_idx] = {'pdb':pdb, 'distance':1000}     #there is no interface!
            else:
                continue
        
        #get a random segment for now
        else:
            pdb = best_s.split('_')[2]
            mesures[pep_idx] = {'pdb':pdb, 'distance': pep_results[pep_idx][best_s][mesure_type][0]}        
        
    return mesures


# def get_n_chains_complexes_distances(pep_results, mesure_type, n_c):
#     mesures = dict()
#     for pep_idx in pep_results.keys():
#         best_s = None
#         best_n_chains = 0
#         #select the structure with most chians
#         for s in pep_results[pep_idx]:
#             n_chains = pep_results[pep_idx][s]['chains']
#             if n_chains==n_c:
#                 best_n_chains = n_chains
#                 best_s = s
#                 break
        
#         if best_n_chains!=n_c:
#             continue
        
#         #get a random segment for now
#         pdb = best_s.split('_')[2]
#         mesures[pep_idx] = {'pdb':pdb, 'distance': pep_results[pep_idx][best_s][mesure_type][0]}       
        
#     return mesures



def get_prediction_scores(peptides, mesures, cut=0):
    distances = list()
    labels = list()

    q = list()
    peptide_list = list()
    protein = list()
    q_scores = list()
    
    # mq = list()
    # mq_scores = list()

    for idx in mesures.keys():
        dist = mesures[idx]['distance']
        #raw data
        distances.append(dist)
        q.append(peptides[idx].q)
        peptide_list.append(peptides[idx].peptide)
        protein.append(peptides[idx].protein)
        #mq.append(peptides[idx].median_q)
        #mq.append(peptides[idx].q)

        if dist <= cut:
            labels.append(+1)
        elif dist > cut:
            labels.append(-1)
        else:
            print('Warning! Unexpected distance: {}'.format(dist))
            return

        q_scores.append(-math.log(peptides[idx].q, 2))
        #mq_scores.append(-math.log(peptides[idx].median_q, 2))
        #mq_scores.append(-math.log(peptides[idx].q, 2))

    return {'peptide' : peptide_list, 'protein': protein, 'distances': distances, 'labels':labels, 'q':q, 'q_scores':q_scores} #, 'mq':mq, 'mq_scores':mq_scores}


def show_data_balance(results):
    n_pos = 0
    n_neg = 0
    for x in results['labels']:
        if x>0:
            n_pos+=1
        elif x<0:
            n_neg+=1
    print('{} positives'.format(n_pos))
    print('{} negatives'.format(n_neg))
    return


def compute_scores(results):
    precision, recall, thresholds = precision_recall_curve(results['labels'], results['q_scores'])
    auc_precision_recall = auc(recall, precision)
    print('   *AUC ROC {} \n   *AUC PR {}'.format(roc_auc_score(results['labels'], results['q_scores']), auc_precision_recall))

    return


def store_scores(results):
    return roc_auc_score(results['labels'], results['q_scores'])




##############################################################################
############################## BEST ACCURACY #################################
##############################################################################


def get_complexes_list(clusters, pdb_distances, peptides, mesure, cut):
    complexes = dict()
    for p in clusters:
        complexes[p] = dict()
        for idx in clusters[p]:
            if idx not in pdb_distances.keys():
                continue      
            if peptides[idx].q<0.05:
                sign = True
            else:
                sign = False
            
            for comp in pdb_distances[idx]:
                if pdb_distances[idx][comp][mesure][0] is None:
                    continue
                
                if comp not in complexes[p].keys():
                    complexes[p][comp] = list()
            
                if pdb_distances[idx][comp][mesure][0]<cut and sign:
                    complexes[p][comp].append(True)
                elif pdb_distances[idx][comp][mesure][0]>cut and not sign:
                    complexes[p][comp].append(True)
                else:
                    complexes[p][comp].append(False)

    #take only the complexes with at least one peptide mapped         
    complexes = {c:complexes[c] for c in complexes if len(complexes[c])>0}        
    return complexes


def get_best_complexes(complexes, clusters):
    b_c = dict()
    for p in complexes:
        best_complex = ''
        best_accuracy = -1
        for c in complexes[p]:
            norm = len(complexes[p][c])
            trues = 0
            for i in complexes[p][c]:
                if i:
                    trues+=1
            accuracy = trues/norm
            if accuracy>best_accuracy:
                best_complex = c
                best_accuracy = accuracy
            
        b_c[p] = best_complex  

    #sort the complexes according to the peptide id
    best_complexes = dict()
    for c in clusters.keys():
        if c not in b_c.keys():
            continue
        for idx in clusters[c]:
            best_complexes[idx] = b_c[c]  

    return best_complexes



def get_best_complex_distances(pdb_distances, best_complexes, mesure):
    mesures = dict()
    for pep_idx in pdb_distances.keys():
        if pep_idx not in best_complexes.keys():
            continue
        comp = best_complexes[pep_idx]
        pdb = comp.split('_')[2]
        #print(pep_idx, comp)
        #print(pdb_distances[pep_idx].keys())
        if comp not in pdb_distances[pep_idx].keys():
            continue
        mesures[pep_idx] = {'pdb':pdb, 'distance': pdb_distances[pep_idx][comp][mesure][0]}        
        
    return mesures


def get_complexes_lengths(complexes):
    lengths = dict()
    for p in complexes:
        l = 0
        for c in complexes[p]:
            norm = len(complexes[p][c])
            if norm>l:
                l=norm
            
        lengths[p] = l

    return lengths


# def big_best_complexes(complexes, clusters, lengths):
#     b_c = dict()
#     for p in complexes:
#         l = lengths[p]
#         best_complex = ''
#         best_accuracy = -1
#         for c in complexes[p]:
#             norm = len(complexes[p][c])
#             if (norm-l)>1 or (norm-l)<-1:
#                 #print('skipping')
#                 continue
        
#             trues = 0
#             for i in complexes[p][c]:
#                 if i:
#                     trues+=1
#             accuracy = trues/norm
#             if accuracy>best_accuracy:
#                 best_complex = c
#                 best_accuracy = accuracy
            
#         b_c[p] = best_complex

#     best_complexes = dict()
#     for c in clusters.keys():
#         if c not in b_c.keys():
#             continue
#         for idx in clusters[c]:
#             best_complexes[idx] = b_c[c]  

    
#     return best_complexes


##############################################################################
############################## RANDOM COMPLEXES #################################
##############################################################################

def get_random_complexes(complexes, clusters, seed=362718):
    ### change the seed for different sampling 
    random.seed(seed)
    selected_complexes = dict()
    n_complexes = dict()
    for p in complexes:
        ### store the number of complexes for each protein id
        n_complexes[p] = len(complexes[p])
        random_idx = random.randint(0, n_complexes[p]-1)

        for idx,c in enumerate(complexes[p]):
            if idx==random_idx:
                selected_complexes[p] = c

        # best_complex = ''
        # best_accuracy = -1
        # for c in complexes[p]:
        #     norm = len(complexes[p][c])
        #     trues = 0
        #     for i in complexes[p][c]:
        #         if i:
        #             trues+=1
        #     accuracy = trues/norm
        #     if accuracy>best_accuracy:
        #         best_complex = c
        #         best_accuracy = accuracy
            
        # b_c[p] = best_complex  

    #sort the complexes according to the peptide id
    random_complexes = dict()
    for c in clusters.keys():
        if c not in selected_complexes.keys():
            #print('Missing {}'.format(c))
            continue
        for idx in clusters[c]:
            random_complexes[idx] = selected_complexes[c]  


    with open('complexes_distribution.pickle', 'wb') as handle:
        pickle.dump(n_complexes, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
    return random_complexes

