import os
import json
import pickle
import statistics
from ost import io
from ost import geom
from var3d import struct_anno
from var3d import pipeline


######SINGLE CHAIN FUNCTIONS#########

# def load_af_models(af_models_path, clusters):
#     '''returns a dictionary containing the single chain structure object for each uniprot id in clusters uniprot_id:structure'''
#     structures = dict()
#     for file in os.listdir(af_models_path):
#         ext = file.split('.')[1]
#         protein = file[3:9]
#         if protein in clusters.keys() and ext=='pdb':
#             file_name = os.path.join(af_models_path, file)
#             structures[protein] = io.LoadPDB(file_name)
    
#     print('loaded {} structures'.format(len(structures)))
    
#     return structures



# def mark_protein_surface(structure, treshold=25):
#     '''add the bool property "surface" to each residue of the structure'''
#     io.mol.alg.Accessibility(structure)
#     for res in structure.residues:
#         if res.GetFloatProp("asaRel")>treshold:
#             res.SetBoolProp('surface', True)
#         else:
#             res.SetBoolProp('surface', False)
#     return


# def get_surface_cleaveged_sites(peptides, structures):
#     '''Returns the list of K cleaved peptides (just the indexes) and a list that says which one of those has the cleavege site in the surface'''

#     peptides_K = list()
#     cleaved_surface = list()
#     complexes_cleaveges = dict()

#     for idx,pep in enumerate(peptides):

#         if pep.tryptics[0]=='Specific-C':
#             res_idx = pep.positions[0]
#             pep_code = pep.peptide[0]

#         elif pep.tryptics[0]=='Specific-N':
#             res_idx = pep.positions[0] + len(pep.peptide)-1
#             pep_code = pep.peptide[-1]

#         else:
#             continue
        
#         if pep.protein not in complexes_cleaveges.keys():
#             complexes_cleaveges[pep.protein] = list()

#         s = structures[pep.protein]
#         mark_protein_surface(s)

#         n_chains = 0
#         try:
#             for chain in s.chains:
#                 res = chain.FindResidue(res_idx)
#                 if res.one_letter_code!=pep_code:
#                     print('Mapping error in, ', res.one_letter_code, pep_code, pep.peptide, pep.tryptics[0])
#                     break
        
#                 peptides_K.append(idx)
#                 cleaved_surface.append(res.GetBoolProp('surface')) 
#                 complexes_cleaveges[pep.protein].append(res.GetBoolProp('surface'))
        
#                 n_chains+=1
#                 if n_chains>1:
#                     print('Error, more chains than expected')
#                     break
#         except:
#             print('Error in, ', idx)

#     return {'peptides_K':peptides_K, 'surface_results':cleaved_surface, 'p_coverage':complexes_cleaveges}




#########################################################################################################
########################################## SMR structures ###############################################
#########################################################################################################


def compute_minimum_distances(peptide, interface):
    min_distances = list()
    for res in peptide:
        res_list = geom.Vec3List()
        res_list.append(res)
        min_distances.append(geom.MinDistance(res_list, interface))
    return min_distances



def peptide_distances(idx_pep, peptides, uniprot_id, structure, annotations):
    n_chains = len(structure.assembly.Select("peptide=true").chains)
    pep_info = {'uniprot': uniprot_id, 'chains':n_chains}
    medians = list()
    averages = list()

    for segment_idx in range(structure.GetNumSegments()):
        interface_positions = geom.Vec3List()
        segment_annotation = annotations['struct_anno']['interface'][structure.hash][segment_idx]
        assert(len(segment_annotation)==len(structure.sequence))
        for r_idx in range(len(structure.sequence)):
            if segment_annotation[r_idx] == 1:
                # IMPORTANT: residue indexing is one based, thus the +1
                ca_pos = structure.GetCAPos(segment_idx, r_idx + 1)
                if ca_pos is not None:
                    interface_positions.append(ca_pos)
            
        #compute it for each segment?
        peptide_positions = geom.Vec3List()
        first = peptides[idx_pep].positions[0]-1
        last = first + len(peptides[idx_pep].peptide)
        for i in range(first,last):
            ca_pos = structure.GetCAPos(segment_idx, i + 1)
            if ca_pos is not None:
                peptide_positions.append(ca_pos)
            
        if len(peptide_positions)==0 or len(interface_positions)==0:
            medians.append(None)
            averages.append(None)

        else:
            min_dist = compute_minimum_distances(peptide_positions, interface_positions)
            medians.append(statistics.median(min_dist))
            averages.append(sum(min_dist)/len(min_dist))
        
        
    pep_info['medians'] = medians
    pep_info['averages'] = averages
    
    return pep_info



def compute_all_oligomeric_distances(peptides, clusters, jsons_path, data_type, output_path):
    pep_results = dict()
    n=0
    for protein in os.listdir(jsons_path):
        try:
            if len(protein.split('_')) > 1:
                continue
            n+=1
            print(protein, ' ', n)
        
            j_path = os.path.join(jsons_path, protein)
            with open(j_path, 'r') as fh:
                data = pipeline.DataContainer.FromJSON(json.load(fh))
    
            structures = list()
            for i in range(data.GetNumStructures()):
                s = data.GetStructure(i)
                if s.structure_identifier.split('_')[1]==data_type:
                    structures.append(s)
        
            if len(structures)>0:
                ann_file = protein.split('.')[0] + '_ann.json'
                j_ann_path = os.path.join(jsons_path, ann_file)
                with open(j_ann_path, 'r') as fh:
                    annotations = json.load(fh)
            else:
                print('No structures for {}'.format(protein))
                continue
        
            uniprot_id = protein.split('.')[0]
            for s in structures:
                #check if segments have interface
                n_chains = len(s.assembly.Select("peptide=true").chains)
                struct_id = s.structure_identifier
                for idx_pep in clusters[uniprot_id]:
                    if idx_pep not in pep_results.keys():
                        pep_results[idx_pep] = dict()
                    
                    if n_chains ==1:
                        pep_results[idx_pep][struct_id] = {'uniprot':uniprot_id, 'chains': n_chains, 'medians':[None], 'averages':[None]}
                    else:
                        pep_results[idx_pep][struct_id] = peptide_distances(idx_pep, peptides, uniprot_id, s, annotations)


            
        except:
            print('failed in: ', protein)
            
    with open(output_path, 'wb') as handle:
        pickle.dump(pep_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    return pep_results
