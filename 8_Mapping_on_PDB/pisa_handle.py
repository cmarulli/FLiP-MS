import os
import re


def get_pisa_info(pisa, pisa_file):
    with open(pisa_file, 'r') as f:
        pdb = ''
        n = 0
        for line in f:
            line = line.strip()
            info = re.split('<|>', line)
            try:
                info[1]
            except:
                #print(line)
                continue
                
            if info[1] == 'pdb_code':
                pdb = info[2]
                pisa[pdb] = dict()
            
            if info[1] == 'ser_no':
                n = info[2]
                pisa[pdb][n] = dict()
            
            if info[1] == 'diss_energy':
                pisa[pdb][n]['diss_energy'] = info[2]
            
            if info[1] == 'diss_area':
                pisa[pdb][n]['diss_area'] =  info[2]

            if info[1] == 'int_energy':
                pisa[pdb][n]['int_energy'] = info[2]
            
            if info[1] == 'bsa':
                pisa[pdb][n]['bsa'] = info[2]
            
            if info[1] == 'asa':
                pisa[pdb][n]['asa'] = info[2]

            if info[1] == 'entropy':
                pisa[pdb][n]['entropy'] = info[2]

    return



def build_pisa_dict(pisa_path):
    pisa = dict()
    n=0
    tot = len(os.listdir(pisa_path))
    for pdb in os.listdir(pisa_path):
        pisa_file = os.path.join(pisa_path,pdb)
        print('storing {}, {} of {}'.format(pdb,n,tot))
        get_pisa_info(pisa, pisa_file)
        n+=1
    
    return pisa



def get_higest_bsa_energies(pisa, cut=0):
    energies = list()
    low_energy_pdb = list()
    high_energy_pdb = list()

    for pdb in pisa:
        energy = None
        bsa = -1000000
        for i in pisa[pdb]:
            if float(pisa[pdb][i]['bsa'])>bsa:
                energy = pisa[pdb][i]['diss_energy']
        energies.append(float(energy))    
        
     
        if float(energy)<=cut:
            low_energy_pdb.append(pdb)
        else:
            high_energy_pdb.append(pdb)
    
    return {'energies':energies, 'low_energy_pdb':low_energy_pdb, 'high_energy_pdb':high_energy_pdb}


def split_bsa_energies(pisa, cut1=0, cut2=1):
    energies = list()
    low_energy_pdb = list()
    middle_energy_pdb = list()
    high_energy_pdb = list()

    for pdb in pisa:
        current_energies = list()
        energy = None
        bsa = -1000000
        for i in pisa[pdb]:
            # iterate over all pisa energy entries
            if float(pisa[pdb][i]['bsa'])>bsa:
                energy = float(pisa[pdb][i]['diss_energy'])
                # if the dissociation energy is not 0 store it
                if energy != 0:
                    current_energies.append(energy)
        
        # if there are non 0 dissociation energies store them and take the mean in case there are multiple onw for a single pdb
        if len(current_energies) > 0:
            energy = sum(current_energies)/len(current_energies)
            energies.append(energy)    
    
            # assign the current pdb as high, middle or low energy
            if float(energy)<=cut1:
                low_energy_pdb.append(pdb)
            elif float(energy)>cut1 and float(energy)<=cut2:
                middle_energy_pdb.append(pdb)
            else:
                high_energy_pdb.append(pdb)
    
    return {'energies':energies, 'low_energy_pdb':low_energy_pdb, 'middle_energy_pdb':middle_energy_pdb,'high_energy_pdb':high_energy_pdb}




def pisa_split(mesures, pisa_info):
    m_low = {k:mesures[k] for k in mesures.keys() if mesures[k]['pdb'] in pisa_info['low_energy_pdb']}
    m_middle = {k:mesures[k] for k in mesures.keys() if mesures[k]['pdb'] in pisa_info['middle_energy_pdb']}
    m_high = {k:mesures[k] for k in mesures.keys() if mesures[k]['pdb'] in pisa_info['high_energy_pdb']}
    
    return {'low_energies':m_low, 'middle_energies':m_middle, 'high_energies':m_high}
