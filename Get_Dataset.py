import os
import glob
import numpy as np
import sys
from ase.io import read
from ase.units import Rydberg
import pandas as pd
from mendeleev import element
from mendeleev.fetch import fetch_table, fetch_ionization_energies
from ase import Atoms
from PSI import get_psi
from GCN import get_gcn
from mp_api.client import MPRester
import pickle
import shutil
import json

df_elements = fetch_table('elements')
df_ionenergy  = fetch_ionization_energies()

Pure_surfs = {'111' : ['Ag','Al','Au', 'Cu','Ir','Ni','Pb','Pd', 'Pt','Rh'],
              '0001': ['Ca','Cd', 'Co', 'Hf','Mg','Os', 'Re', 'Ru', 'Tc','Zn','Zr'],
              '110' : ['Cr','Fe','Mo', 'Nb','Ta','W'],
              '100' : ['Ga','V'],
              '1010': ['Sc','Y'],
              '221' : ['Sn'],
              '1121': ['Ti'],
             }

Abnormal_syst = ['AB_bcc101', 'A3B_hcp0001']

_sites = {'A3B_fcc111': {'fccAAA': [(2./3., 1./6.)], 
                            'fccAAB': [(1./6.,1./6.), (2./3., 2./3.), (1./6., 2./3.)], 
                            'hcpAAA': [(1./3.,5./6.)], 
                            'hcpAAB': [(5./6., 5./6.), (5./6., 1./3.), (1./3., 1./3.)], 
                            'ontopA': [(0,0), (1./2., 0), (1./2., 1./2.)], 
                            'ontopB': [(1, 1./2.), (0, 1./2.)],
                           }, 
              'Pure': {'fcc_111' : { 'fcc': [(1./6.,1./6.),(2./3., 1./6.),(2./3., 2./3.), (1./6., 2./3.)],
                                     'hcp': [(1./3., 1./3.), (1./3.,5./6.),(5./6., 5./6.), (5./6., 1./3.)], 
                                     'ontop': [(0,0), (1./2., 0), (1./2., 1./2.), (0, 1./2.)],
                                     'bridge' : [(1./4., 0), (3./4., 0), (0, 1./4.), (1./4., 1./4.), (1./2., 1./.4), (3./4., 1./4.), (0, 3./4.),(1./4., 3./4.), (1./2., 3./4.), (3./4., 3./4.)],
                                   },
                       'bcc_110' : {'shortbridge': [(0, 1./4.), (1./2., 1./4.), (1./4., 1./4.), (3./4., 1./4.), (0, 3./4.), (1./4., 3./4.), (1./2., 3./4.), (3.4, 3./4.)], 
                                    'longbridge': [(1./4., 0), (3./4., 0), (1./4., 1./2.), (3./4., 1./2.)], 
                                    'hollow': [(1./6.,1./6.), (2./3., 1./6.), (1./3., 1./3.), (1./6., 2./3.), (1./3., 5./6), (5./6., 1./3.), (2./3., 2./3.), (5./6., 5./6.)],
                                    'ontop':[(0,0), (1./2., 0), (0, 1./2.), (1./2., 1./2.)],
                                   },
                       'bcc_100' : {'bridge':[(1./4., 0), (3./4.,0), (0, 1./4.), (1./2., 1./4.), (0,3./4.), (1./2., 3./4.)],
                                    'hollow':[(1./4., 1./4.), (3./4., 1./4.), (1./4., 3./4.), (3./4., 3./4.)],
                                    'ontop':[(0,0), (1./2., 0), (0, 1./2.), (1./2., 1./2.)],
                                   },
                       'hcp_0001': {'ontop': [(0, 0), (1./2., 0), (0, 1./2.), (1./2., 1./2.)], 
                                    'fcc': [(1./6.,1./6.),(2./3., 1./6.),(2./3., 2./3.), (1./6., 2./3.)],
                                    'hcp': [(1./3., 1./3.), (1./3.,5./6.),(5./6., 5./6.), (5./6., 1./3.)],
                                   },
                       'hcp_1010': {'ontop': [(0, 0), (1./2., 0), (0, 1./2.), (1./2., 1./2.)],
                                   },
                      },
              'AB_bcc101': {'ontopA': [(1./2.,0), (1./2, 1./2.)], #For all these coordinates to be correct, the atom in the position (1./2., 0) must be indeed A
                           'ontopB': [(0,1./4.), (0, 3./4.)], 
                           'shortbridge': [(1./4.,1./8.), (3./4., 1./8.), (1./4., 3./8.),(1./4.,5./8.), (1./4., 7./8.),(3./4., 3./8.), (3./4., 5./8.), (3./4, 7./8.)], 
                           'longbridgeB': [(0,1./2.), (0,0)], 
                           'longbridgeA': [(1./2.,1./4.), (1./2., 3./4.)], 
                           'threefoldAAB': [(1./3.,1./4.), (2./3., 1./4.), (1./3., 3./4.), (2./3., 3./4.)], 
                           'threefoldABB': [(1./6.,1./2.), (5./6., 1./2.),(1./6., 0), (5./6., 0)],
                          },
              'AB_fcc101': {'fccAAB': [(1./6., 2./3.), (2./3., 1./6.)], 
                           'fccABB' : [(1./6., 1./6.), (2./3., 2./3.)], 
                           'hcpAAB' : [(5./6., 1./3.), (1./3., 5./6.)], 
                           'hcpABB' : [(1./3., 1./3.), (5./6., 5./6.)],
                           'ontopA' : [(0, 0), (1./2., 1./2.)],
                           'ontopB' : [(1./2., 0), (0, 1./2.)],
                          },
              'A2B_hcp0001': {'fccAAB_1' : [(11./25., 11./20.), (11./25., 19./20.), (2./25., 11./20.)],
                              'fccAAB_2' : [(1./3., 1./6.), (2./25., 1./6.), (2./25, 1./3.)],
                              'hcpAAA_A' : [(0, 0)],
                              'hcpAAA_B' : [(1./3., 2./3.)],
                              'ontopA' : [(0.17, 0.83), (0.17, 0.34), (0.66, 0.83)],
                              'ontopB' : [(2./3., 1./3.)],
                             },
              'A3B_hcp0001': {'fccAAA' : [(0,0)],
                              'fccAAB' : [(1./2., 1./2.), (1./2., 0), (0, 1./2.)],
                              'hcpAAA' : [(1./3., 2./3.)],
                              'hcpAAB' : [(0.862, 0.724), (0.862, 0.138), (0.276, 0.138)], #Approximate values, given that the center of atoms might slightly differ
                              'ontopA' : [(0.138, 0.862), (0.724, 0.862), (0.138, 0.276)], #Approximate values, given that the center of atoms might slightly differ
                              'ontopB' : [(2./3., 1./3.)],
                             },

             }


_special_sites = {'AB_bcc101': {'ontopB': [(1./2.,0), (1./2, 1./2.)],
                           'ontopA': [(0,1./4.), (0, 3./4.)], 
                           'shortbridge': [(1./4.,1./8.), (3./4., 1./8.), (1./4., 3./8.),(1./4.,5./8.), (1./4., 7./8.),(3./4., 3./8.), (3./4., 5./8.), (3./4, 7./8.)], 
                           'longbridgeA': [(0,1./2.), (0,0)], 
                           'longbridgeB': [(1./2.,1./4.), (1./2., 3./4.)], 
                           'threefoldABB': [(1./3.,1./4.), (2./3., 1./4.), (1./3., 3./4.), (2./3., 3./4.)], 
                           'threefoldAAB': [(1./6.,1./2.), (5./6., 1./2.),(1./6., 0), (5./6., 0)],
                          },

                      'A3B_hcp0001': {'fccAAA' : [(0,0)],
                                        'fccAAB' : [(1./2., 1./2.), (1./2.,0), (0, 1./2.)],
                                        'hcpAAA' : [(2./3., 1./3.)],
                                        'hcpAAB' : [(0.138, 0.862), (0.724, 0.862), (0.138, 0.276)],
                                        'ontopA' : [(0.862, 0.724), (0.862, 0.138), (0.276, 0.138)],
                                        'ontopB' : [(1./3., 2./3.)],
                                       },
                      }

missing_eng_pauling = {#2 He n
                       #10 Ne
                       #18 Ar
                       '36' : 3.00, #Kr
                       '61' : 1.13, #Pm
                       '63' : 1.2,  #Eu
                       '65' : 1.1,  #Tb
                       '70' : 1.1,  #Yb
                       '86' : 2.2,  #Rn
                       '95' : 1.13, #Am
                       '96' : 1.28, #Cm
                       '97' : 1.3,  #Bk
                       '98' : 1.3,  #Cf
                       '99' : 1.3,  #Es
                       '100' : 1.3, #Fm
                       '101' : 1.3, #Md
                       '102' : 1.3, #No
                       '103' : 1.3, #Lr
                      }


E_mol = {'H2': -2.3332373616, #Energies are given in Ry
        'O2': -64.4520718752, #-64.3978649935,
        'H2O': -34.7828771047
        } #All these values were obtained with QuantumEspresso

Ry2eV = 13.605684958731

PGtogeom = {'Fm-3m' : 'fcc', 'Im-3m' : 'bcc', 'P6_3/mmc' : 'hcp'}


def get_features(sample_info, formula_ok = False, HinMat = False, OinMat = False):

    features = get_weighted_features(sample_info['system'])

    if formula_ok:
        formula = sample_info['system']
    else:
        if '_' in sample_info['system']:
            formula = sample_info['system'].split('_')[0]+sample_info['system'].split('_')[1]
        else:
            formula = sample_info['system']
    chem_symb = Atoms(formula).get_chemical_symbols()

    Done = []
    if not HinMat:
        Done.append('H')
    if not OinMat:
        Done.append('O')

    letters_dict = {'0':'A', '1':'B', '2':'C', '3':'D'}
    i = 0
    RE_dict = {'La' : 3.,
               'Ce' : 4.,
               'Pr' : 5.,
               'Nd' : 6.,
               'Pm' : 7.,
               'Sm' : 8.,
               'Eu' : 9.,
               'Gd' : 10.,
               'Tb' : 11.,
               'Dy' : 12.,
               'Ho' : 13.,
               'Er' : 14.,
               'Tm' : 15.,
               'Yb' : 16.,
               'Lu' : 17.,
               'Pa' : 5.,
               'Np' : 7.,
               'Th' : 4.,
               'U' : 6.,

              }
    letter_to_symb = {}
    for symb in chem_symb:
        if symb not in Done:
            group = df_elements['group_id'][df_elements['symbol'] == symb].values[0]
            dict_key = 'out_e'+letters_dict[str(i)]
            letter_to_symb[letters_dict[str(i)]] = symb
            if symb in RE_dict.keys():
                features[dict_key] = RE_dict[symb]
            else:
                features[dict_key] = group if group < 13 else group - 10

            i += 1
            Done.append(symb)


    if 'out_eB' not in features.keys():
        features['out_eB'] = features['out_eA']

    features['PSI'] = get_psi(sample_info, features)
    features['GCN'] = get_gcn(sample_info)

    for letter in letter_to_symb.keys():
        if 'chem_symb' not in sample_info.keys():
            features['magnetic_moment_'+letter] = 0
        else:
            symb_loc = np.where(np.array(sample_info['chem_symb']) == letter_to_symb[letter])[0]
            features['magnetic_moment_'+letter] = np.mean(np.array(sample_info['magnetic_moment'])[symb_loc])

    if 'magnetic_moment_B' not in features.keys():
        features['magnetic_moment_B'] = features['magnetic_moment_A']

    adsorbates = ['H','O']

    features['magnetic_moment_ads'] = 0
    if 'chem_symb' in sample_info.keys():
        n_ads = 0
        for atom in adsorbates:
            if atom in sample_info['chem_symb']:
                symb_loc = np.where(np.array(sample_info['chem_symb']) == atom)[0]
                features['magnetic_moment_ads'] += np.array(sample_info['magnetic_moment'])[symb_loc][0]
                n_ads += 1
     
        if n_ads > 1:
            features['magnetic_moment_ads'] /= n_ads

    return features

def get_weighted_features(system, formula_ok =False):
    properties = ['en_pauling', 'atomic_radius', 'ionenergy']
    names = ['WEN', 'WAR', 'WIE']

    if formula_ok:
        formula =system
    else:
        formula = system.split('_')[0]+system.split('_')[1] if '_' in system else system
    atoms = Atoms(formula)
    numbers = list(atoms.get_atomic_numbers()) # Use ASE Atoms to get atomic number list

    weighted_features = {}
    for name in names:
        weighted_features[name] = 0

    for n in set(numbers): # loop over unique atom numbers
        loc = df_elements['atomic_number'] == n
        data = df_elements[loc] # Get data table for element only
        weight = numbers.count(n)/len(numbers) # weight for calculating average property

        for i, prop in enumerate(weighted_features.keys()):
            feature_value = data[properties[i]].values[0] if properties[i] != 'ionenergy' else df_ionenergy['IE1'].values[n-1]
            if properties[i] == 'en_pauling':
                #print('ionenergy', feature_value)
                if np.isnan(feature_value):
                    feature_value = missing_eng_pauling[str(n)]
            weighted_features[prop] += weight * feature_value

    return weighted_features

def check_pseudopotential(INPUT):

    input_file = open(INPUT, 'r')
    
    i = 0
    is_valid =  True
    pseudo_tags = ['rrkjus', 'uspp']
    mol = read(INPUT)
    n_symb = np.unique(mol.get_chemical_symbols()).shape[0]
    j = 0
    for line in input_file:
        if j == n_symb:
            break
        if i == 1:
            j +=1
            pseudo = line.split()[2]
            if 'pbe' not in pseudo:
                print(('Warning: PBE functional was not used for %s...') % (line.split()[0]))
                is_valid = False
            n_tags = 0
            for tag in pseudo_tags:
                if tag in pseudo:
                    n_tags += 1
            if n_tags == 0:
                print(('Warning: Pseudopotential was not ultrasoft for %s...') % (line.split()[0]))
                is_valid =False

        if 'ATOMIC_SPECIES' in line:
            i = 1
    input_file.close()

    return is_valid


def get_slab_E(sample_info):

    is_valid = check_pseudopotential(sample_info['dir']+sample_info['input'])

    if is_valid:

        FILE = sample_info['dir']+sample_info['output']
        with open(FILE, 'r') as file:
            lines = file.readlines()

        if not '   JOB DONE.\n' in lines:
            print('\tJob did not end properly. Please check...')
            is_valid = False
            
            return is_valid, None
        elif '     history already reset at previous step: stopping\n' in lines:
            print('\tHistory reset message found. Please check...')
            for n_line, LINE in enumerate(lines):
                if 'bfgs converged' in LINE:
                    warning_idx = n_line #np.where(np.array(lines) == '     history already reset at previous step: stopping\n')[0][-1]
                    break
            E_err = float(lines[warning_idx-3].split()[-2]) if 'Ry' in lines[warning_idx-3] else float(lines[warning_idx-3].split()[-1])
            grad_err = float(lines[warning_idx-2].split()[-2]) if 'Ry' in lines[warning_idx-2] else float(lines[warning_idx-2].split()[-1]) 
            E_thres = float(lines[warning_idx+1].split()[3]) #6e-05
            grad_thres = float(lines[warning_idx+1].split()[7]) #1e-04
            if E_err > E_thres or grad_err > grad_thres:
                print('\tJob finished but it was not converged:')
                print('\t\tEnergy: error is %s (thres = %s)' % (str(E_err), str(E_thres)))
                print('\t\tForces: error is %s (thres = %s)' % (str(grad_err), str(grad_thres)))
                F = open('Wrong_convergence_slab.txt', 'a')
                F.write(('%s\t%s\t%s\t%s\t%s\n') % ( sample_info['system'], sample_info['Type'],'slab', sample_info['strain'], 'slab'))
                F.close()
                is_valid = False
                return is_valid, None


        else:
            E = float(list(filter(lambda x: 'Final energy' in x, lines))[0].split()[3])
     
            return is_valid, E
    else:

        return is_valid, None

def get_bulk_data(API_key, sample_info):

    mol_in = read(sample_info['dir']+sample_info['input'])
    symbols = np.unique(mol_in.get_chemical_symbols())
    not_valid = ['H', 'O']
    elements = [S for S in symbols if S not in not_valid]
    formula = sample_info['system'].split('_')[0]+sample_info['system'].split('_')[1] if '_' in sample_info['system'] else sample_info['system']

    if sample_info['Type'] == 'Pure':
        if sample_info['system'] == 'Cd':
            valid_structures = ['Hexagonal']
        else:
            valid_structures = ['Cubic', 'Hexagonal']
    else:
        valid_structures = ['Cubic'] if 'fcc' in sample_info['Type'] or 'bcc' in sample_info['Type'] else ['Hexagonal']

    with MPRester(str(API_key)) as mpr:
        docs = mpr.summary.search(elements=elements,formula=formula, fields=["material_id","symmetry","structure","energy_above_hull","volume"])
        hulls = np.array([doc.energy_above_hull for doc in docs])

        if docs[np.argmin(hulls)].symmetry.crystal_system in valid_structures:
            good_idx = np.argmin(hulls)
        else:
            for idx in np.argsort(hulls):
                if docs[idx].symmetry.crystal_system in valid_structures:
                    good_idx = idx
                    break

        Volume = docs[good_idx].volume
        Point_group = docs[good_idx].symmetry.symbol

    return Volume, Point_group

def calculate_Eads(E, E_slab, adsorbate):

    if adsorbate == 'OH':
        Eads = E - (E_mol['H2O'] - (0.5 * E_mol['H2'])) - E_slab
    else:
        Eads = E - (E_slab + (0.5 * E_mol[str(adsorbate)+'2']))

    return Eads * Ry2eV

def _canonical_site(out_pos):

    for i in range(out_pos.shape[0]):

        if out_pos[i] > 0.94 and out_pos[i] < 1.05:
            out_pos[i] = 0
        else:
            if out_pos[i] < 0:
                out_pos[i] += 1.
            if out_pos[i] > 1.:
                out_pos[i] -= 1.

    return out_pos

def get_site(out_pos, Type, thres = 0.08, point_group = None, system = None, ref_atom = None):

    if out_pos.shape[0] == 2:
        for i in range(out_pos.shape[0]):
            out_pos[i,:] = _canonical_site(out_pos[i,:])
        out_pos = np.squeeze(np.mean(out_pos,axis = 0)) #For OH we take the mean fractional position considering both the positions of O and H atoms
    else:
        out_pos = np.squeeze(out_pos)
        out_pos = _canonical_site(out_pos)

    if Type == 'Pure':
        geom = PGtogeom[point_group]
        for surf in Pure_surfs.keys():
            if system in Pure_surfs[surf]:
                pure_type = geom+'_'+surf
                break

    if Type in Abnormal_syst:
        
        ref_coordinates = [(0,0), (0, 1./2.)] if Type == 'AB_bcc101' else [(1./3., 2./3.)]
        ref_atom = _canonical_site(ref_atom)
        for i, coord in enumerate(ref_coordinates):
            upper_x = coord[0] + 0.05
            bottom_x = coord[0] - 0.05
            upper_y = coord[1] + 0.05
            bottom_y = coord[1] - 0.05
            if (ref_atom[0] > bottom_x and ref_atom[0] < upper_x) and (ref_atom[1] > bottom_y and ref_atom[1] < upper_y):
                special = False
                break
            else:
                special = True
        
        if special:
            print('The slab was constructed following a different procedure than typical. A and B are shifted...')
            print('Position of A or B in the second layer is:', ref_atom)
            valid_sites = _special_sites[Type]
        else:
            valid_sites = _sites[Type] 
    else:
        valid_sites = _sites[Type] if Type != 'Pure' else _sites[Type][pure_type]
        special = False
    
    was_found = False
    for i, key in enumerate(valid_sites.keys()):
        for j, coord in enumerate(valid_sites[key]):
            upper_x = coord[0] + thres
            bottom_x = coord[0] - thres
            upper_y = coord[1] + thres
            bottom_y = coord[1] - thres
            if (out_pos[0] > bottom_x and out_pos[0] < upper_x) and (out_pos[1] > bottom_y and out_pos[1] < upper_y):
                was_found = True
                return key, special

    if not was_found:
        return None, special

def get_data_from_pickle(sample_info, inout_info, E_slab):

    sample_data = {}
    sample_data['is_valid'] = inout_info['is_valid']

    E = inout_info['E']
    sample_data['Eads'] = calculate_Eads(E, E_slab, sample_info['adsorbate'])
        
    sample_data['frac_coord'] = inout_info['frac_coord']
    sample_data['cart_coord'] = inout_info['cart_coord']
    sample_data['Z'] = inout_info['Z']
    sample_data['chem_symb'] = inout_info['chem_symb']
    sample_data['cell'] = inout_info['cell']
    sample_info['chem_symb'] = inout_info['chem_symb']

    sample_data['magnetic_moment'] = inout_info['magnetic_moment']
    sample_info['magnetic_moment'] = inout_info['magnetic_moment']

    sample_data['initial_ad_pos'] = inout_info['initial_ad_pos']
    sample_data['final_ad_pos'] = inout_info['final_ad_pos']
    sample_info['final_ad_pos'] = inout_info['final_ad_pos']

    sample_data['init_frac_coord'] = inout_info['init_frac_coord']
    sample_data = get_additional_data(sample_info, sample_data)

    return sample_data


def get_data(sample_info, E_slab):

    sample_data = {}
    inout_info = {}
    mol_in = read(sample_info['dir']+sample_info['input'])

    sample_data['is_valid'] = check_pseudopotential(sample_info['dir']+sample_info['input'])
    inout_info['is_valid'] = check_pseudopotential(sample_info['dir']+sample_info['input'])

    if not sample_data['is_valid']:
        print('\tWrong pseudopotential used. Please check...')
        return sample_data, inout_info

    with open(sample_info['dir']+sample_info['output'], 'r') as file:
        lines = file.readlines()


    if not '   JOB DONE.\n' in lines:
        print('\tJob did not end properly. Please check...')
        F = open('No_JobEnd.txt', 'a')
        F.write(('%s\t%s\t%s\t%s\t%s\n') % ( sample_info['system'],sample_info['Type'], sample_info['adsorbate'], sample_info['strain'], sample_info['site']))
        F.close()
        sample_data['is_valid'] = False
        inout_info['is_valid'] = False
        return sample_data, inout_info

    if '     convergence NOT achieved after 200 iterations: stopping\n' in lines:
        print('\tJob did not converged. Please check...')
        F = open('No_converged.txt', 'a')
        F.write(('%s\t%s\t%s\t%s\t%s\n') % ( sample_info['system'],sample_info['Type'], sample_info['adsorbate'], sample_info['strain'], sample_info['site']))
        F.close()
        sample_data['is_valid'] = False
        inout_info['is_valid'] = False
        return sample_data, inout_info

    if '     history already reset at previous step: stopping\n' in lines:
        print('\tHistory reset message found. Please check...')
        for n_line, LINE in enumerate(lines):
            if 'bfgs converged' in LINE:
                warning_idx = n_line #np.where(np.array(lines) == '     history already reset at previous step: stopping\n')[0][-1]
                break
        E_err = float(lines[warning_idx-3].split()[-2]) if 'Ry' in lines[warning_idx-3] else float(lines[warning_idx-3].split()[-1])
        grad_err = float(lines[warning_idx-2].split()[-2]) if 'Ry' in lines[warning_idx-2] else float(lines[warning_idx-2].split()[-1]) 
        E_thres = float(lines[warning_idx+1].split()[3]) #6e-05
        grad_thres = float(lines[warning_idx+1].split()[7]) #1e-04
        if E_err > E_thres or grad_err > grad_thres:
            print('\tJob finished but it was not converged:')
            print('\t\tEnergy: error is %s (thres = %s)' % (str(E_err), str(E_thres)))
            print('\t\tForces: error is %s (thres = %s)' % (str(grad_err), str(grad_thres)))
            F = open('Wrong_convergence.txt', 'a')
            F.write(('%s\t%s\t%s\t%s\t%s\n') % ( sample_info['system'], sample_info['Type'],sample_info['adsorbate'], sample_info['strain'], sample_info['site']))
            F.close()
            sample_data['is_valid'] = False
            inout_info['is_valid'] = False
            return sample_data, inout_info

    mol_out = read(sample_info['dir']+sample_info['output'])

    if sample_data['is_valid'] and 'End final coordinates\n' in lines:
        E = float(list(filter(lambda x: 'Final energy' in x, lines))[0].split()[3])
        sample_data['Eads'] = calculate_Eads(E, E_slab, sample_info['adsorbate'])
        inout_info['E'] = E
        inout_info['final_coord_mess'] = True
    else:
        sample_data['Eads'] = None
        print('Eads not available. Please check...')
        F = open('No_FinalCooord_Message.txt', 'a')
        F.write(('%s\t%s\t%s\t%s\t%s\n') % ( sample_info['system'], sample_info['Type'],sample_info['adsorbate'], sample_info['strain'], sample_info['site']))
        F.close()
        sample_data['is_valid'] = False
        inout_info['is_valid'] = False

        return sample_data, inout_info

    ad_pos_in = np.array(mol_in.get_scaled_positions())[-2:,:] if sample_info['adsorbate'] == 'OH' else np.array(mol_in.get_scaled_positions())[-1,:][None,:]
    ad_pos_out = np.array(mol_out.get_scaled_positions())[-2:,:] if sample_info['adsorbate'] == 'OH' else np.array(mol_out.get_scaled_positions())[-1,:][None,:]

    sample_data['frac_coord'] = mol_out.get_scaled_positions()
    sample_data['cart_coord'] = mol_out.get_positions()
    sample_data['Z'] = mol_out.get_atomic_numbers()
    sample_data['chem_symb'] = mol_out.get_chemical_symbols()
    sample_info['chem_symb'] = mol_out.get_chemical_symbols()
    sample_data['cell'] = mol_out.cell[:]
    sample_data['init_frac_coord'] = mol_in.get_scaled_positions()

    inout_info['frac_coord'] = mol_out.get_scaled_positions()
    inout_info['cart_coord'] = mol_out.get_positions()
    inout_info['Z'] = mol_out.get_atomic_numbers()
    inout_info['chem_symb'] = mol_out.get_chemical_symbols()
    inout_info['cell'] = mol_out.cell[:]
    inout_info['init_frac_coord'] = mol_in.get_scaled_positions()

    try:
        sample_data['magnetic_moment'] = mol_out.get_magnetic_moments()
        sample_info['magnetic_moment'] = mol_out.get_magnetic_moments()
        inout_info['magnetic_moment'] = mol_out.get_magnetic_moments()
    except:
        sample_data['magnetic_moment'] = np.zeros(np.array(sample_info['chem_symb']).shape[0])
        sample_info['magnetic_moment'] = np.zeros(np.array(sample_info['chem_symb']).shape[0])
        inout_info['magnetic_moment'] = np.zeros(np.array(sample_info['chem_symb']).shape[0])

    sample_data['initial_ad_pos'] = np.copy(ad_pos_in)
    sample_data['final_ad_pos'] = np.copy(ad_pos_out)
    sample_info['final_ad_pos'] = np.copy(ad_pos_out)

    inout_info['initial_ad_pos'] = np.copy(ad_pos_in)
    inout_info['final_ad_pos'] = np.copy(ad_pos_out)

    sample_data = get_additional_data(sample_info, sample_data)

    return sample_data, inout_info

def get_additional_data(sample_info, sample_data):

    thres = 0.05 if sample_info['site'] != 'hollow' and 'threefold' not in sample_info['site'] else 0.09 
    site_name = sample_info['site'] if '-d' not in sample_info['site'] else sample_info['site'].split('-')[0]

    if sample_info['Type'] == 'Pure':
        initial_site, sample_info['special'] = get_site( sample_data['initial_ad_pos'], sample_info['Type'], thres = 0.08, point_group = sample_info['Point_group'], system = sample_info['system'])
    else:
        if sample_info['Type'] in Abnormal_syst:
            cut = -2 if sample_info['adsorbate'] == 'OH' else -1
            ref_positions = np.array(sample_data['init_frac_coord'])[:cut,:]
            z_lastlayer = np.argsort(ref_positions[:,2])[:4]
            ref_symbols = np.array(sample_data['chem_symb'])[z_lastlayer]
            ref_positions = ref_positions[z_lastlayer,:2]
            A_atom, B_atom = sample_info['system'].split('_')
            target_atom = A_atom if sample_info['Type'] == 'AB_bcc101' else B_atom

            for i, member in enumerate(ref_symbols):
                if member in target_atom:
                    ref_atom = ref_positions[i,:]
                    break
            initial_site, sample_info['special'] = get_site(sample_data['initial_ad_pos'], sample_info['Type'], thres = 0.08, ref_atom = ref_atom)
        else:
            ref_atom = None
            initial_site, sample_info['special'] = get_site(sample_data['initial_ad_pos'], sample_info['Type'], thres = 0.08)

    if initial_site != site_name and initial_site is not None:
        print('Warning: initial site does not correspond to input name. Please check...')
        F = open('Wrong_initialsitename.txt', 'a')
        F.write(('%s\t%s\t%s\t%s\t%s\t%s\n') % ( sample_info['system'],sample_info['Type'], sample_info['adsorbate'], sample_info['strain'], sample_info['site'], initial_site))
        F.close()
        wrong_site = True
    else:
        wrong_site = False
    
    diff_x = np.abs(sample_data['initial_ad_pos'][:,0] - sample_data['final_ad_pos'][:,0])
    diff_y =  np.abs(sample_data['initial_ad_pos'][:,1] - sample_data['final_ad_pos'][:,1])

    if sample_info['special']:
        token_y = thres + 0.05
        token_x = thres + 0.05
    else:
        token_y = np.mean(diff_y)
        token_x = np.mean(diff_x)

    
    if token_y > thres or token_x > thres:
    
        if (token_y > thres) and (token_y < 0.95 or token_y > 1.05):
            new_site, _ = get_site( sample_data['final_ad_pos'], sample_info['Type'], ref_atom = ref_atom) if sample_info['Type'] != 'Pure' else get_site( sample_data['final_ad_pos'], sample_info['Type'], point_group = sample_info['Point_group'], system = sample_info['system'])
    
            site2check = initial_site if wrong_site else site_name #sample_info['site']
    
            if new_site != site2check:
                if new_site is None:
                    new_site = site2check+'-d'
                print(('Site changed from %s to %s after optimization') % (site2check, new_site))
            else:
                if wrong_site:
                    new_site = initial_site
                    print(('Site changed from %s to %s after optimization') % (site_name, new_site))
    
            if new_site != sample_info['site']:
                sample_info['new_site'] = new_site
                sample_data['new_site'] = new_site
    
    
        elif (token_x > thres) and (token_x < 0.95 or token_x > 1.05):
            new_site, _ = get_site( sample_data['final_ad_pos'], sample_info['Type'], ref_atom = ref_atom ) if sample_info['Type'] != 'Pure' else get_site( sample_data['final_ad_pos'], sample_info['Type'], point_group = sample_info['Point_group'], system = sample_info['system'])
    
            site2check = initial_site if wrong_site else site_name #sample_info['site']
    
            if new_site != site2check:
                if new_site is None:
                    new_site = site2check+'-d'
                print(('Site changed from %s to %s after optimization') % (site2check, new_site))
            else:
                if wrong_site:
                    new_site = initial_site
                    print(('Site changed from %s to %s after optimization') % (site_name, new_site))
    
            if new_site != sample_info['site']:
                sample_info['new_site'] = new_site
                sample_data['new_site'] = new_site
    else:
        if wrong_site:
            sample_info['new_site'] = initial_site
            sample_data['new_site'] = initial_site

    features = get_features(sample_info) 

    for i, key in enumerate(features.keys()):
        if key not in sample_data.keys():
            sample_data[key] = features[key]

    for i, key in enumerate(sample_info.keys()):
        if key not in sample_data.keys():
            sample_data[key] = sample_info[key]

    Material = sample_info['system'].split('_')[0]+sample_info['system'].split('_')[1] if sample_info['Type'] != 'Pure' else sample_info['system']
    
    if sample_info['Type'] == 'Pure':
        geom = PGtogeom[sample_info['Point_group']]
    else:
        geom = sample_info['Type'].split('_')[-1]

    sample_data['ID'] = Material+'_'+geom+'_'+sample_info['strain']+'_'+sample_info['site']+'_'+sample_info['adsorbate']
    sample_data['geom'] = geom

    return sample_data

def check_file(sample_info, slab_data=None):

    FILE = open(sample_info['dir']+sample_info['input'], 'r')
    i = 0
    j = 0
    n_symb = 10
    if sample_info['Type'] == 'Pure':
        geom =''
    else:
        geom = sample_info['Type'].split('_')[-1]

    if slab_data is None:
        save_slab = True
        slab_data = {'H': 'H.pbe-rrkjus_psl.1.0.0.UPF',
                     'O': 'O.pbe-n-rrkjus_psl.1.0.0.UPF'
                    }
    else:
        save_slab = False
        symbs = []
        pseudos = []
    for line in FILE:
        if 'ntyp' in line:
            n_symb = int(line.split()[2])
        if j == n_symb:
            break
        if i == 1:
            j +=1
            symb = line.split()[0]
            pseudo = line.split()[2]
            if save_slab:
                slab_data[symb] = pseudo
            else:
                symbs.append(symb)
                pseudos.append(pseudo)

        if 'ATOMIC_SPECIES' in line:
            i = 1
        
    FILE.close()
    
    if not save_slab:
        error_file = open('pseudopotential_issue.txt', 'a')
        for i in range(len(symbs)):
            if pseudos[i] != slab_data[symbs[i]]:
                print(symbs[i], 'Used:', pseudos[i], 'Needed:', slab_data[symbs[i]])
                error_file.write(('%s\t%s\t%s\t%s\t%s\t%s\t%s\n') % (sample_info['system'], sample_info['adsorbate'], sample_info['strain'], sample_info['site'],geom, slab_data[symbs[i]], pseudos[i]))
        error_file.close()
        print()
    
    return slab_data


def select_better_site(n_site1, n_site2, ads_pos1, ads_pos2, Type, system, is_special, point_group = None):
    thres = 0.08

    site1 = n_site1.split('-')[0] if '-d' in n_site1 else n_site1
    site2 = n_site2.split('-')[0] if '-d' in n_site2 else n_site2

    pos_diff = np.mean(np.abs(ads_pos1[:,:2] - ads_pos2[:,:2]))
    if pos_diff > thres:
        print('The differences between adsorption positions are larger than %0.2f . Are you sure they are the same?')

    if Type == 'Pure':
        geom = PGtogeom[point_group]
        for surf in Pure_surfs.keys():
            if system in Pure_surfs[surf]:
                pure_type = geom+'_'+surf
                break

    if Type == 'Pure':
        pos_site1 = _sites[Type][pure_type][site1]
        pos_site2 = _sites[Type][pure_type][site2]
    else:
        pos_site1 = _special_sites[Type][site1] if is_special else _sites[Type][site1]
        pos_site2 = _special_sites[Type][site2] if is_special else _sites[Type][site2]

    mean_diff = 999
    if ads_pos1.shape[0] == 2:
        ads_pos1 = np.mean(ads_pos1, axis = 0)
        ads_pos2 = np.mean(ads_pos2, axis = 0)
    for pos in pos_site1:
        diff1 = np.mean(np.abs(pos - ads_pos1[:2]))
        diff2 = np.mean(np.abs(pos - ads_pos2[:2]))
        mean_tmp = (diff1 + diff2) / 2.
        if mean_tmp < mean_diff:
            mean_diff = mean_tmp
            site = n_site1
    for pos in pos_site2:
        diff1 = np.mean(np.abs(pos - ads_pos1[:2]))
        diff2 = np.mean(np.abs(pos - ads_pos2[:2]))
        mean_tmp = (diff1 + diff2) / 2.
        if mean_tmp < mean_diff:
            mean_diff = mean_tmp
            site = n_site2

    return site


if __name__ == '__main__':


    OS = sys.platform
    if OS == 'win32' or OS == 'cygwin':
        API_key_name = sys.argv[1] #Path to the txt file containing the API_key for using Materials Project within this script
        API_key=np.loadtxt(API_key_name,dtype='str')
        folder_lim = '\\'
    else:
        API_key=np.loadtxt(API_key_name,dtype='str')
        folder_lim ='/'

    if len(sys.argv) > 2: #In case you already have a pikle file with Dataset information you can load it an avoid collecting data you have
        overwrite = False
        with open(sys.argv[2], 'rb') as FILE:
            Available_dataset = pickle.load(FILE)
    else:
        overwrite = True
        Available_dataset ={'dir' :[]}


    valid_strains = ['Biaxial']
    valid_main = ['Bulk', 'No_strain', 'Strains']
    valid_adds = ['O','H','OH']
    tags_for_files = ['O','H','OH','slab']#,'Bulk']
    
    invalid_surfaces = ['Sc', 'Ti', 'Cr']
    invalid_ontop = ['Mg', 'Ca', 'Al','Ga', 'Sn', 'Ti']
    sites4ontop = ['ontopA','ontop', 'ontopB', 'fccAAA', 'hcpAAA', 'hcpAAA_A', 'hcpAAA_B']

    Bulks_Vol = {}
    Point_groups = {}
    Dataset='./'
    Token = True
    DIRS = sorted(glob.glob('*/'))
    print(DIRS)
    print()
    E_slabs = {}
    Dataset = []
    slab_pseudo = {}
    use_binaries = True
    while Token:
        DIRS_tmp = []
        Dataset_paths_tmp = []
        for i,d in enumerate(DIRS):
            sp = False
            valid_dirs = []
            dir_tag = d.split(folder_lim)[-2]
            if dir_tag == 'Strains':
                valid_dirs = valid_strains
                sp = True
            if dir_tag == 'No_strain' or '%' in dir_tag:
                system_tag = d.split(folder_lim)[-6] if '%' in dir_tag else d.split(folder_lim)[-3]
                valid_dirs = valid_adds + [system_tag+'slab']#+'/'
                sp = True
            if 'Bimetallic' in d and len(d.split(folder_lim)) == 5: 
                valid_dirs = valid_main
                sp = True
            if 'Pure_metals' in d and len(d.split(folder_lim)) == 3:
                valid_dirs = valid_main
                sp = True
            if 'Bimetallic' in d and 'No_strain' in d:
                if len(d.split(folder_lim)) == 8:
                    valid_dirs = []
                    sp = True
            if 'Bimetallic' in d and 'Strains' in d:
                if len(d.split(folder_lim)) == 11:
                    valid_dirs = []
                    sp = True
            if 'Pure_metals' in d and 'No_strain' in d:
                if len(d.split(folder_lim)) == 6:
                    valid_dirs = []
                    sp = True
            if 'Pure_metals' in d and 'Strains' in d:
                if len(d.split(folder_lim)) == 9: 
                    valid_dirs = []
                    sp = True
            if dir_tag == 'Bulk' or 'slab' in dir_tag:
                valid_dirs = []
                sp = True 
            if d != 'Gases'+folder_lim:
                DIRS_tmp2 = np.array([DIR if DIR.split(folder_lim)[-2] in valid_dirs else None for DIR in sorted(glob.glob(d+'*/'))]) if sp else [DIR for DIR in sorted(glob.glob(d+'*/'))]
                if sp:
                    DIRS_tmp2 = list(DIRS_tmp2[DIRS_tmp2 != None])
                DIRS_tmp += DIRS_tmp2
                if len(DIRS_tmp2) == 0:
                    inputs = sorted(glob.glob(d+'*.in'))
                    outputs = sorted(glob.glob(d+'*.out'))


                    if len(inputs) == 0:
                        print('WARNING: No inputs found in %s' % (d))
                    if len(outputs) == 0:
                        print('WARNING: No outputs found in %s' % (d))


                    if len(inputs) > 0 and len(outputs) > 0:
                        folder_tokens = np.array([True if string in d.split(folder_lim)[-2] else False for string in tags_for_files])
                        if folder_tokens.any():
                            system = d.split(folder_lim)[3] if 'Bimetallic' in d else d.split(folder_lim)[1]
                            Type = 'Pure' if 'Pure_metals' in d else d.split(folder_lim)[2]+'_'+d.split(folder_lim)[1].split('_')[0]+d.split(folder_lim)[1].split('_')[1]

                            if Type == 'Pure':
                                Binaries_path += folder_lim+'Pure_metals'
                            else:
                                Binaries_path += folder_lim+'Bimetallic'+folder_lim+d.split(folder_lim)[1]+folder_lim+Type.split('_')[0]

                            Binaries_path += folder_lim+system


                           
                            if not folder_tokens[-1]:# and not folder_tokens[-1]:
                                if 'OH' in inputs[0]:
                                    adsorbate = 'OH'
                                elif 'O.' in inputs[0] or 'O_' in inputs[0]:
                                    adsorbate = 'O'
                                else:
                                    adsorbate = 'H'
                               
                                site = inputs[0].split(folder_lim)[-2].split('.')[0].split(system)[-1].split(adsorbate)[0]
                            else:
                                adsorbate = 'Slab'
                                site = ''

                            if 'Strains' in d:
                                strain = d.split(folder_lim)[-3].split('%')[0] if folder_tokens[-1] else d.split(folder_lim)[-4].split('%')[0]
                                Binaries_path += folder_lim+'Strains'+folder_lim+'Biaxial'
                                if int(strain) > 8:
                                    continue
                                if 'Compression' in d:
                                    Binaries_path += folder_lim+'Compression'+folder_lim+strain+'%'
                                    strain = '-'+strain
                                else:
                                    Binaries_path += folder_lim+'Tension'+folder_lim+strain+'%'
                            else:
                                strain = '0'
                                Binaries_path += folder_lim+'No_strain'

                            print(system, Type, strain, adsorbate, site)

                            in_file = inputs[0].split(folder_lim)[-1]
                            out_file = outputs[-1].split(folder_lim)[-1]

                            sample_info = {'input': in_file,
                                    'output': out_file,
                                    'dir': d,
                                    'system': system,
                                    'adsorbate': adsorbate,
                                    'strain': strain,
                                    'site': site,
                                    'Type': Type,
                                    }
                            

                            if system in invalid_surfaces or system in invalid_ontop:
                                print('Sample not valid...\n')
                                F = open('Invalid_surfaces.txt', 'a')
                                F.write(('%s\t%s\t%s\t%s\t%s\n') % ( sample_info['system'], sample_info['Type'], sample_info['adsorbate'], sample_info['strain'], sample_info['site']))
                                F.close()
                                continue

                            if adsorbate == 'Slab':
                                is_valid, E_slab = get_slab_E(sample_info)

                                if is_valid:
                                    slab_key = system+'_'+Type+'_'+strain
                                    if slab_key not in E_slabs.keys():
                                        E_slabs[slab_key] = E_slab
                                    else:
                                        print('Warning: Possible repeated slab...')
                            else:
                                if not overwrite and sample_info['dir'] in Available_dataset['dir']:
                                    print('Sample already in dataset...')
                                    sample_data = {}
                                    sample_idx = np.where(np.array(Available_dataset['dir']).ravel() == sample_info['dir'])[0]
                                    for key in Available_dataset.keys():
                                        sample_data[key] = Available_dataset[key][sample_idx[0]]
                                    Dataset.append(sample_data)
                                else:
                                    print('New sample for dataset...')
                                    sample_data = {}
                                    if system not in Bulks_Vol.keys():
                                        Bulks_Vol[system], Point_groups[system] = get_bulk_data(API_key, sample_info)
                                    
                                    sample_info['Volume'] = Bulks_Vol[system]
                                    sample_info['Point_group'] = Point_groups[system]
                                    slab_key = system+'_'+Type+'_'+strain
                                    if slab_key not in E_slabs.keys():# or system not in Bulks_Vol.keys():
                                        print('Slab energy is not available. Cannot compute Eads. Please check...')
                                        F = open('No_slabs_av.txt', 'a')
                                        F.write(('%s\t%s\t%s\t%s\t%s\n') % ( sample_info['system'], sample_info['Type'], sample_info['adsorbate'], sample_info['strain'], sample_info['site']))
                                        F.close()
                                        sample_data['is_valid'] = False
                                    else:
                                        full_sample_name = system+site+adsorbate

                                        sample_data, inout_info = get_data(sample_info, E_slabs[slab_key])
     
                                    if sample_data['is_valid']:
                                        Dataset.append(sample_data)

                            print()

        DIRS = DIRS_tmp
        if len(DIRS) == 0:
            Token = False

    MP_data = {'Volumes' : Bulks_Vol, 'Point_groups' : Point_groups}
    
    print(('Dataset size: %i') % (len(Dataset)))
    None_sites = []
    Changed_sites = []
    available_sites = {}
    Eads_dict = {}
    for i in range(len(Dataset)):
        if 'new_site' in Dataset[i].keys():
            if Dataset[i]['new_site'] is None:
                None_sites.append(Dataset[i])
            else:
                Changed_sites.append(Dataset[i])
        else:
            system = Dataset[i]['system']
            strain = Dataset[i]['strain']
            adsorbate = Dataset[i]['adsorbate']
            site = Dataset[i]['site']
            tag_sites = system+'_'+strain+'_'+adsorbate
            tag_Eads = system+'_'+strain+'_'+adsorbate+'_'+site
            Eads = Dataset[i]['Eads']
            if tag_sites not in available_sites.keys():
                available_sites[tag_sites] = []
            available_sites[tag_sites].append(site)
            if tag_sites not in Eads_dict.keys():
                Eads_dict[tag_sites] = {}
            Eads_dict[tag_sites][site] = Eads

    del_IDs = []
    change_sites = {'Changed' : [],
                    'None': [],
                }
    d_sites =[]
    print('Size of None sites:', len(None_sites))
    for i in range(len(Changed_sites)):
        tag = Changed_sites[i]['system']+'_'+Changed_sites[i]['strain']+'_'+Changed_sites[i]['adsorbate']
        if tag in available_sites.keys():
            if Changed_sites[i]['new_site'] in available_sites[tag]:
                del_IDs.append(Changed_sites[i]['ID'])
                Rep_file = open('Ad_changed_invalid.txt', 'a')
                Rep_file.write(('%s\t%s\t%s\t%s\t%s\t%s\n') % (Changed_sites[i]['system'], Changed_sites[i]['Type'], Changed_sites[i]['adsorbate'], Changed_sites[i]['strain'], Changed_sites[i]['site'], Changed_sites[i]['new_site']))
                Rep_file.close()
            else:
                available_sites[tag].append(Changed_sites[i]['new_site'])
                if tag not in Eads_dict.keys():
                    Eads_dict[tag] = {}
                Eads_dict[tag][Changed_sites[i]['new_site']] = Changed_sites[i]['Eads']
                if '-d' in Changed_sites[i]['new_site']:
                    d_sites.append(Changed_sites[i])
                else:
                    change_sites['Changed'].append(Changed_sites[i]['ID'])
                    Rep_file = open('Ad_changed_valid.txt', 'a')
                    Rep_file.write(('%s\t%s\t%s\t%s\t%s\t%s\n') % (Changed_sites[i]['system'], Changed_sites[i]['Type'], Changed_sites[i]['adsorbate'], Changed_sites[i]['strain'], Changed_sites[i]['site'], Changed_sites[i]['new_site']))
                    Rep_file.close()
        else:
            available_sites[tag] = [Changed_sites[i]['new_site']]
            if tag not in Eads_dict.keys():
                Eads_dict[tag] = {}
            Eads_dict[tag][Changed_sites[i]['new_site']] = Changed_sites[i]['Eads']
            if '-d' in Changed_sites[i]['new_site']:
                d_sites.append(Changed_sites[i])
            else:
                change_sites['Changed'].append(Changed_sites[i]['ID'])
                Rep_file = open('Ad_changed_valid.txt', 'a')
                Rep_file.write(('%s\t%s\t%s\t%s\t%s\t%s\n') % (Changed_sites[i]['system'], Changed_sites[i]['Type'], Changed_sites[i]['adsorbate'], Changed_sites[i]['strain'], Changed_sites[i]['site'], Changed_sites[i]['new_site']))
                Rep_file.close()

    thres = 0.014 # This is equal to 0.001 Ry
    print('Number of displaced sites:', len(d_sites))
    d_rm_tmp = []
    for i in range(len(d_sites)):
        if i in d_rm_tmp:
            continue
        tag_i = d_sites[i]['system']+'_'+d_sites[i]['strain']+'_'+d_sites[i]['adsorbate']
        for j in range(i+1, len(d_sites)):
            if j in d_rm_tmp:
                continue
            tag_j = d_sites[j]['system']+'_'+d_sites[j]['strain']+'_'+d_sites[j]['adsorbate']
            if tag_i == tag_j:
                diff = np.abs(Eads_dict[tag_i][d_sites[i]['new_site']] - Eads_dict[tag_j][d_sites[j]['new_site']])
                if diff < thres:
                    good_site = select_better_site(d_sites[i]['new_site'], d_sites[j]['new_site'], d_sites[i]['final_ad_pos'], d_sites[j]['final_ad_pos'], d_sites[i]['Type'], d_sites[i]['system'], d_sites[i]['special'], point_group = d_sites[i]['Point_group'])
                    idx2rm = i if good_site == d_sites[i]['new_site'] else j
                    other_idx = j if idx2rm == i else i
                    del_IDs.append(d_sites[idx2rm]['ID'])
                    d_rm_tmp.append(idx2rm)
                    Rep_file = open('Ad_changed_invalid.txt', 'a')
                    Rep_file.write(('%s\t%s\t%s\t%s\t%s\t%s\n') % (d_sites[idx2rm]['system'], d_sites[idx2rm]['Type'], d_sites[idx2rm]['adsorbate'], d_sites[idx2rm]['strain'], d_sites[idx2rm]['site'], d_sites[other_idx]['new_site']))
                    Rep_file.close()

        is_i_valid = True if i not in d_rm_tmp else False
        if is_i_valid:
            for tag_k in Eads_dict.keys():
                if tag_i == tag_k:
                    d_sites_removed = [d_sites[j]['new_site'] for j in d_rm_tmp]
                    for site_k in Eads_dict[tag_k].keys():
                        if site_k in d_sites_removed:
                            continue
                        if site_k != d_sites[i]['new_site']:
                            diff = np.abs(Eads_dict[tag_i][d_sites[i]['new_site']] - Eads_dict[tag_k][site_k])
                            if diff < thres:
                                is_i_valid = False
                                del_IDs.append(d_sites[i]['ID'])
                                d_rm_tmp.append(i)
                                Rep_file = open('Ad_changed_invalid.txt', 'a')
                                Rep_file.write(('%s\t%s\t%s\t%s\t%s\t%s\n') % (d_sites[i]['system'], d_sites[i]['Type'], d_sites[i]['adsorbate'], d_sites[i]['strain'], d_sites[i]['site'], site_k))
                                Rep_file.close()
                                break

        if is_i_valid:
            change_sites['Changed'].append(d_sites[i]['ID'])
            Rep_file = open('Ad_changed_valid.txt', 'a')
            Rep_file.write(('%s\t%s\t%s\t%s\t%s\t%s\n') % (d_sites[i]['system'], d_sites[i]['Type'], d_sites[i]['adsorbate'], d_sites[i]['strain'], d_sites[i]['site'], d_sites[i]['new_site']))
            Rep_file.close()

    for i in range(len(None_sites)):
        tag =  None_sites[i]['system']+'_'+None_sites[i]['strain']+'_'+None_sites[i]['adsorbate']
        Token = True
        if tag in Eads_dict.keys():
            for site in Eads_dict[tag].keys():
                diff = np.abs(None_sites[i]['Eads'] - Eads_dict[tag][site])
                if diff < thres:
                    Token = False
                    del_IDs.append(None_sites[i]['ID'])
                    Rep_file = open('Ad_changed_invalid.txt', 'a')
                    Rep_file.write(('%s\t%s\t%s\t%s\t%s\t%s\n') % (None_sites[i]['system'], None_sites[i]['Type'], None_sites[i]['adsorbate'], None_sites[i]['strain'], None_sites[i]['site'], site))
                    Rep_file.close()
                    break
            if Token:
                Eads_dict[tag][None_sites[i]['site']+'-d'] = None_sites[i]['Eads']
                change_sites['None'].append(None_sites[i]['ID'])
                Rep_file = open('Ad_changed_valid.txt', 'a')
                Rep_file.write(('%s\t%s\t%s\t%s\t%s\t%s\n') % (None_sites[i]['system'], None_sites[i]['Type'], None_sites[i]['adsorbate'], None_sites[i]['strain'], None_sites[i]['site'], None_sites[i]['site']+'-d'))
                Rep_file.close()
        else:
            Eads_dict[tag] = {None_sites[i]['site']+'-d' : None_sites[i]['Eads']}
            change_sites['None'].append(None_sites[i]['ID'])
            Rep_file = open('Ad_changed_valid.txt', 'a')
            Rep_file.write(('%s\t%s\t%s\t%s\t%s\t%s\n') % (None_sites[i]['system'], None_sites[i]['Type'], None_sites[i]['adsorbate'], None_sites[i]['strain'], None_sites[i]['site'], None_sites[i]['site']+'-d'))
            Rep_file.close()
    
    print((' %i samples will be deleted from dataset...') % (len(del_IDs)))
    Final_dataset = {}
    H_dataset = {}
    O_dataset = {}
    OH_dataset = {}
    O_OH_dataset = {}
    all_keys = [key if key != 'new_site' else None for key in Dataset[0].keys()]
    all_keys = np.array(all_keys)[np.array(all_keys) != None]
    for key in all_keys:
        Final_dataset[key] = []
        H_dataset[key] = []
        O_dataset[key] = []
        OH_dataset[key] = []
        O_OH_dataset[key] = []


    keys2change = ['site', 'ID', 'dir','input', 'output']

    for i in range(len(Dataset)):
        if Dataset[i]['ID'] not in del_IDs:
            adsorbate = Dataset[i]['adsorbate']
            for key in Final_dataset.keys():
                if key in keys2change:
                    if Dataset[i]['ID'] in change_sites['Changed']:
                        valid_site = Dataset[i]['new_site']
                    elif Dataset[i]['ID'] in change_sites['None']:
                        valid_site = Dataset[i]['site']+'-d'
                    else:
                        valid_site = Dataset[i]['site']
                    if key == 'ID': 
                        Material = Dataset[i]['system'].split('_')[0]+Dataset[i]['system'].split('_')[1] if Dataset[i]['Type'] != 'Pure' else Dataset[i]['system']     
                        new_ID = Material+'_'+Dataset[i]['geom']+'_'+Dataset[i]['strain']+'_'+valid_site+'_'+Dataset[i]['adsorbate']
                        to_append = new_ID

                    if key == 'dir' or key == 'input' or key == 'output':
                        if valid_site != Dataset[i]['site']:
                            new_name = Dataset[i]['system']+valid_site+Dataset[i]['adsorbate']
                            if key == 'dir':
                                dir_parts = Dataset[i][key].split(folder_lim)[:-2]
                                new_dir = ''
                                for dir_part in dir_parts:
                                    new_dir += dir_part+folder_lim
                                new_dir += new_name+folder_lim
                                to_append = new_dir
                            else:
                                to_append = new_name+'.in' if key == 'input' else new_name+'.out'
                        else:
                            to_append = Dataset[i][key]
                    if key == 'site':
                        to_append = valid_site

                    Final_dataset[key].append(to_append)
                    if adsorbate == 'H':
                        H_dataset[key].append(to_append)
                    else:
                        O_OH_dataset[key].append(to_append)
                        if adsorbate == 'O':
                            O_dataset[key].append(to_append)
                        else:
                            OH_dataset[key].append(to_append)
                else:
                    Final_dataset[key].append(Dataset[i][key])
                    if adsorbate == 'H':
                        H_dataset[key].append(Dataset[i][key])
                    else:
                        O_OH_dataset[key].append(Dataset[i][key])
                        if adsorbate == 'O':
                            O_dataset[key].append(Dataset[i][key])
                        else:
                            OH_dataset[key].append(Dataset[i][key])
        else:
            print(('Removing sample %s from dataset...') % (Dataset[i]['ID']))

    print()
    print(('Final Dataset size: %i') % (len(Final_dataset['ID'])))
    print(('H Dataset size: %i') % (len(H_dataset['ID'])))
    print(('O Dataset size: %i') % (len(O_dataset['ID'])))
    print(('OH dataset size: %i') % (len(OH_dataset['ID'])))
    print(('O-OH dataset size: %i') % (len(O_OH_dataset['ID'])))
   
    with open('Full_dataset_H-O-OH.pickle','wb') as f:
        pickle.dump(Final_dataset, f)

    Final_dataset_json = {}
    for key in Final_dataset.keys():
        Final_dataset_json[key] = [Final_dataset[key][i].tolist() for i in range(len(Final_dataset[key]))] if isinstance(Final_dataset[key][0],np.ndarray) else Final_dataset[key]
    with open( "Full_dataset_H-O-OH.json" , "w" ) as f:
        json.dump( Final_dataset_json , f )
