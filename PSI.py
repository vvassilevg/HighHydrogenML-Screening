from mendeleev import element
from mendeleev.fetch import fetch_table
from ase import Atoms
import numpy as np

#YOU NEED TO INCLUDE THE OUTER ELECTRONS OF A AND B, THE ELECTRONEGATIVITY OF A AND B, THE ESTEQUIOMETRY, GEOMETRY AND ADSORPTION SITE

#OUTER ELECTRONS
#outereA_zero= 12
#outereB_zero= 5

#ELECTRONEGATIVITY (Pauli)
#engA_zero= 1.65
#engB_zero= 1.6

#estequiom= 'A2B'   #options: A3B, AB, A (pure), A2B
#geom= 'hcp'       #options: fcc, bcc, hcp
#site= 'fccAAB_2'    #options: fcc, hcp, ontop, fccAAA, fccAAB, fccABB, hcpAAA, hcpAAB, hcpABB ontopA, ontopB, shortbridge, longbridgeAA, longbridgeBB, threefoldAAB, threefoldABB, hcpAAA_A, hcpAAA_B, fccAAB_1, fccAAB_2.


df_elements = fetch_table('elements')

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



def get_psi(sample_info, features, formula_ok = False):

 #DO NOT CHANGE ANYTHING BELOW THIS LINE
 #--------------------------------------------------------------------------------------------

 #proportions of A and B depending on stoichiometry and adsorption site
 #n es átomos en total (A+B)

 PGtogeom = {'Fm-3m' : 'fcc', 'Im-3m' : 'bcc', 'P6_3/mmc' : 'hcp'}


 outere_zero_dict = {}
 for name in features.keys():
     if 'out_e' in name:
         letter = name[-1]
         outere_zero_dict['outere'+letter+'_zero'] = features[name]

 print('outer e Dict for PSI:')
 print(outere_zero_dict)
 print()

 #outereA_zero = features['out_eA']
 #outereB_zero = features['out_eB']
 

 if formula_ok:
     formula = sample_info['system']
 else:
     formula = sample_info['system'].split('_')[0]+sample_info['system'].split('_')[1] if '_' in sample_info['system'] else sample_info['system']

 atoms = Atoms(formula)
 numbers = list(atoms.get_atomic_numbers())
 #print(formula, numbers)
 #print()
 n_elements = np.unique(numbers).shape[0]

 letters_dict = {0 : 'A', 1 : 'B', 2 : 'C', 3 : 'D'}

 eng_zero_dict = {}
 Done = []
 i = 0
 for n in numbers:
     if n not in Done:
         letter = letters_dict[i]
         eng_tmp = df_elements[df_elements['atomic_number'] == n]['en_pauling'].values[0]
         if np.isnan(eng_tmp):
             eng_tmp = missing_eng_pauling[str(n)]

         eng_zero_dict['eng'+letter+'_zero'] = eng_tmp

         #print(i, n)
         #if i == 0:
         #    engA_zero = df_elements[df_elements['atomic_number'] == n]['en_pauling'].values[0]
         #    if np.isnan(engA_zero):
         #        engA_zero = missing_eng_pauling[str(n)]
         #if i == 1:
         #    engB_zero = df_elements[df_elements['atomic_number'] == n]['en_pauling'].values[0]
         #    if np.isnan(engB_zero):
         #        engB_zero = missing_eng_pauling[str(n)]
         Done.append(n)
         i += 1
         if i == n_elements:
             break

 print('Eng dict for PSI:')
 print(eng_zero_dict)
 print()

 if sample_info['Type'] == 'Pure':
     #engB_zero = engA_zero
     eng_zero_dict['engB_zero'] = eng_zero_dict['engA_zero']
     #eng_zero_dict['engC_zero'] = eng_zero_dict['engA_zero'] #NEEDED IF DATASET INCLUDES TRIMETALIC
     estequiom = 'A'
     geom = PGtogeom[sample_info['Point_group']]
 else:
     estequiom = sample_info['Type'].split('_')[0]
     if 'hcp' in sample_info['Type']:
         geom = 'hcp'
     if 'fcc' in sample_info['Type']:
         geom = 'fcc'
     if 'bcc' in sample_info['Type']:
         geom = 'bcc'

     #NEED TO ADD A ENG_C FOR BIMETALIC SYSTEMS

 if 'new_site' in sample_info.keys():
     if sample_info['new_site'] is None:
         site = sample_info['site']
     else:
         site = sample_info['new_site'] if '-d' not in sample_info['new_site'] else sample_info['new_site'].split('-')[0]
 else:
     site = sample_info['site'] if '-d' not in sample_info['site'] else sample_info['site'].split('-')[0]

 
 Proportions_dict = {}

 #ABC hcp 
 if (estequiom=='ABC' and geom=='hcp' and site=='fccABC'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] =range(0)
    Proportions_dict['C']=range(0)
    n=3
 elif (estequiom=='ABC'and geom=='hcp' and site=='bridgeBC'):
    Proportions_dict['A']= 0
    Proportions_dict['B']= range(0)
    Proportions_dict['C']=range(0)
    n=2
 elif (estequiom=='ABC'and geom=='hcp' and site=='ontopA'):
    Proportions_dict['A']= range(0)
    Proportions_dict['B']=range(2)
    Proportions_dict['C']=range(2)
    n=7
 elif (estequiom=='ABC'and geom=='hcp' and site=='ontopB'):
    Proportions_dict['A']= 0
    Proportions_dict['B']= range(0)
    Proportions_dict['C']=range(2)
    n=4
 elif (estequiom=='ABC'and geom=='hcp' and site=='ontopC'):
    Proportions_dict['A']= 0
    Proportions_dict['B']= range(2)
    Proportions_dict['C']= range(0)
    n=4


 #A3BC fcc 
 if (estequiom=='A3BC' and geom=='fcc' and site=='fccAAB'):
    Proportions_dict['A'] = range(1)
    Proportions_dict['B'] =range(0)
    Proportions_dict['C'] =range(0)
    n=4
 elif (estequiom=='A3BC'and geom=='fcc' and site=='hcpAAA'):
    Proportions_dict['A'] = range(2)
    Proportions_dict['B'] = 0
    Proportions_dict['C'] =range(0)
    n=4
 elif (estequiom=='A3BC'and geom=='fcc' and site=='hcpAAB'):
    Proportions_dict['A'] = range(2)
    Proportions_dict['B'] =range(0)
    Proportions_dict['C'] =range(0)
    n=5
 elif (estequiom=='A3BC'and geom=='fcc' and site=='ontopA'):
    Proportions_dict['A'] = range(6)
    Proportions_dict['B'] = range(2)
    Proportions_dict['C'] =range(1)
    n=12
 elif (estequiom=='A3BC'and geom=='fcc' and site=='ontopB'):
    Proportions_dict['A'] = range(2)
    Proportions_dict['B'] = range(0)
    Proportions_dict['C'] = 0
    n=4
 elif (estequiom=='A3BC'and geom=='fcc' and site=='ontopC'):
    Proportions_dict['A'] = 0
    Proportions_dict['B'] = 0
    Proportions_dict['C'] =range(0)
    n=1

 #A3B fcc 
 if (estequiom=='A3B' and geom=='fcc' and site=='fccAAA'):
    Proportions_dict['A'] = range(5)
    Proportions_dict['B'] =0
    n=6
 elif (estequiom=='A3B'and geom=='fcc' and site=='fccAAB'):
    Proportions_dict['A'] = range(3)
    Proportions_dict['B'] = range(1)
    n=6
 elif (estequiom=='A3B'and geom=='fcc' and site=='hcpAAA'):
    Proportions_dict['A'] = range(2)
    Proportions_dict['B'] =0
    n=3
 elif (estequiom=='A3B'and geom=='fcc' and site=='hcpAAB'):
    Proportions_dict['A'] = range(1)
    Proportions_dict['B']= range(0)
    n=3
 elif (estequiom=='A3B'and geom=='fcc' and site=='ontopA'):
    Proportions_dict['A'] = range(6)
    Proportions_dict['B'] = range(2)
    n=10
 elif (estequiom=='A3B'and geom=='fcc' and site=='ontopB'):
    Proportions_dict['A'] = range(8)
    Proportions_dict['B'] = range(0)
    n=10
	
 #A3B hcp 
 if (estequiom=='A3B' and geom=='hcp' and site=='fccAAA'):
    Proportions_dict['A'] = range(5)
    Proportions_dict['B'] =0
    n=6
 elif (estequiom=='A3B'and geom=='hcp' and site=='fccAAB'):
    Proportions_dict['A'] = range(3)
    Proportions_dict['B'] = range(1)
    n=6
 elif (estequiom=='A3B'and geom=='hcp' and site=='hcpAAA'):
    Proportions_dict['A'] = range(2)
    Proportions_dict['B'] =0
    n=3
 elif (estequiom=='A3B'and geom=='hcp' and site=='hcpAAB'):
    Proportions_dict['A'] = range(1)
    Proportions_dict['B'] = range(0)
    n=3
 elif (estequiom=='A3B'and geom=='hcp' and site=='ontopA'):
    Proportions_dict['A'] = range(6)
    Proportions_dict['B'] = range(2)
    n=10
 elif (estequiom=='A3B'and geom=='hcp' and site=='ontopB'):
    Proportions_dict['A'] = range(8)
    Proportions_dict['B'] = range(0)
    n=10	
		

 #PURE (A)	
 elif (estequiom=='A'and geom=='bcc' and site=='fcc'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1
 elif (estequiom=='A' and geom=='bcc'and site=='hcp'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1
 elif (estequiom=='A' and geom=='bcc'and site=='ontop'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1	
 elif (estequiom=='A' and geom=='bcc'and site=='hollow'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1
 elif (estequiom=='A' and geom=='bcc'and site=='shortbridge'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1
 elif (estequiom=='A' and geom=='bcc'and site=='longbridge'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1
 elif (estequiom=='A' and geom=='bcc'and site=='bridge'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1

 elif (estequiom=='A'and geom=='fcc' and site=='fcc'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1
 elif (estequiom=='A' and geom=='fcc'and site=='hcp'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1
 elif (estequiom=='A' and geom=='fcc'and site=='ontop'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1	
 elif (estequiom=='A' and geom=='fcc'and site=='hollow'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1
 elif (estequiom=='A' and geom=='fcc'and site=='bridge'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1
	
 elif (estequiom=='A'and geom=='hcp' and site=='fcc'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1
 elif (estequiom=='A' and geom=='hcp'and site=='hcp'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1
 elif (estequiom=='A' and geom=='hcp'and site=='ontop'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1	
 elif (estequiom=='A' and geom=='hcp'and site=='hollow'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = 0
    n=1

 #AB bcc	 	
 elif (estequiom=='AB' and geom=='bcc' and site=='shortbridge'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = range(0)
    n=2
 elif (estequiom=='AB' and geom=='bcc' and site=='ontopA'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = range(5)
    n=7
 elif (estequiom=='AB' and geom=='bcc' and site=='ontopB'):
    Proportions_dict['A'] = range(5)
    Proportions_dict['B'] = range(0)
    n=7
 elif (estequiom=='AB' and geom=='bcc' and site=='longbridgeA'):
    Proportions_dict['A'] = range(1)
    Proportions_dict['B'] = range(0)
    n=3
 elif (estequiom=='AB' and geom=='bcc' and site=='longbridgeB'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = range(1)
    n=3
 elif (estequiom=='AB' and geom=='bcc' and site=='threefoldAAB'):
    Proportions_dict['A'] = range(1)
    Proportions_dict['B'] = range(0)
    n=3
 elif (estequiom=='AB' and geom=='bcc' and site=='threefoldABB'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = range(1)
    n=3

 #AB fcc
 elif (estequiom=='AB' and geom=='fcc' and site=='fccAAB'):
    Proportions_dict['A'] = range(3)
    Proportions_dict['B'] = range(1)
    n=6		
 elif (estequiom=='AB' and geom=='fcc' and site=='fccABB'):
    Proportions_dict['A'] = range(1)
    Proportions_dict['B'] = range(3)
    n=6		
 elif (estequiom=='AB' and geom=='fcc' and site=='hcpAAB'):
    Proportions_dict['A'] = range(1)
    Proportions_dict['B'] = range(0)
    n=3	
 elif (estequiom=='AB' and geom=='fcc' and site=='hcpABB'):
    Proportions_dict['A'] = range(0)
    Proportions_dict['B'] = range(1)
    n=3
 elif (estequiom=='AB' and geom=='fcc' and site=='ontopA'):
    Proportions_dict['A'] = range(3)
    Proportions_dict['B'] = range(5)
    n=10	
 elif (estequiom=='AB' and geom=='fcc' and site=='ontopB'):
    Proportions_dict['A'] = range(5)
    Proportions_dict['B'] = range(3)
    n=10
	
 #A2B hcp
 elif (estequiom=='A2B' and geom=='hcp' and site=='fccAAB_1'):
    Proportions_dict['A'] = range(1)
    Proportions_dict['B'] = range(1)
    n=4		
 elif (estequiom=='A2B' and geom=='hcp' and site=='fccAAB_2'):
    Proportions_dict['A'] = range(2)
    Proportions_dict['B'] = range(1)
    n=5		
 elif (estequiom=='A2B' and geom=='hcp' and site=='hcpAAA_A'):
    Proportions_dict['A'] = range(3)
    Proportions_dict['B'] = 0
    n=4	
 elif (estequiom=='A2B' and geom=='hcp' and site=='hcpAAA_B'):
    Proportions_dict['A'] = range(2)
    Proportions_dict['B'] = range(0)
    n=4
 elif (estequiom=='A2B' and geom=='hcp' and site=='ontopA'):
    Proportions_dict['A'] = range(4)
    Proportions_dict['B'] = range(1)
    n=7	
 elif (estequiom=='A2B' and geom=='hcp' and site=='ontopB'):
    Proportions_dict['A'] = range(5)
    Proportions_dict['B'] = range(0)
    n=7
	

 #outereA= outereA_zero
 #outereB= outereB_zero
 #engA= engA_zero
 #engB= engB_zero

 outere_dict = {}
 eng_dict = {}
 for letter in Proportions_dict.keys():
     if Proportions_dict[letter] == 0:
         outere_tmp = 1
         eng_tmp = 1
     else:
         outere_tmp = outere_zero_dict['outere'+letter+'_zero']
         eng_tmp = eng_zero_dict['eng'+letter+'_zero']
         for i in Proportions_dict[letter]:
             outere_tmp *= outere_zero_dict['outere'+letter+'_zero']
             eng_tmp *= eng_zero_dict['eng'+letter+'_zero']

     outere_dict[letter] =  outere_tmp
     eng_dict[letter] = eng_tmp

 p1_e = 1
 p2_e = 1
 p1_eng = 1
 for letter in outere_dict.keys():
     p1_e *= outere_dict[letter]
     p2_e *= outere_dict[letter]
     p1_eng *= eng_dict[letter]

 p1_e = p1_e ** (1. / n)
 p2_e = p2_e ** (1. / n)
 
 geometric_mean_outere_squared = p1_e * p2_e

 geometric_mean_eng = p1_eng ** (1. / n)

 PSI = geometric_mean_outere_squared / geometric_mean_eng

 print('PSI:', PSI)

 """
 #new variables depending on quantity
 if Proportions_dict['A'] == 0:
    outereA = 1
 else:
    for i in Proportions_dict['A']:
        outereA = outereA_zero * outereA
        #print (outereA)

 if Proportions_dict['B'] == 0:
    outereB = 1
 else:
    for i in Proportions_dict['B']:
        outereB = outereB_zero * outereB
        #print (outereB)

 geometric_mean_outere_squared_old = ((outereA*outereB)**(1/n)) * ((outereA*outereB)**(1/n))
 #print (geometric_mean_outere_squared)
 
 
 #new variables depending on quantity
 if Proportions_dict['B'] == 0:
    engA = 1
 else:
    for i in A:
        engA = engA_zero * engA
        #print (engA)

 if B == 0:
    engB = 1
 else:
    for i in B:
        engB = engB_zero * engB
        #print (engB)
 
 geometric_mean_eng_old = (engA*engB)**(1/n)
 #print (geometric_mean_eng)



 #-------------------------------------------------------------------------------------------

 PSI_old = geometric_mean_outere_squared_old/geometric_mean_eng_old
 #print ('PSI is equal to', (PSI))
 """

 return PSI






