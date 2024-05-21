######################################################################################################
#Code written by Valentin Vassilev-Galindo & Carmen Martínez-Alonso
"""
Copyright 2024 Valentin Vassilev-Galindo & Carmen Martínez-Alonso
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE
"""

#####################################################################################################

from mendeleev import element
from mendeleev.fetch import fetch_table
from ase import Atoms
import numpy as np

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

 PGtogeom = {'Fm-3m' : 'fcc', 'Im-3m' : 'bcc', 'P6_3/mmc' : 'hcp'}


 outere_zero_dict = {}
 for name in features.keys():
     if 'out_e' in name:
         letter = name[-1]
         outere_zero_dict['outere'+letter+'_zero'] = features[name]

 if formula_ok:
     formula = sample_info['system']
 else:
     formula = sample_info['system'].split('_')[0]+sample_info['system'].split('_')[1] if '_' in sample_info['system'] else sample_info['system']

 atoms = Atoms(formula)
 numbers = list(atoms.get_atomic_numbers())
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

         Done.append(n)
         i += 1
         if i == n_elements:
             break

 if sample_info['Type'] == 'Pure':
     eng_zero_dict['engB_zero'] = eng_zero_dict['engA_zero']
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

 return PSI






