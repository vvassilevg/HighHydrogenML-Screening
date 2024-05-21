######################################################################################################
#Code written by Valentin Vassilev-Galindo & Carmen Martínez-Alonso
"""
Copyright 2024 Valentin Vassilev-Galindo & Carmen Martínez-Alonso
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE
"""

#####################################################################################################



def get_gcn(sample_info):

 PGtogeom = {'Fm-3m' : 'fcc', 'Im-3m' : 'bcc', 'P6_3/mmc' : 'hcp'}
 bcc_surf = {'Cr': '110', 'Fe' :'110', 'Mo' : '110', 'Nb' : '110', 'Ta' : '110', 'V' : '100', 'W' :'110'}

 Type = sample_info['Type']
 estequiom = Type.split('_')[0]

 if estequiom == 'Pure':
     if sample_info['Point_group'] in PGtogeom.keys():
         geom = PGtogeom[sample_info['Point_group']]
     else:
         print(('Point group not in keys: %s - Please check...') % (sample_info['Point_group']))
 else:
     if 'hcp' in Type:
         geom = 'hcp'
     if 'fcc' in Type:
         geom = 'fcc' 
     if 'bcc' in Type:
         geom = 'bcc'

 if Type != 'Pure':
     surface = Type.split('_')[-1].split(geom)[-1] # '111'   #options: 111, 101, 0001, 110, 100
 else:
     if geom == 'fcc':
         surface = '111'
     if geom == 'hcp':
         surface = '0001'
     if geom == 'bcc':
         if sample_info['system'] in bcc_surf.keys():
             surface = bcc_surf[sample_info['system']]
         else:
             print('Element not available. Setting to 110...')
             surface = '110'

 if 'new_site' in sample_info.keys():
     if sample_info['new_site'] is None:
         site = sample_info['site']
     else:
         site = sample_info['new_site'] if '-d' not in sample_info['new_site'] else sample_info['new_site'].split('-')[0]
 else:
     site = sample_info['site'] if '-d' not in sample_info['site'] else sample_info['site'].split('-')[0]
 
 if 'bridge' in site:
     if 'long' in site:
         site = 'longbridge'
     elif 'short' in site:
         site ='shortbridge'
     else:
         site = 'bridge'

 if 'ontop' in site:
     site = 'ontop'

 if 'fcc' in site:
     site = 'fcc'

 if 'hcp' in site:
     site = 'hcp'

 if 'threefold' in site:
     site = 'threefold'


 #FCC 111 = FCC 101 = HCP 0001
 if (surface=='111' and geom=='fcc' and site=='fcc'):
    GCN=5.25
 elif (surface=='111' and geom=='fcc' and site=='hcp'):
    GCN=3.25
 elif (surface=='111' and geom=='fcc' and site=='ontop'):
    GCN=0.75
 elif (surface=='111' and geom=='fcc' and site=='bridge'):
    GCN = 0
	
 elif (surface=='101' and geom=='fcc' and site=='fcc'):
    GCN=5.25
 elif (surface=='101' and geom=='fcc' and site=='hcp'):
    GCN=3.25
 elif (surface=='101' and geom=='fcc' and site=='ontop'):
    GCN=0.75
	
 elif (surface=='0001' and geom=='hcp' and site=='fcc'):
    GCN=5.25
 elif (surface=='0001' and geom=='hcp' and site=='hcp'):
    GCN=3.25
 elif (surface=='0001' and geom=='hcp' and site=='ontop'):
    GCN=0.75
 elif (surface=='0001' and geom == 'hcp' and site == 'bridge'):
     GCN = 1.5

 #BCC 101 = BCC 110
 elif (surface=='101' and geom=='bcc' and site=='shortbridge'):
    GCN=1.5
 elif (surface=='101' and geom=='bcc' and site=='ontop'):
    GCN=0.75
 elif (surface=='101' and geom=='bcc' and site=='longbridge'):
    GCN=2.5
 elif (surface=='101' and geom=='bcc' and site=='threefold'):
    GCN=3.25
 elif (surface=='101' and geom=='bcc' and site=='hollow'):
    GCN=3.25
	
 elif (surface=='110' and geom=='bcc' and site=='shortbridge'):
    GCN=1.5
 elif (surface=='110' and geom=='bcc' and site=='ontop'):
    GCN=0.75
 elif (surface=='110' and geom=='bcc' and site=='longbridge'):
    GCN=2.5
 elif (surface=='110' and geom=='bcc' and site=='threefold'):
    GCN=3.25	
 elif (surface=='110' and geom=='bcc' and site=='hollow'):
    GCN=3.25	


 #BCC 100
 elif (surface=='100' and geom=='bcc' and site=='hollow'):
    GCN=3
 elif (surface=='100' and geom=='bcc' and site=='ontop'):
    GCN=0.5
 elif (surface=='100' and geom=='bcc' and site=='bridge'):
    GCN=1
	
 return GCN





