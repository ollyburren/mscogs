#source /home/participant/Documents/mscogs/bin/activate
#exec(open("/home/participant/Documents/mscogs/mscogs-repo/py_code_nfg.py").read())

import rpy2.robjects as robjects
import numpy as np
import pandas as pd
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import re
from pybedtools import BedTool

#------------------------------------------------------------------------
##FUNCTIONS
#take a pandas df object and convert it to a BedTool object
def pnd_to_bt(df):
    return(BedTool(df.to_csv(index = False, header = False, sep = ' '), from_string = True))
#------------------------------------------------------------------------

#load SM2 RData object
robjects.r['load']("snpmod-fixed2.RData")
#get list names
names = robjects.r('names(SM2)')

#iterate through diseases and extract models`
models = {}
for a in names:
    p = "as.data.frame(attributes(SM2[['%s']])$models)" %a
    models[a] = robjects.r(p)

#example disease
z = models['GRAVES']
#which models have 'INS' in their snp names
np.where(list(map(lambda x: re.search('INS', x) != None, z['str'])))
#(array([139]),)

z1 = z['str'] #140

snps = list()
snps_flat = list()

for a in z1:
    snps0=list()
    for x in a.split('%'):
        x0 = x.split('.')
        snps0.append(x0)
        snps_flat.append(x0)
    snps.append(snps0)

#more than four elements in a snp name
np.where(list(map(lambda x: sum(list(map(lambda w: len(w) != 4, x))) > 0, snps)))
#(array([139]),)
#index of models with at least one snp without the rs ID
ind = np.where(list(map(lambda x: sum(list(map(lambda w: re.search('rs[0-9]*', w[0]) == None, x))) > 0, snps)))
#convert to vector
ind = np.concatenate(ind) + 1
#drop rows at ind
z_new = z.drop(ind)
z_new.reset_index()

#new snps, flat
snps_new = list()
for z in z_new['str']:
    for x in z.split('%'):
        x0 = x.split('.')
        snps_new.append(x0)

#snp info as a df (need to fix exceptions first)
snps_df = pd.DataFrame(snps_new, columns = ['rs_id', 'pos', 'a1', 'a2'])
#drop duplicate snps
snps_df = snps_df.drop_duplicates()
#adding chrm column
snps_df['chr'] = pd.Series([10] * len(snps_df), index = snps_df.index)

snps_str = snps_df.ix[:, ['chr', 'pos', 'pos', 'rs_id']].to_csv(index = False, header = False, sep = ' ')
snps_bed = BedTool(snps_str, from_string = True)

##genomic ranges a-la python
chic = pd.read_csv('merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab', sep = "\t")

bt = chic.ix[:,['baitChr', 'baitStart', 'baitEnd']].to_csv(index = False, header = False, sep = ' ')
chic_bed = BedTool(bt, from_string = True)

np.concatenate(np.where(chic['Erythroblasts'] > 5))

#which snps intersect any of the chic elements?
ww = snps_bed.intersect(chic_bed, u = True)
#data frame
snps_bed_df = pd.read_table(ww.fn, names = ['chr', 'start', 'stop', 'rs_id'])

snps_df['bait'] = pd.Series(np.isin(np.array(snps_df['rs_id']), np.array(snps_bed_df['rs_id'])), index = snps_df.index)

#chic = {pandas df with 4 columns in BedTools format + last column CHICAGO score}, snps = {BedTool object}
def boolSeries(chic, snps_bed):
    snps_df = pd.read_table(snps_bed.fn, names = ['chr', 'start', 'stop', 'rs_id'])
    chic_bed = chic.iloc[np.concatenate(np.where(chic.ix[:,[3]] > 5))]
    chic_bed = pnd_to_bt(chic_bed.ix[:, [0, 1, 2]])
    ww = snps_bed.intersect(chic_bed, u = True)
    snps_int_df = pd.read_table(ww.fn, names = ['chr', 'start', 'stop', 'rs_id'])
    return(pd.Series(np.isin(np.array(snps_df['rs_id']), np.array(snps_int_df['rs_id'])), index = snps_df.index))


ct = 'Erythroblasts'
chic_df = chic.ix[:, ['oeChr', 'oeStart', 'oeEnd', ct]]

ff = boolSeries(chic_df, snps_bed)




















if False:
    #junk
    #alternative for format
    [snps.append(x.split('.')) for x in z1.split('%')]
    #double map instead of double for loop
    pp = list(map(lambda z:
            list(map(lambda x: x.split('.'), z.split('%'))), z1
            )
        )

    w = pandas2ri.ri2py(z)
    u = robjects.r['SM2']
    #subscripting
    z.ix[ind]['lPP']


    def match(a, b):
        indx = [0] * len(a)
        for i in range(len(a)):
            if a[i] in b:
                 indx[i] = 1
        return(indx)

    match = lambda a, b: [ 1 if x in b else 0 for x in a ]
    match = lambda a, b: [ b.index(x)+1 if x in b else None for x in a ]

    chic_bed = chic.iloc[np.concatenate(np.where(chic['Erythroblasts'] > 5))]
    chic_bed = pnd_to_bt(chic_bed.ix[:, ['oeChr', 'oeStart', 'oeEnd']])
    ww = snps_bed.intersect(chic_bed, u = True)
    snps_bed_df = pd.read_table(ww.fn, names = ['chr', 'start', 'stop', 'rs_id'])
    return(pd.Series(np.isin(np.array(snps_df['rs_id']), np.array(snps_bed_df['rs_id'])), index = snps_df.index))






























##
