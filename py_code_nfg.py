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

def bt_to_pnd(bt):
    return(pd.read_table(bt.fn, names = ['chr', 'start', 'stop', 'id']))

#chic = {pandas df with 4 columns in BedTools format + last column CHICAGO score}, snps = {BedTool object}
def boolSeries(chic, snps_bed):
    snps_df = bt_to_pnd(snps_bed)
    chic_bed = chic.iloc[np.concatenate(np.where(chic.ix[:,[3]] > 5))]
    chic_bed = pnd_to_bt(chic_bed.ix[:, [0, 1, 2]])
    ww = snps_bed.intersect(chic_bed, u = True)
    snps_int_df = bt_to_pnd(ww)
    return(np.isin(np.array(snps_df['id']), np.array(snps_int_df['id'])))

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


for ct in chic.columns[15:]:
    chic_df = chic.ix[:, ['oeChr', 'oeStart', 'oeEnd', ct]]
    snps_df[ct] = pd.Series(boolSeries(chic_df, snps_bed), index = snps_df.index)
    print(ct)

oe_df = chic.ix[:, ['oeChr', 'oeStart', 'oeEnd', 'ensg']]
bait_df = chic.ix[:, ['baitChr', 'baitStart', 'baitEnd', 'ensg']]

region_info = '10p-6030243-6169685'.split('-')
region_info[0] = re.sub("[pq]", '', region_info[0])
ri_bed = BedTool(' '.join(region_info), from_string = True)

oe_ol = pnd_to_bt(oe_df).intersect(ri_bed, u = True)
bait_ol = pnd_to_bt(bait_df).intersect(ri_bed, u = True)
oe_ol = bt_to_pnd(oe_ol)
bait_ol = bt_to_pnd(bait_ol)
uensg = pd.concat([oe_ol['id'], bait_ol['id']]).unique()
uensg_ind = np.isin(np.array(chic['ensg']), uensg)
chic_filt = chic[uensg_ind]


for ct in chic_filt.columns[15:]:
    chic_df = chic_filt.ix[:, ['oeChr', 'oeStart', 'oeEnd', ct]]
    snps_df[ct] = pd.Series(boolSeries(chic_df, snps_bed), index = snps_df.index)
    print(ct)

tmp = [None] * len(uensg)
for i in range(l):
    a = uensg[i]
    foo = snps_df.ix[:, [0, 1]]
    foo['ensg'] = pd.Series([a] * len(snps_df), index = foo.index)
    foo['chr'] = pd.Series([10] * len(snps_df), index = foo.index)
    foo = foo.ix[:, ['chr', 'pos', 'pos', 'rs_id']]
    tmp[i] = foo
    print(i)

snps_df_stack = pd.concat(tmp)

snps_str_stack = snps_df_stack.to_csv(index = False, header = False, sep = ' ')
snps_bed_stack = BedTool(snps_str_stack, from_string = True)


for ct in chic_filt.columns[15:]:
    chic_df = chic_filt.ix[:, ['oeChr', 'oeStart', 'oeEnd', ct]]
    snps_df_stack[ct] = pd.Series(boolSeries(chic_df, snps_bed_stack), index = snps_df_stack.index)
    print(ct)











##
