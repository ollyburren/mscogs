#source /home/participant/Documents/mscogs/bin/activate
#exec(open("/home/participant/Documents/mscogs/mscogs-repo/py_code_nfg.py").read())

import rpy2.robjects as robjects
import numpy as np
import pandas as pd
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import re
import os
import os.path
from pybedtools import BedTool



## put the location of the data in an environment variable so code is more portable

data_dir = os.environ['MSCOGS_DATA']
model_r_obj_file = os.path.join(data_dir,'snpmod-fixed2.RData')
chic_file = os.path.join(data_dir,'merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab')

region_info = '10p-6030243-6169685'.split('-')

#------------------------------------------------------------------------
##FUNCTIONS
#take a pandas df object and convert it to a BedTool object
def pnd_to_bt(df):
    return(BedTool(df.to_csv(index = False, header = False, sep = ' '), from_string = True))

#take BedTools object and covert it to a pandas df object
def bt_to_pnd(bt):
    return(pd.read_table(bt.fn, names = ['chr', 'start', 'stop', 'id']))

#given a dataframe df and a list of columns lc, drop lc from df and return the resulting data frame
def drop_c(df, lc):
    cols = df.columns.tolist()
    blinx = ~np.isin(cols, lc)
    return(df.ix[:,blinx])

#------------------------------------------------------------------------
#load SM2 RData object
robjects.r['load'](model_r_obj_file)
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
chic = pd.read_csv(chic_file, sep = "\t")


oe_df = chic.ix[:, ['oeChr', 'oeStart', 'oeEnd', 'ensg']]
bait_df = chic.ix[:, ['baitChr', 'baitStart', 'baitEnd', 'ensg']]

#BedTool'ing region info
region_info[0] = re.sub("[pq]", '', region_info[0])
ri_bed = BedTool(' '.join(region_info), from_string = True)

#determining overlaps between bait and oe and the region of interest
oe_ol = pnd_to_bt(oe_df).intersect(ri_bed, u = True)
bait_ol = pnd_to_bt(bait_df).intersect(ri_bed, u = True)
oe_ol = bt_to_pnd(oe_ol)
bait_ol = bt_to_pnd(bait_ol)
#pulling together id's of genes overlapping the region; only keeping unique id's
uensg = pd.concat([oe_ol['id'], bait_ol['id']]).unique()
uensg_ind = np.isin(np.array(chic['ensg']), uensg)
#creating filtered chic object, only keeping genes overlapping region of interest
chic_filt = chic[uensg_ind]

#stacking snps_df's (1 for each unique gene in the region)
#creating list of snps_df

tmp = [None] * len(uensg)
for i in range(len(uensg)):
    a = uensg[i]
    indx = np.isin(chic_filt['ensg'], a)
    chic_tmp = chic_filt[indx]
    #looping throug cell types
    snps_df_stub = snps_df
    for ct in chic_filt.columns[15:]:
        chic_df = chic_tmp.ix[:, ['oeChr', 'oeStart', 'oeEnd', ct]]
        snps_df_stub[ct] = pd.Series(boolSeries(chic_df, snps_bed), index = snps_df_stub.index)
        print(ct)

    #foo = snps_df_stub.ix[:, [0, 1]]
    snps_df_stub = snps_df_stub.assign(ensg = pd.Series([a] * len(snps_df_stub), index = snps_df_stub.index))
    snps_df_stub = snps_df_stub.assign(chr = pd.Series([10] * len(snps_df_stub), index = snps_df_stub.index))
    #foo = foo.ix[:, ['chr', 'pos', 'pos', 'rs_id', 'ensg']]
    tmp[i] = snps_df_stub
    print(i)

#stacking all the snps_df in the list
snps_df_stack = pd.concat(tmp, ignore_index = True)

snps_bed_stack = pnd_to_bt(snps_df_stack.ix[:, [0, 1, 2, 3]])

#multistack case
for ct in chic_filt.columns[15:]:
    chic_df = chic_filt.ix[:, ['oeChr', 'oeStart', 'oeEnd', ct]]
    snps_df_stack[ct] = pd.Series(boolSeries(chic_df, snps_bed_stack), index = snps_df_stack.index)
    print(ct)



foo = cell_indicator(snps_df_stack, ['Neutrophils'])
ind = np.concatenate(np.where(foo == 1))

snps_df_stack.ix[ind,['rs_id', 'ensg']]

snps_df_stack.ix[ind]['ensg'].value_counts()
snps_df_stack['ensg'].value_counts()


#---------------------------------------------------------
#melting chic_filt dataframe
keys = chic_filt.columns[15:].tolist()
id_vars = chic_filt.columns[:15].tolist()
chic_melt = pd.melt(chic_filt, id_vars = id_vars, value_vars = keys)
#combine ensg id and varible (cell) name into one id
chic_melt['uid'] = pd.Series(["%s-%s" % t for t in zip(chic_melt['ensg'], chic_melt['variable'])], index = chic_melt.index)
chic_melt_bt = pnd_to_bt(chic_melt.ix[:, ['oeChr', 'oeStart', 'oeEnd', 'uid', 'value']])
res = snps_bed.intersect(chic_melt_bt, wb = True)
res = pd.read_table(ww.fn, names = ['snp_pos1', 'snp_pos2', 'snp_id', 'chr', 'ensg_start', 'ensg_stop', 'ensg_id', 'score'])
#filter according to score strength
res = res.loc[res['score'] > 5]
res.index = range(len(res))
##splitting ensg_id into a df with two columns
#uu = pd.DataFrame([x.split('-') for x in res['ensg_id']], columns = ['ensg', 'cell_type'])
##'cbinding' res and uu
#res = pd.concat([res, uu], axis = 1)
##dropping ensg_id as not needed anymore
#res = drop_c(res, 'ensg_id')
snps_by_ensg = {ensg: rows['snp_id'].unique() for ensg, rows in res.groupby('ensg_id')}



for key in snps_by_ensg:
    print "key: %s" % (key)


z1 = z['str'] #140

uu = pd.DataFrame([x.split('%') for x in z['str']])





tmp = [None] * len(uensg)
for i in range(len(uensg)):

a = uensg[i]
indx = np.isin(chic_filt['ensg'], a)
chic_tmp = chic_filt[indx]
#looping throug cell types
snps_df_stub = snps_df
for ct in chic_filt.columns[15:]:
    chic_df = chic_tmp.ix[:, ['oeChr', 'oeStart', 'oeEnd', ct]]
    snps_df_stub[ct] = pd.Series(boolSeries(chic_df, snps_bed), index = snps_df_stub.index)
    print(ct)

foo = snps_df_stub.ix[:, [0, 1]]
foo = foo.assign(ensg = pd.Series([a] * len(snps_df_stub), index = foo.index))
foo = foo.assign(chr = pd.Series([10] * len(snps_df_stub), index = foo.index))


foo.loc[:,'ensg'] = pd.Series([a] * len(snps_df_stub), index = foo.index)
foo.loc[:,'chr'] = pd.Series([10] * len(snps_df_stub), index = foo.index)
foo = foo.ix[:, ['chr', 'pos', 'pos', 'rs_id', 'ensg']]
tmp[i] = foo
print(i)










##
