#source /home/participant/Documents/mscogs/bin/activate

import rpy2.robjects as robjects
import numpy as np
import pandas as pd
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import re
import BedTool from pybedtools

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

for z in z1:
    snps0=list()
    for x in z.split('%'):
        x0 = x.split('.')
        snps0.append(x0)
        snps_flat.append(x0)
    snps.append(snps0)

#more than four elements in a snp name
np.where(list(map(lambda x: sum(list(map(lambda w: len(w) != 4, x))) > 0, snps)))
#(array([139]),)
##--> need to decide how to fix this!

#snp info as a df (need to fix exceptions first)
#pa.DataFrame(snps, columns = 'rs_is', 'pos', 'a1', 'a2')

#snp positions
snps_pos = list(map(lambda x: x[1], snps_flat))


##genomic ranges a-la python
chic = pd.read_csv('merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab', sep="\t")
w = chic.ix[:,['ensg','baitChr','baitStart','baitEnd']]
bt = chic.ix[:,['baitChr','baitStart','baitEnd']].to_csv(index = False, header = False, sep = ' ')
foo = BedTool(bt, from_string=True)













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




































##
