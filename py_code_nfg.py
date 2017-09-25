#source /home/participant/Documents/mscogs/mscogs/bin/activate

import rpy2.robjects as robjects
import numpy as np
import pandas as pa
from rpy2.robjects import pandas2ri
pandas2ri.activate()

robjects.r['load']("snpmod-fixed2.RData")
u = robjects.r['SM2']

for i in range(len(u))
z = robjects.r('as.data.frame(attributes(SM2[[1]])$models)')
#dataframe
w = pandas2ri.ri2py(z)





























##
