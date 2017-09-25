#source /home/participant/Documents/mscogs/bin/activate

import rpy2.robjects as robjects
import numpy as np
import pandas as pa
from rpy2.robjects import pandas2ri
pandas2ri.activate()

robjects.r['load']("snpmod-fixed2.RData")
z = robjects.r('as.data.frame(attributes(SM2[[1]])$models)')








































##
