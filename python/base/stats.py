import numpy
import rpy2.robjects 
R = rpy2.robjects.r

# =============================================================================
def correlate(x, y, method="pearson"):
  """
  performs a correlation between two vectors (assumed floats) and a given 
  correlation method. returns cor.coefficient and p-value.
  """
  xr = rpy2.robjects.FloatVector(x)
  yr = rpy2.robjects.FloatVector(y)
  res = R['cor.test'](xr, yr, method=method)
  #for i in range(len(res)):
  #  k = res.names[i]
  #  v = res[i]
  #  print i, "|", k, "=", v
  p = res.subset('p.value')[0][0]
  cor = res.subset('estimate')[0][0]
  return cor, p
  
# =============================================================================
def average(array):
  return numpy.average(array)
  
# =============================================================================
def median(array):
  return numpy.median(array)
  

# =============================================================================
def stdev(array):
  return numpy.std(array)
  
