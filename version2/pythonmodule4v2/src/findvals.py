import numpy as np




def index_and_value(value, array, array2=[]):
  ''' function for finding index where 
    a given array is closest 
    to a specific value. Returns either 
    index and value of given array at index,
    and (optionally) value of a different array,
    called array2, at same index '''

  index = np.argmin(abs(array-value))
  
  if isinstance(array2, np.ndarray):
    return index, array[index], array2[index] 

  else:
    return index, array[index]
  