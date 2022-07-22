import numpy as np


######## function(s) for converting c++ values into python ones  ########

def read_cpp_into_floats(filename, is_print=True):
  """make dictionary of value: float for 
  doubles and ints in c++ file. Also make
  dictionary of notfloats for values that
  couldn't be converted"""

  constants = {}
  notfloats = {}
  with open(filename) as file:
      rlines=[]
      filelines = file.readlines()
      for line in filelines:
        ind1 = line.find("=")
        ind2 = line.find(";")
        if ind2 != -1 and ind1 != -1:
          line = line[:ind2+1]
          if "double" in line or "int" in line: 
            rlines.append(line)
      
      for line in rlines:
        if line[:13] == "const double ": x=13

        elif line[:10] == "const int ": x=10
        elif line[:4] == "int ": x=4
        elif line[:7] == "double ": x=7
        ind1 = line.find("=")
        ind2 = line.find(";")
        indx = line.find(" ", x) 
        symbol = line[x:indx]
        try: 
          float(line[ind1+2:ind2])
          constants[symbol] = float(line[ind1+2:ind2])
        except ValueError:
          notfloats[symbol] = line[ind1+2:ind2]
  
  if is_print:
    print_dict_statement(filename, "constants (CONSTS)", constants)
    print_dict_statement(filename, "not floats", notfloats)


  return constants, notfloats



def print_dict_statement(filename, dictname, dict):

  print("\n---- "+dictname+" from ", filename, "-----")
  for c in dict:
    print(c, "=", dict[c])
  print("---------------------------------------------\n")
  


def inits_dict(CONSTS, is_print=True):
  '''return INITS dictionary containing
    some specific key,values from 
    CONSTS dictionary'''

  INITS = {
    "iW"      : CONSTS["iW"],
    "DROPVOL" : CONSTS["DROPVOL"],
    "nsupers" : int(CONSTS["NSUPERS"]),
  }
  
  if is_print:
    print_dict_statement("CONSTS dict", "initial conditions (INITS)", INITS)

  return INITS



def mconsts_dict(CONSTS, is_print=True):
  '''return MCONSTS dictionary containing
    some derived key,values from values in
    CONSTS dictionary'''

  MCONSTS = {
    "RGAS_DRY"   : CONSTS["RGAS_UNIV"]/CONSTS["MR_DRY"],
    "RGAS_V"     : CONSTS["RGAS_UNIV"]/CONSTS["MR_WATER"],
    "CP0"        : CONSTS["CP_DRY"],
    "MR0"        : CONSTS["MR_DRY"],
    "Mr_ratio"   : CONSTS["MR_WATER"]/CONSTS["MR_DRY"],
  }
  MCONSTS["RHO0"]       = CONSTS["P0"]/(MCONSTS["RGAS_DRY"]*CONSTS["TEMP0"])

  if is_print:
    print_dict_statement("CONSTS dict", "derived constants (MCONSTS)", MCONSTS)

  return MCONSTS


