
def read_cpp_into_floats(filename):
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
  print("---- Constants read from ", filename, "-----")
  for c in constants:
    print(c, "=", constants[c])
  print("---------------------------------------------")
  print("---- Not floats read from ", filename, "-----")
  for st in notfloats:
    print(st, "=", notfloats[st])
  print("---------------------------------------------")

  return constants, notfloats

