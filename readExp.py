import sys
from pyexcel_ods import get_data  

exp = str(sys.argv[1])
print('Codename of experiment:', exp)

data = get_data('../exps/experiments.ods')

test = exp in (item for sublist in data['Sheet1'] for item in sublist)
if test==False:
    print('ERROR: No experiment "{}"'.format(exp))
    quit()

exp, minradVoid, rmin, rmax, sec, fxa = next(row for row in data['Sheet1'] if row[0] == exp) 

print('minradVoid =',minradVoid)
print('rmin =',rmin)
print('rmax =',rmax)
print('sec =',sec)
print('fxa =',fxa)

