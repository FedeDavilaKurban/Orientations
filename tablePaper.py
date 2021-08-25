#%%
import numpy as np
from tabulate import tabulate
from astropy.io import ascii
from astropy.table import Table


# %%
minradV, maxradV = 7.0, 0.0


a = []
r = []
s = []
spin = []
mass = []
vrad = []
s5 = []
for vfilter in ['-','hi','lo']:

    for sfilter in ['-','hi','lo']:

        for sec in [0,123,1,3,456,4,6,14,36]:

            zeta = []
            zeta_std = []

            for vtype in ['a','r','s']:

                if vfilter!='-':
                    filename = '../data/eta/eta_vfilter{}_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
                        .format(vfilter,minradV,maxradV,sec,vtype)
                else: 
                    filename = '../data/eta/eta_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
                        .format(minradV,maxradV,sec,vtype)
                etaTable = ascii.read(filename)

                eta0 = 1./(np.sqrt(2)-1)
                eta_ran_std = np.sqrt(28.1421/etaTable['N'][0])

                zeta.append( (etaTable['eta'][0]-eta0)/eta_ran_std )
                zeta_std.append(etaTable['eta_std'][0]/eta_ran_std)
                #if zeta[-1]<3: break

            if sec==0: 
                spin.append('-')
                mass.append('-')
            if sec==1: 
                spin.append('H')
                mass.append('L')
            if sec==3: 
                spin.append('H')
                mass.append('H')
            if sec==4: 
                spin.append('L')
                mass.append('L')
            if sec==6: 
                spin.append('L')
                mass.append('H')
            if sec==14: 
                spin.append('-')
                mass.append('L')
            if sec==36: 
                spin.append('-')
                mass.append('H')
            if sec==123: 
                spin.append('H')
                mass.append('-')
            if sec==456: 
                spin.append('L')
                mass.append('-')
            if vfilter == 'hi': vrad.append('H')
            if vfilter == 'lo': vrad.append('L')
            if vfilter == '-': vrad.append('-')

            if sfilter == 'hi': s5.append('H')
            if sfilter == 'lo': s5.append('L')
            if sfilter == '-': s5.append('-')

            a.append( '{:.1f}+-{:.1f}'.format(zeta[0],zeta_std[0]) )
            r.append( '{:.1f}+-{:.1f}'.format(zeta[1],zeta_std[1]) )
            s.append( '{:.1f}+-{:.1f}'.format(zeta[2],zeta_std[2]) )


info = {'Spin': [s for s in spin], 'Mass': [m for m in mass], r'$\Sigma_5$': [s for s in s5], 'Vrad': [v for v in vrad],\
        'All': [i for i in a], 'R-Void': [i for i in r], 'S-Void': [i for i in s]}
print(tabulate(info, headers='keys'))



# %%
