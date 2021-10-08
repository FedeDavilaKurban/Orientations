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
Na = []
Nr = []
Ns = []

for vfilter in ['-','hi','lo']:

    for sfilter in ['-','hi','lo']:

        for sec in [0,123,1,3,456,4,6,14,36]:

            zeta = []
            zeta_std = []
            ns = []

            for vtype in ['a','r','s']:

                if (vfilter=='-')&(sfilter=='-'):
                    filename = '../data/eta/eta_minradV{}_maxradV{}_sec{}_vtype{}_forTable.txt'\
                        .format(minradV,maxradV,sec,vtype)
                elif (vfilter=='-')&(sfilter!='-'):
                    filename = '../data/eta/eta_sfilter{}_minradV{}_maxradV{}_sec{}_vtype{}_forTable.txt'\
                        .format(sfilter,minradV,maxradV,sec,vtype)
                elif (vfilter!='-')&(sfilter=='-'):
                    filename = '../data/eta/eta_vfilter{}_minradV{}_maxradV{}_sec{}_vtype{}_forTable.txt'\
                        .format(vfilter,minradV,maxradV,sec,vtype)
                elif (vfilter!='-')&(sfilter!='-'):
                    filename = '../data/eta/eta_sfilter{}_vfilter{}_minradV{}_maxradV{}_sec{}_vtype{}_forTable.txt'\
                        .format(sfilter,vfilter,minradV,maxradV,sec,vtype)


                etaTable = ascii.read(filename)

                eta0 = 1./(np.sqrt(2)-1)
                eta_ran_std = np.sqrt(28.1421/etaTable['N'][0])

                zeta.append( (etaTable['eta'][0]-eta0)/eta_ran_std )
                zeta_std.append(etaTable['eta_std'][0]/eta_ran_std)
                #if zeta[-1]<3: break

                ns.append(etaTable['N'])

            if sec==0: 
                spin.append('& -')
                mass.append('& -')
            if sec==1: 
                spin.append('& H')
                mass.append('& L')
            if sec==3: 
                spin.append('& H')
                mass.append('& H')
            if sec==4: 
                spin.append('& L')
                mass.append('& L')
            if sec==6: 
                spin.append('& L')
                mass.append('& H')
            if sec==14: 
                spin.append('& -')
                mass.append('& L')
            if sec==36: 
                spin.append('& -')
                mass.append('& H')
            if sec==123: 
                spin.append('& H')
                mass.append('& -')
            if sec==456: 
                spin.append('& L')
                mass.append('& -')
            if vfilter == 'hi': vrad.append('& H')
            if vfilter == 'lo': vrad.append('& L')
            if vfilter == '-': vrad.append('& -')

            if sfilter == 'hi': s5.append('& H')
            if sfilter == 'lo': s5.append('& L')
            if sfilter == '-': s5.append('& -')

            #a.append( '& {:.1f}+-{:.1f}'.format(zeta[0],zeta_std[0]) )
            
            if zeta[0]>=3.:
                a.append( f'& \\resaltar{{{zeta[0]:4.1f} \\textpm {zeta_std[0]:3.1f}}}' )
            else: 
                a.append( f'& {zeta[0]:4.1f} \\textpm {zeta_std[0]:3.1f}' )
            
            if zeta[1]>=3.:
                r.append( f'& \\resaltar{{{zeta[1]:4.1f} \\textpm {zeta_std[1]:3.1f}}}' )
            else: 
                r.append( f'& {zeta[1]:4.1f} \\textpm {zeta_std[1]:3.1f}' )

            if zeta[2]>=3.:
                s.append( f'& \\resaltar{{{zeta[2]:4.1f} \\textpm {zeta_std[2]:3.1f}}}' )
            else: 
                s.append( f'& {zeta[2]:4.1f} \\textpm {zeta_std[2]:3.1f}' )
                
            Na.append('& {:.0f}'.format(ns[0][0]))
            Nr.append('& {:.0f}'.format(ns[1][0]))
            Ns.append('& {:.0f} \\\\'.format(ns[2][0]))

subsample = []
for i in range(len(a)):
    subsample.append('Ss'+str(i)) 

info = {'Subsample': [s for s in subsample],\
        'Spin': [s for s in spin], 'Mass': [m for m in mass], r'$\Sigma_5$': [s for s in s5], 'Vrad': [v for v in vrad],\
        'All': [i for i in a], 'R-Void': [i for i in r], 'S-Void': [i for i in s],
        'N_a': [i for i in Na], 'N_r': [i for i in Nr],'N_s': [i for i in Ns]}
print(tabulate(info, headers='keys'))



# %%
