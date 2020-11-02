"""
Genera tabla con deltaMax, rVoid, a1, a2, p1, p2
Hace el T-test para la diferencia de medias (<a1> - <a2> o <p1> - <p2>)
"""
#voidsProfile = ascii.read('voidsProfile.txt',format='commented_header')
#converters = {'x': [ascii.convert_numpy(np.float64)],
#	      'y': [ascii.convert_numpy(np.float64)],
#	      'z': [ascii.convert_numpy(np.float64)],
#	      'r': [ascii.convert_numpy(np.float64)]}
#print 'Reading data...'
#voidsFinal = ascii.read('voidsFinal.txt',format='commented_header',converters=converters)

#-------------------DELTA MAX----------------------------

#rho_m = 1.8284091259259258e-08 # Para Halos
rho_m = 2.3566933333333334e-10 # Para gxs cortadas por masa
deltaMax = []
for i in range(len(voidsFinal)):
	r = voidsFinal[i]['r']
	dRho = voidsProfile['Void{}'.format(i)]/rho_m-1.
	deltaMax.append( dRho[(voidsProfile['r']>3*r)&(voidsProfile['r']>r)].max() )
	
col_dmax = Table.Column(deltaMax,'dmax')
#----------------------------------------------

table = ascii.read('pvalue_{0}-{1}-{2}.txt'.format(rmin,rmax,perc))

col_id = Table.Column([i for i in range(len(voidsFinal))],'ID')

table.add_columns([col_dmax,voidsFinal['r'],col_id],indexes=[0,0,4])
#--------------------------------------------------

param = 'a1'
mask = (table['dmax']<.1)&(table['dmax']>-.1)
#plt.scatter(table['dmax'],table[param])
#plt.scatter(table['dmax'][mask],table[param][mask],color='gray')
#plt.xlabel('$\delta Max$',size=18)
#plt.ylabel(param,size=18)
#plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)

ax1.scatter(table['dmax'],table['a1'])
ax1.scatter(table['dmax'][mask],table['a1'][mask],color='gray')

ax2.scatter(table['dmax'],table['a2'])
ax2.scatter(table['dmax'][mask],table['a2'][mask],color='gray')

ax3.scatter(table['dmax'],table['pvalue1'])
ax3.scatter(table['dmax'][mask],table['pvalue1'][mask],color='gray')

ax4.scatter(table['dmax'],table['pvalue2'])
ax4.scatter(table['dmax'][mask],table['pvalue2'][mask],color='gray')

ax1.set_title('a1',size=10)
ax2.set_title('a2',size=10)
ax3.set_title('pvalue 1',size=10)
ax4.set_title('pvalue 2',size=10)

plt.show()
#-------------------------------------------------

# T-Test
svoids = table[table['dmax']>.05]
rvoids = table[table['dmax']<-.05]

scipy.stats.ttest_ind(svoids['a1'],rvoids['a1'],equal_var='False')
scipy.stats.ttest_ind(svoids['a2'],rvoids['a2'],equal_var='False')

