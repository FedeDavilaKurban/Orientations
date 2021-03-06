#-------------------------------------------------
"""
Creo una tabla 'jvm' con Masa, J, SFR, ID
"""

jvm = Table()

groups = il.groupcat.loadHalos(basePath,135)

gxs = Table(groups['GroupPos'],names=['x','y','z'],dtype=['float64','float64','float64'])
#col_sp = Table.Column(groups['GroupSpin'],name='sp',dtype='float64')
col_id = Table.Column([i for i in range(len(gxs))],name='ID')
gxs.add_columns([col_id])

col_m = Table.Column(groups['GroupMass'],name='Mass')
col_sfr = Table.Column(groups['GroupSFR'],name='SFR')
#col_spn = Table.Column(np.sum(np.abs(gxs['sp'])**2,axis=-1)**(1./2),name='sp_n',dtype='float64')

jvm.add_columns([col_m,col_sfr,col_id])
#------------------------------------------------
"""
Graficos J vs Mass
"""


#m2,m3 = .7,.9 
#b2,b3 = 2.1,1.1
x = np.array([-1.,0.,1.,2.,3.,4.])

jvm1 = jvm[(np.log10(jvm['Mass'])<3.)&(np.log10(jvm['Mass'])>-1.)] # Corto la tabla a un rango de masa razonable

m,b,rvalue,pvalue,std=scipy.stats.linregress(np.log10(jvm1['Mass']),np.log10(jvm1['sp_n']))

jvm2 = jvm1[np.log10(jvm1['sp_n'])>(np.log10(jvm1['Mass'])*m+b)] # Gxs de alto spin siguiendo la relacion J vs M
jvm3 = jvm1[np.log10(jvm1['sp_n'])<(np.log10(jvm1['Mass'])*m+b)] # Gxs de bajo spin siguiendo la relacion J vs M

print len(jvm2),len(jvm3)

#plt.scatter(np.log10(jvm2['Mass']),np.log10(jvm2['sp_n']))
#plt.scatter(np.log10(jvm3['Mass']),np.log10(jvm3['sp_n']))
plt.scatter(np.log10(jvm1['Mass']),np.log10(jvm1['sp_n']))

#plt.plot(x,x*m2+b2)
plt.plot(x,x*m+b,ls='--')
plt.plot(x,x*(-1./m)+b)
z = np.array([0.,.8,1.6])
plt.scatter(z,z*m+b,c='r',s=100)
plt.plot(x,x*-1./m+.8*(m+1./m)+b)
plt.plot(x,x*-1./m+1.6*(m+1./m)+b)
plt.xlim([-6,6])
plt.ylim([-6,6])
#plt.plot(x,x*m3+b3)
#plt.plot(x,x*0.66+1.6,ls=':')
plt.show()
#-------------------------------------------------------------
"""
Graficos J vs SFR/M* 
"""

fig = plt.figure()
#fig.suptitle()
gs = gridspec.GridSpec(2, 1)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[1,0],sharex=ax1,sharey=ax1)

#ax1.scatter(np.log10(jvm2['SFR']/jvm2['Mass']),np.log10(jvm2['sp_n']))
#ax2.scatter(np.log10(jvm3['SFR']/jvm3['Mass']),np.log10(jvm3['sp_n']))

ax1.hist(jvm2['SFR']/jvm2['Mass'],bins=1000,normed='False')
ax2.hist(jvm3['SFR']/jvm3['Mass'],bins=1000,normed='False')

ax1.set_xlim([0.,0.03])
ax2.set_xlim([0.,0.03])


ax1.set_title('Alto Spin',size=14)
ax2.set_title('Bajo Spin',size=14)

plt.show()




