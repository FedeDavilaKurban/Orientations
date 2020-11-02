
rmin, rmax = .7, 1.2
perc = 75
nv = 13
print 'nv=',nv
#----------------------------------------------------
"""
+ Lectura de Datos
+ Creacion de Tabla 'gxs' (restringida en Masa) con Pos., Spin, Masa, SFR
+ Creacion del periodicCKDTree
"""

#a_ran = rdr.read('ran_coeffs_{0}-{1}-{2}.txt'.format(rmin,rmax,perc))#,format='commented_header')
#a_ran = rdr.read('ran_coeffs.txt')
#print a_ran
print 'Reading data...'
voidsFinal = rdr.read('voidsFinal.txt')#,format='commented_header')#,converters=converters)
subgroups = il.groupcat.loadSubhalos(basePath,135)

gxs = Table(subgroups['SubhaloPos'],names=['x','y','z'],dtype=['float64','float64','float64'])
col_id = Table.Column([i for i in range(len(gxs))],name='ID')
col_sp = Table.Column(subgroups['SubhaloSpin'],name='sp',dtype='float64')
col_sfr = Table.Column(subgroups['SubhaloSFR'],name='SFR')
col_m = Table.Column(subgroups['SubhaloMass'],name='mass',dtype='float64')
gxs.add_columns([col_sp,col_m,col_sfr,col_id])


col_spn = Table.Column(np.sum(np.abs(gxs['sp'])**2,axis=-1)**(1./2),name='sp_n',dtype='float64')
gxs.add_column(col_spn,index=5)

gxs = gxs[(np.log10(gxs['mass'])<3.)&(np.log10(gxs['mass'])>-1.)]

#converters = {'x': [ascii.convert_numpy(np.float64)],
#	      'y': [ascii.convert_numpy(np.float64)],
#	      'z': [ascii.convert_numpy(np.float64)],
#	      'r': [ascii.convert_numpy(np.float64)]}

lbox = 75000. #kpc
bounds = [lbox,lbox,lbox]
tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)


#---------------------------------------------------------
"""
+ 'gals' son las galaxias en el shell del void
"""

x,y,z,r = voidsFinal[nv]

idx1 = tree.query_ball_point([x,y,z],r*rmax)
idx2 = tree.query_ball_point([x,y,z],r*rmin)
shell = [g for g in idx1 if g not in idx2]

gals=gxs[shell]

#plt.close()
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.quiver(gals['x'],gals['y'],gals['z'],gals['sp_x'],gals['sp_y'],gals['sp_z'],length=1000.,color='Tomato')
#ax.set_title('3D Vector Field')             # title
#ax.view_init(elev=18, azim=30)              # camera elevation and angle
#ax.dist=8                                   # camera distance
#plt.show()


#-------------------------------------------------------
"""
+ Calculo de cosenos de angulos 3D 
+ Histograma de la Distribucion de Cosenos
+ Histograma de SFR para gxs con Cosenos cercanos y lejanos a 0
"""
cos=[]
for i in range(len(gals)):
	u = [gals[i]['x']-x,gals[i]['y']-y,gals[i]['z']-z]
	v = gals[i]['sp']
	c = ( np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)) )
	cos.append( c )
	#cos.append( -c )
	#cos.append( np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)) )
	
plt.close()
plt.hist(cos,bins=20,histtype='stepfilled')#,cumulative=True)
plt.xlabel('$|cos \\theta |$',size=15)
plt.show()

col_c = Table.Column(cos,name='cos',dtype='float64')
gals.add_column(col_c)

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

n1, bins, patches = ax1.hist(gals['SFR'][(gals['cos']>-.5)&(gals['cos']<.5)],bins=15,normed=False,histtype='stepfilled',cumulative=False,alpha=30.)

n2, bins1, patches1 = ax2.hist(gals['SFR'][(gals['cos']<-.5)|(gals['cos']>.5)],bins=15,normed=False,histtype='stepfilled',cumulative=False,alpha=30.)

plt.show()


#----------------------------------------------------
"""
+ Separacion entre Alto y Bajo J segun la relacion J vs Mass
+ Histogramas de Cosenos para las dos poblaciones
"""

#m,b = 0.65202943253484447, 1.840361622864827 # Sacados con stats.linregress() de todas las gxs con 0.<log10(Mass)<3. en un grafico J vs Mass
m,b = 0.70249417767041611, 1.79943480137090320 # Sacados con stats.linregress() de todas las gxs con -1.<log10(Mass)<3. en un grafico J vs Mass

gals_h = gals[np.log10(gals['sp_n'])>(np.log10(gals['mass'])*m+b)] # Gxs de alto spin siguiendo la relacion J vs M
gals_l = gals[np.log10(gals['sp_n'])<(np.log10(gals['mass'])*m+b)] # Gxs de bajo spin siguiendo la relacion J vs M

#execfile('spin.py')
	
#plt.close()

#fig = plt.figure()
#ax1 = fig.add_subplot(1, 2, 1)
#ax2 = fig.add_subplot(1, 2, 2)

#n1, bins, patches = ax1.hist(gals_h['cos'],bins=15,normed=False,histtype='stepfilled',cumulative=False,alpha=30.)
#ax1.set_xlabel('$cos \\theta (High Spin)$',size=18)

#n2, bins1, patches1 = ax2.hist(gals_l['cos'],bins=15,normed=False,histtype='stepfilled',cumulative=False,alpha=30.)
#ax2.set_xlabel('$cos \\theta (Low Spin)$',size=18)

#plt.show()
#-----------------------------------------------------------
"""
+ Calculo y Grafico de la ECDF
+ Fitteado de la resta ECDF-Recta Lineal Random
+ Coeffs a1 a2

"""

cos_h = np.sort(gals_h['cos'])

ecdf = ECDF(cos_h) # Empirical cumulated distribution function

#y = ecdf(cos_h) - cos_h # Diferencia entre la acumulada y la recta 
y = ecdf(cos_h)-(cos_h+1.)/2. # Para cuando tomamos el valor verdadero de cos (no el absoluto) 
							   # Cambia tmb la 'func' de abajo

#plt.step(cos_h,y)
#plt.show()

#def func(x,a,b,c,d,e):
#	return a + b*np.sin( np.pi*x ) + c*np.sin( 2.*np.pi*x ) + d*np.sin( 3.*np.pi*x ) + e*np.sin( 4.*np.pi*x )
#def dfunc(x,b,c,d,e):
#	return np.pi*b*np.cos( np.pi*x ) + 2.*np.pi*c*np.cos( 2.*np.pi*x ) + 3.*np.pi*d*np.cos( 3.*np.pi*x ) + 4.*np.pi*e*np.cos( 4.*np.pi*x )

def func(x,a,b,c,d,e):
	return a + b*np.sin( np.pi*(x+1.)/2. ) + c*np.sin( 2.*np.pi*(x+1.)/2. ) + d*np.sin( 3.*np.pi*(x+1.)/2. ) + e*np.sin( 4.*np.pi*(x+1.)/2. )
def dfunc(x,b,c,d,e):
	return np.pi/2.*b*np.cos( np.pi*(x+1.)/2. ) + np.pi*c*np.cos( 2.*np.pi*(x+1.)/2. ) + np.pi*3./2.*d*np.cos( 3.*np.pi*(x+1.)/2. ) + np.pi*2.*e*np.cos( 4.*np.pi*(x+1.)/2. )


x = np.array(cos_h)

coeffs, cov = curve_fit(func, x, y)
yfit = func(x,coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4])
d_yfit = dfunc(x,coeffs[1],coeffs[2],coeffs[3],coeffs[4])

a1 = coeffs[1]
a2 = coeffs[2]
print a1,a2

plt.close()

#plt.plot(x,y)
#plt.plot(x,yfit,'r--')
#plt.show()

#------------------------------------------------------------------------
gs = gridspec.GridSpec(3, 1)

fig = plt.figure()
fig.suptitle('',size=10)

ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[1,0])
ax3 = fig.add_subplot(gs[2,0])

ax1.hist(gals_h['cos'],bins=15,normed=True,histtype='stepfilled',cumulative=False,alpha=30.)
ax1.plot(x,d_yfit+.5,'r--')
#ax1.set_xlabel('$cos \\theta (High Spin)$',size=18)
#ax1.legend()

ax2.plot(cos_h,ecdf(cos_h))
ax2.plot(cos_h,(cos_h+1.)/2.,'b--')
ax2.plot(cos_h,yfit+(cos_h+1.)/2.,'r--')
#ax2.set_xlabel('cos',size=18)

ax3.step(x,y)
ax3.plot(x,yfit,'r--')
ax3.set_xlabel('cos',size=18)


ax1.set_title('cos (high spin)',size=14)
ax2.set_title('ECDF',size=14)

plt.show()
#plt.save('')
#------------------------------------------------------------------------
#a1_ran = a_ran['a1 {}'.format(nv)] 
#a2_ran = a_ran['a2 {}'.format(nv)] 

pvalue1 =  float(len([x for x in a1_ran if x>a1])) / float(len(a1_ran)) 
pvalue2 =  float(len([x for x in a2_ran if x<a2])) / float(len(a2_ran)) 

print pvalue1, pvalue2

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

fig.suptitle('Rmin, Rmax: '+str(rmin)+', '+str(rmax)+'. Percentile: '+str(perc) )

ax1.set_title('$a1$ Coefficient',size=18)
ax2.set_title('$a2$ Coefficient',size=18)

#ax1.text(left,top,'pvalue:'+str(pvalue1))#,ha='center', va='center',transform=ax.transAxes)

ax1.annotate('P value: '+str(round(pvalue1,2)), xy=(0.05, 0.95), xycoords='axes fraction')
ax2.annotate('P value: '+str(round(pvalue2,2)), xy=(0.05, 0.95), xycoords='axes fraction')

n1, bins, patches = ax1.hist(a1_ran,bins=15,normed=True,histtype='stepfilled')
n2, bins2, patches2 = ax2.hist(a2_ran,bins=15,normed=True,histtype='stepfilled')

ax1.axvline(a1,color='r',linestyle='dashed',label='a1')
ax2.axvline(a2,color='r',linestyle='dashed',label='a2')

ax1.legend()
ax2.legend()
plt.show()





