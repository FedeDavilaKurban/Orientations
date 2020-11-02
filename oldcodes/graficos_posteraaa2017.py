"""
Programa para ver las orientaciones de galaxias de alto y bajo spin en un solo grafico

Estos graficos se usaron en el poster de la AAA2017

"""
#---------------------------------


rmin, rmax = .7, 1.1
#perc = 75
#nv = 0
#sec = 123
#print 'nv=',nv
#----------------------------------------------------
"""
+ Lectura de Datos
+ Creacion de Tabla 'gxs' (restringida en Masa) con Pos., Spin, Masa, SFR
+ Creacion del periodicCKDTree
"""

print 'Reading data...'
#voidsFinal = ascii.read('voidsFinal_gxstracers.txt')#,format='commented_header')#,converters=converters)
profile = ascii.read('voidsProfileFinal.txt')
#voidsFinal = ascii.read('voids_ill_polaco.dat',names=['r','x','y','z'])#,format='commented_header')#,converters=converters)
#for col in voidsFinal.colnames:
#	voidsFinal[col]*=1000. #Paso a kpc
#flag=1
voidsFinal = ascii.read('voidsFinal_gxstracers.txt')
flag=0

gxs = Table(subgroups['SubhaloPos'],names=['x','y','z'],dtype=['float64','float64','float64'])
#N = len(gxs)
N = len(subgroups['SubhaloPos'][(np.log10(subgroups['SubhaloMass'])>1.)])
rho_m = float(N) / (75000.**3)

col_id = Table.Column([i for i in range(len(gxs))],name='ID')
col_sp = Table.Column(subgroups['SubhaloSpin'],name='sp',dtype='float64')
col_sfr = Table.Column(subgroups['SubhaloSFR'],name='SFR')
col_m = Table.Column(subgroups['SubhaloMass'],name='mass',dtype='float64')
gxs.add_columns([col_sp,col_m,col_sfr,col_id])


col_spn = Table.Column(np.sum(np.abs(gxs['sp'])**2,axis=-1)**(1./2),name='sp_n',dtype='float64')
gxs.add_column(col_spn,index=5)

g_m,g_b,g_rvalue,g_pvalue,g_std=scipy.stats.linregress(np.log10(gxs['mass']),np.log10(gxs['sp_n'])) # Este ajuste es de las galaxias en total (ajuste global "g_...")

gxs = gxs[(np.log10(gxs['mass'])>-1.)]

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

if flag==1: r,x,y,z = voidsFinal[nv]
else: x,y,z,r = voidsFinal[nv]

idx1 = tree.query_ball_point([x,y,z],r*rmax)
idx2 = tree.query_ball_point([x,y,z],r*rmin)
shell = [g for g in idx1 if g not in idx2]

gals = gxs[shell]
print 'N of galaxies:',len(gals)

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
Plot del perfil de densidad del Void
"""
gs = gridspec.GridSpec(3, 1)
fig = plt.figure(figsize=(20,10))#,dpi=100)
fig.suptitle('Void{}_{:.3}Mpc'.format(nv,(voidsFinal['r'][nv])/1000.),size=14)

ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[1,0])
ax3 = fig.add_subplot(gs[2,0])
#ax4 = fig.add_subplot(gs[0,1])

ax2.set_xlabel('Distance [kpc]')
ax2.set_ylabel('Density/MeanDensity - 1')
ax2.plot(profile['r'],profile['Void{}'.format(nv)]/rho_m-1.)

#-------------------------------------------------------
"""
+ Seleccion de la muestra segun la relacion J vs Mass
+ Histogramas de Cosenos para las dos poblaciones (#)
"""
m,b,rvalue,pvalue,std=scipy.stats.linregress(np.log10(gxs['mass']),np.log10(gxs['sp_n'])) # El ajuste tiene que ser con las 'gxs' (no con las 'gals')
m1 = -1./m
b1 = -.7*(m-m1)+b
b2 = -.3*(m-m1)+b

b_s = b
b_i = b


M = np.log10(gals['mass'])
S = np.log10(gals['sp_n'])

#::::::::::::::::::::::::::::::::
# Limites en el espacio Spin-Masa

#Alto Spin
#if sec == 123: gals_h = gals[(S > M*m+b)] # Solo limite en Spin

#Bajo Spin
#if sec == 456: gals_h = gals[(S < M*m+b)] # Solo limite en Spin

#gals_h = gals
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ax3 = fig.add_subplot(gs[1,1])
ax3.scatter(np.log10(gals['mass']),np.log10(gals['sp_n']),c='blue', label='Galaxies Being Analyzed')
#ax3.plot(M,M*m1+b1,ls=':')
#ax3.vlines(-.7,0,3,color='r',linestyle='dashed')
#ax3.plot(M,M*m1+b2,ls=':')
#ax3.vlines(-.3,0,3,color='r',linestyle='dashed')
ax3.plot(M,M*m+b,ls=':', label = 'Void Galaxies Low/High Spin Linear Regression')
ax3.plot(M,M*m+b_s,ls=':')
ax3.plot(M,M*m+b_i,ls=':')

ax3.plot(M,M*g_m+g_b,ls='--', label = 'Box Galaxies Low/High Spin Linear Regression')
ax3.set_xlabel('log10( Mass/SolarMass )')
ax3.set_ylabel('Spin')
ax3.legend(loc=4)


#print 'N of final galaxies:',len(gals_h)
#-----------------------------------------------------------
for sec in [123,456]:
	x,y,z,r = voidsFinal[nv]
	if sec == 123: gals_h = gals[(S > M*m+b_s)]
	else: gals_h = gals[(S < M*m+b_i)]
	print sec
	"""
	+ Calculo de cosenos de angulos 3D 
	+ Histograma de la Distribucion de Cosenos
	+ Histograma de SFR para gxs con Cosenos cercanos y lejanos a 0
	"""
	cos=[]
	cos1=[]
	for i in range(len(gals_h)):
		u = [gals_h[i]['x']-x,gals_h[i]['y']-y,gals_h[i]['z']-z]
		v = gals_h[i]['sp']
		#print 'U:',u,'V:',v
		c = abs( np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)) )
		cos.append( c )
		cos1.append( c )
		cos1.append( -c )
		#cos.append( np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)) )
	
	#----------------------------------------------------
	"""
	+ Calculo y Grafico de la ECDF
	+ Fitteado de la resta ECDF-Recta Lineal Random
	+ Coeffs a1 a2

	"""
	cos_h = np.sort(cos1)

	ecdf = ECDF(cos_h) # Empirical cumulated distribution function
	y = ecdf(cos_h)-(cos_h+1.)/2. # Para cuando tomamos el valor verdadero de cos (no el absoluto) 
								   # Cambia tmb la 'func' de abajo

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
	a3 = coeffs[3]
	a4 = coeffs[4]
	print 'a1=',a1,'; a2=',a2,'; a3=',a3,'; a4=',a4
	
	if sec==123: ls,label=['b--','High Spin']
	else: ls,label=['g--','Low Spin']
	ax1.step(x,y)
	ax1.plot(x,yfit,ls,label=label)

#------------------------------------------------------------------------


#ax1.hist(cos1,bins=40,normed=True,histtype='stepfilled',cumulative=False,alpha=30.)
#ax1.plot(x,d_yfit+.5,'r--')
#ax1.set_xlabel('$cos \\theta (High Spin)$',size=18)
#ax1.legend()

#ax2.plot(cos_h,ecdf(cos_h))
#ax2.plot(cos_h,(cos_h+1.)/2.,'b--')
#ax2.plot(cos_h,yfit+(cos_h+1.)/2.,'r--')
#ax2.set_xlabel('cos',size=18)

#ax3.set_xlabel('cos',size=18)


#ax1.set_title('Cos (Mirrored)',size=12)
#ax2.set_title('ECDF',size=12)
ax1.legend()
ax1.set_title('Residues',size=12)
ax1.set_ylim(-.05,.05)

#plt.savefig('PlotsAlignments/alignment_Void{}_Sec{}.pdf'.format(nv,sec))
#plt.savefig('PlotsAlignments/alignment_Void{}_Sec{}.png'.format(nv,sec))
#plt.savefig('posterAAA2017_Void{}_LowAndHighSpin.svg'.format(nnv))
plt.show()
plt.close()

