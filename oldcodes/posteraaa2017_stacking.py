"""

"""
#---------------------------------


rmin, rmax = 0.7, 1.1

voidsFinal = ascii.read('voidsFinal_gxstracers.txt')#,format='commented_header')#,converters=converters)
#voidsFinal = ascii.read('voids_ill_polaco.dat')#,format='commented_header')#,converters=converters)
#for col in voidsFinal.colnames:
#	voidsFinal[col]*=1000. #Paso a kpc

N = range(len(voidsFinal))
#sec = 456

#----------------------------------------------------
"""
+ Lectura de Datos
+ Creacion de Tabla 'gxs' (restringida en Masa) con Pos., Spin, Masa, SFR
+ Creacion del periodicCKDTree
"""

print 'Reading data...'

#profile = ascii.read('voidsProfileFinal.txt')

gxs = Table(subgroups['SubhaloPos'],names=['x','y','z'],dtype=['float64','float64','float64'])
col_id = Table.Column([i for i in range(len(gxs))],name='ID')
col_sp = Table.Column(subgroups['SubhaloSpin'],name='sp',dtype='float64')
col_sfr = Table.Column(subgroups['SubhaloSFR'],name='SFR')
col_m = Table.Column(subgroups['SubhaloMass'],name='mass',dtype='float64')
gxs.add_columns([col_sp,col_m,col_sfr,col_id])


col_spn = Table.Column(np.sum(np.abs(gxs['sp'])**2,axis=-1)**(1./2),name='sp_n',dtype='float64')
gxs.add_column(col_spn,index=5)

g_m,g_b,g_rvalue,g_pvalue,g_std=scipy.stats.linregress(np.log10(gxs['mass']),np.log10(gxs['sp_n'])) # Este ajuste es de las galaxias en total (ajuste global "g_...")

gxs = gxs[(np.log10(gxs['mass'])>-1.)]


lbox = 75000. #kpc
bounds = [lbox,lbox,lbox]
tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)


"""
+ Calculo de cosenos de angulos 3D 
+ Histograma de la Distribucion de Cosenos
+ Histograma de SFR para gxs con Cosenos cercanos y lejanos a 0
"""
cosh=[]
cosl=[]
for nv in N:
	print nv
	#---------------------------------------------------------
	"""
	+ 'gals' son las galaxias en el shell del void
	"""

	x,y,z,r = voidsFinal[nv]

	idx1 = tree.query_ball_point([x,y,z],r*rmax)
	idx2 = tree.query_ball_point([x,y,z],r*rmin)
	shell = [g for g in idx1 if g not in idx2]

	gals = gxs[shell]
	print 'N of galaxies:',len(gals)

	#-------------------------------------------------------
	"""
	Plot del perfil de densidad del Void
	"""
	#plt.figure(1)
	#plt.plot(profile['r'],profile['Void{}'.format(nv)])

	#-------------------------------------------------------
	"""
	+ Seleccion de la muestra segun la relacion J vs Mass
	+ Histogramas de Cosenos para las dos poblaciones (#)
	"""
	m,b,rvalue,pvalue,std=scipy.stats.linregress(np.log10(gxs['mass']),np.log10(gxs['sp_n'])) # El ajuste tiene que ser con las 'gxs' (no con las 'gals')
	m1 = -1./m
	b1 = -.7*(m-m1)+b
	b2 = -.3*(m-m1)+b

	M = np.log10(gals['mass'])
	S = np.log10(gals['sp_n'])

	#::::::::::::::::::::::::::::::::
	# Limites en el espacio Spin-Masa
	#Alto Spin
	gals_h = gals[(S > M*m+b)] # Solo limite en Spin
	for i in range(len(gals_h)):
		u = [gals_h[i]['x']-x,gals_h[i]['y']-y,gals_h[i]['z']-z]
		v = gals_h[i]['sp']
		c = abs( np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)) )
		#cos.append( c )
		cosh.append( c )
		cosh.append( -c )

	#Bajo Spin
	gals_h = gals[(S < M*m+b)] # Solo limite en Spin
	for i in range(len(gals_h)):
		u = [gals_h[i]['x']-x,gals_h[i]['y']-y,gals_h[i]['z']-z]
		v = gals_h[i]['sp']
		c = abs( np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)) )
		#cos.append( c )
		cosl.append( c )
		cosl.append( -c )

	#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	
#----------------------------------------------------
"""
+ Calculo y Grafico de la ECDF
+ Fitteado de la resta ECDF-Recta Lineal Random
+ Coeffs a1 a2

"""

#------
# Genero datos artificiales para plot del poster AAA2017
#gauss = np.random.normal(0.,.6,100000)
#gauss = gauss[abs(gauss)<1.]
#gauss1 = gauss[(gauss>0.)]*-1.+1.
#gauss2 = gauss[(gauss<0.)]*-1.-1.
#gaussinv = np.concatenate((gauss1,gauss2)) 

#unif = np.random.uniform(-1.,1.,100000)


#cos1 = unif
#----

gs = gridspec.GridSpec(1, 1)
fig = plt.figure(figsize=(20,10))#,dpi=100)
ax3 = fig.add_subplot(gs[0,0])
ax3.set_xlim(-1.,1.)
ax3.set_ylim(-.015,.015)

for cos1 in [cosh,cosl]:
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

	if cos1==cosh: ls,label=['b--','High Spin']
	else: ls,label=['g--','Low Spin']
	ax3.step(x,y)
	ax3.plot(x,yfit,ls,label=label)




#if N == range(len(voidsFinal)):
#	fig.suptitle('AllVoids_Sec{}'.format(sec),size=14)
#else:
#	fig.suptitle('Voids{}_Sec{}'.format(N,sec),size=14)


#ax1 = fig.add_subplot(gs[0,0])
#ax2 = fig.add_subplot(gs[1,0])


#errN = ax1.hist(cos1,bins=40,normed=True,histtype='stepfilled',cumulative=False,alpha=30.)
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

ax3.text(0.6,-0.013,'rmin={0}, rmax={1}'.format(rmin,rmax),fontsize=18)
ax3.legend()
ax3.set_title('Residues',size=12)

plt.savefig('posterAAA2017_AllVoids_HighAndLowSpin.svg')
plt.show()



