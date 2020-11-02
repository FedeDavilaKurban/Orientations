"""
+ Lectura de Datos
+ Creacion de Tabla 'gxs' (restringida en Masa) con Pos., Spin, Masa, SFR
+ Creacion del periodicCKDTree
"""

print 'Reading data...'
voidsFinal = rdr.read('voidsFinal.txt')#,format='commented_header')#,converters=converters)
#subgroups = il.groupcat.loadSubhalos(basePath,135)

gxs = Table(subgroups['SubhaloPos'],names=['x','y','z'],dtype=['float64','float64','float64'])
col_id = Table.Column([i for i in range(len(gxs))],name='ID')
col_sp = Table.Column(subgroups['SubhaloSpin'],name='sp',dtype='float64')
col_m = Table.Column(subgroups['SubhaloMass'],name='mass',dtype='float64')
gxs.add_columns([col_sp,col_m,col_id])


col_spn = Table.Column(np.sum(np.abs(gxs['sp'])**2,axis=-1)**(1./2),name='sp_n',dtype='float64')
gxs.add_column(col_spn,index=5)

gxs = gxs[(np.log10(gxs['mass'])<3.)&(np.log10(gxs['mass'])>-1.)]

lbox = 75000. #kpc
bounds = [lbox,lbox,lbox]
tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)
#---------------------------------------------------------

rmin, rmax = .7, 1.2
perc = 75
#nv = 9
#paramList = []
rantable = Table()
a1_ran = []
a2_ran = []
a3_ran = []
a4_ran = []

sec = 123

for _ in range(500):
	print _
	#Nt = []

	#---------------------------------------------------------
	"""
	+ 'gals' son las galaxias en el shell del void
	"""

	xv,yv,zv,r = np.random.uniform(0,75000),np.random.uniform(0,75000),np.random.uniform(0,75000),np.random.uniform(4000,8000)
#	print xv,yv,zv,r
	
	idx1 = tree.query_ball_point([xv,yv,zv],r*rmax)
	idx2 = tree.query_ball_point([xv,yv,zv],r*rmin)
	shell = [g for g in idx1 if g not in idx2]

	gals=gxs[shell]
#	print 'N of galaxies:',len(gals)


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

	#gals_h = gals[(S > M*m1+b1)&(S < M*m1+b2)&(S < M*m+b)] # Limite superior e inferior en el espacio Spin-Masa
	#gals_h = gals[(S < M*m1+b1)&(S < M*m+b)] # Solo limite superior
	#gals_h = gals[(S > M*m1+b2)&(S < M*m+b)] # Solo limite inferior
	#gals_h = gals[(S < M*m+b)] # Solo limite en Spin
	#gals_h = gals
	#-------------------
	# Limites lineales

	#Alto Spin
	if sec == 1: gals_h = gals[(M < -.7)&(S > M*m+b)] # Solo limite superior
	if sec == 2: gals_h = gals[(M > -.7)&(M < -.3)&(S > M*m+b)] # Limite superior e inferior en el espacio Spin-Masa
	if sec == 3: gals_h = gals[(M > -.3)&(S > M*m+b)] # Solo limite inferior
	if sec == 12: gals_h = gals[(M < -.3)&(S > M*m+b)] # Solo limite superior
	if sec == 23: gals_h = gals[(M < -.7)&(S > M*m+b)] # Solo limite inferior
	if sec == 123: gals_h = gals[(S > M*m+b)] # Solo limite en Spin

	#Bajo Spin
	if sec == 4: gals_h = gals[(M < -.7)&(S < M*m+b)] # Solo limite superior
	if sec == 5: gals_h = gals[(M > -.7)&(M < -.3)&(S < M*m+b)] # Limite superior e inferior en el espacio Spin-Masa
	if sec == 6: gals_h = gals[(M > -.3)&(S < M*m+b)] # Solo limite inferior
	if sec == 45: gals_h = gals[(M < -.3)&(S < M*m+b)] # Solo limite superior
	if sec == 56: gals_h = gals[(M < -.7)&(S < M*m+b)] # Solo limite inferior
	if sec == 456: gals_h = gals[(S < M*m+b)] # Solo limite en Spin

	if sec == 'all': gals_h = gals
	#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	#plt.figure(2)
	#plt.scatter(M,S,c='gray')
	#plt.scatter(np.log10(gals_h['mass']),np.log10(gals_h['sp_n']),c='blue')
	#plt.vlines(-.7,0,3,color='r',linestyle='dashed')
	#plt.vlines(-.3,0,3,color='r',linestyle='dashed')
	#plt.plot(M,M*m+b,ls=':')
	#plt.plot(M,M*g_m+g_b,ls='--')
	#plt.show()

	print 'N of final galaxies:',len(gals_h)
	
#----------------------------------------------------

	"""
	+ Calculo de cosenos de angulos 3D 
	+ Histograma de la Distribucion de Cosenos
	+ Histograma de SFR para gxs con Cosenos cercanos y lejanos a 0
	"""
	cos=[]
	cos1=[]
	for i in range(len(gals_h)):
		u = [gals_h[i]['x']-xv,gals_h[i]['y']-yv,gals_h[i]['z']-zv]
		v = gals_h[i]['sp']
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

	#y = ecdf(cos_h) - cos_h # Diferencia entre la acumulada y la recta 
	y = ecdf(cos_h)-(cos_h+1.)/2. # Para cuando tomamos el valor verdadero de cos (no el absoluto) 
								   # Cambia tmb la 'func' de abajo

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

	a1_ran.append( coeffs[1] )
	a2_ran.append( coeffs[2] )
	a3_ran.append( coeffs[3] )
	a4_ran.append( coeffs[4] )
#	print 'a1=',a1[j],'; a2=',a2[j],'; a3=',a3[j],'; a4=',a4[j]

	#------------------------------------------------------------------------
#	gs = gridspec.GridSpec(3, 1)

#	fig = plt.figure()
#	fig.suptitle('',size=10)

#	ax1 = fig.add_subplot(gs[0,0])
#	ax2 = fig.add_subplot(gs[1,0])
#	ax3 = fig.add_subplot(gs[2,0])

#	ax1.hist(cos1,bins=40,normed=True,histtype='stepfilled',cumulative=False,alpha=30.)
#	ax1.plot(x,d_yfit+.5,'r--')

#	ax2.plot(cos_h,ecdf(cos_h))
#	ax2.plot(cos_h,(cos_h+1.)/2.,'b--')
#	ax2.plot(cos_h,yfit+(cos_h+1.)/2.,'r--')

#	ax3.step(x,y)
#	ax3.plot(x,yfit,'r--')
#	ax3.set_xlabel('cos',size=18)

#	ax1.set_title('cos (high spin)',size=14)
#	ax2.set_title('ECDF',size=14)

#	plt.show()
	#plt.save('')

	#------------------------------------------------------------------------

col_a1 = Table.Column(a1_ran,name='a1_ran',dtype='float64')
col_a2 = Table.Column(a2_ran,name='a2_ran',dtype='float64')
col_a3 = Table.Column(a3_ran,name='a3_ran',dtype='float64')
col_a4 = Table.Column(a4_ran,name='a4_ran',dtype='float64')

rantable.add_columns([col_a1,col_a2,col_a3,col_a4])
#ascii.write(rantable,'ran_coeffs.txt')

