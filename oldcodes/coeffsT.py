"""
Hago un bineado de las galaxias en el espacio Spin-Masa, calculo los parametros a1...a4, y los guardo en una tabla

Escrito a partir del coeffs1.py
"""

#----------------------------------------------------
"""
+ Lectura de Datos
+ Creacion de Tabla 'gxs' (restringida en Masa) con Pos., Spin, Masa, SFR
+ Creacion del periodicCKDTree
"""

#a_ran = rdr.read('ran_coeffs_{0}-{1}-{2}.txt'.format(rmin,rmax,perc))#,format='commented_header')
a_ran = rdr.read('ran_coeffs.txt') # Escrito con el programa a_ran.py
#print a_ran
print 'Reading data...'
voidsFinal = rdr.read('voidsFinal.txt')#,format='commented_header')#,converters=converters)
#subgroups = il.groupcat.loadSubhalos(basePath,135)

gxs = Table(subgroups['SubhaloPos'],names=['x','y','z'],dtype=['float64','float64','float64'])
col_id = Table.Column([i for i in range(len(gxs))],name='ID')
col_sp = Table.Column(subgroups['SubhaloSpin'],name='sp',dtype='float64')
col_sfr = Table.Column(subgroups['SubhaloSFR'],name='SFR')
col_m = Table.Column(subgroups['SubhaloMass'],name='mass',dtype='float64')
gxs.add_columns([col_sp,col_m,col_sfr,col_id])


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
paramList = []
a2_ran = a_ran['a2_ran'] 
a4_ran = a_ran['a4_ran'] 

for nv in range(len(voidsFinal)):
	print 'nv=',nv
	Nt = []

	table = Table()
	#---------------------------------------------------------
	"""
	+ 'gals' son las galaxias en el shell del void
	"""

	xv,yv,zv,r = voidsFinal[nv]

	idx1 = tree.query_ball_point([xv,yv,zv],r*rmax)
	idx2 = tree.query_ball_point([xv,yv,zv],r*rmin)
	shell = [g for g in idx1 if g not in idx2]

	gals=gxs[shell]
	#print 'N of galaxies:',len(gals)


	#-------------------------------------------------------
	"""
	+ Seleccion de la muestra segun la relacion J vs Mass
	+ Histogramas de Cosenos para las dos poblaciones (#)
	"""

	m,b,rvalue,pvalue,std=scipy.stats.linregress(np.log10(gxs['mass']),np.log10(gxs['sp_n'])) # El ajuste tiene que ser con las 'gxs' (no con las 'gals')
	m1 = -1./m
	#beta = np.logspace(-1.,np.log10(3.1)/np.log10(10.),6)-1.1
	beta = [-1.,-.7,0.,2.]
	
	M = np.log10(gals['mass'])
	S = np.log10(gals['sp_n'])

	N = []
	beta1 = []
	beta2 = []
	a1 = []
	a2 = []
	a3 = []
	a4 = []
	pvalue2 = []
	pvalue4 = []

	for j in range (len(beta)-1):

		b1 = beta[j]*(m-m1)+b
		b2 = beta[j+1]*(m-m1)+b
			
		gals_h = gals[(S > M*m1+b1)&(S < M*m1+b2)&(S > M*m+b)] # Limite superior e inferior en el espacio Spin-Masa
		
		if len(gals_h) <= 2. : 
			print 'BREAK: N of gxs in bin <= 2' 
			break
		 
		#plt.scatter(M,S,c='gray')
		#plt.scatter(np.log10(gals_h['mass']),np.log10(gals_h['sp_n']),c='blue')
		#plt.plot(M,M*m1+b1,ls=':')
		#plt.plot(M,M*m1+b2,ls=':')
		#plt.plot(M,M*m+b,ls=':')
		#plt.show()

		#print 'N of final galaxies:',len(gals_h)
		
		N.append(len(gals_h))
		Nt.append(len(gals))
		
		beta1.append('%.2f'%beta[j])
		beta2.append('%.2f'%beta[j+1])

		#-----------------------------------------------------------

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
			
		#plt.close()
		#plt.hist(cos,bins=20,histtype='stepfilled')#,cumulative=True)
		#plt.xlabel('$|cos \\theta |$',size=15)
		#plt.show()

		#col_c = Table.Column(cos,name='cos',dtype='float64')
		#gals.add_column(col_c)

		#----------------------------------------------------
		"""
		+ Calculo y Grafico de la ECDF
		+ Fitteado de la resta ECDF-Recta Lineal Random
		+ Coeffs a1 a2

		"""
		cos_h = np.sort(cos1)
		
		# Empirical cumulated distribution function
		ecdf = ECDF(cos_h) 
		
		# Diferencia entre la acumulada y la recta
		#y = ecdf(cos_h) - cos_h  
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

		#a1.append( coeffs[1] )
		a2.append( coeffs[2] )
		#a3.append( coeffs[3] )
		a4.append( coeffs[4] )

		#------------------------------------------------------------------------
		#gs = gridspec.GridSpec(3, 1)

		#fig = plt.figure()
		#fig.suptitle('',size=10)

		#ax1 = fig.add_subplot(gs[0,0])
		#ax2 = fig.add_subplot(gs[1,0])
		#ax3 = fig.add_subplot(gs[2,0])

		#ax1.hist(cos1,bins=40,normed=True,histtype='stepfilled',cumulative=False,alpha=30.)
		#ax1.plot(x,d_yfit+.5,'r--')

		#ax2.plot(cos_h,ecdf(cos_h))
		#ax2.plot(cos_h,(cos_h+1.)/2.,'b--')
		#ax2.plot(cos_h,yfit+(cos_h+1.)/2.,'r--')

		#ax3.step(x,y)
		#ax3.plot(x,yfit,'r--')
		#ax3.set_xlabel('cos',size=18)

		#ax1.set_title('cos (high spin)',size=14)
		#ax2.set_title('ECDF',size=14)

		#plt.show()
		#plt.save('')

		#------------------------------------------------------------------------

		pvalue2.append( float(len([x for x in a2_ran if x<a2[j]])) / float(len(a2_ran)) )
		pvalue4.append( float(len([x for x in a4_ran if x<a4[j]])) / float(len(a4_ran)) )
		
		#print ' a2=',a2[j],'; a4=',a4[j],'; p2=',pvalue2[j],'; p4=',pvalue4[j]
		#print pvalue1, pvalue2

		#fig = plt.figure()
		#ax1 = fig.add_subplot(1, 2, 1)
		#ax2 = fig.add_subplot(1, 2, 2)

		##fig.suptitle('Rmin, Rmax: '+str(rmin)+', '+str(rmax)+'. Percentile: '+str(perc) )

		#ax1.set_title('$a2$',size=18)
		#ax2.set_title('$a4$',size=18)

		##ax1.text(left,top,'pvalue:'+str(pvalue1))#,ha='center', va='center',transform=ax.transAxes)

		#ax1.annotate('P value: '+'%.3f'%pvalue2[j], xy=(0.05, 0.95), xycoords='axes fraction')
		#ax2.annotate('P value: '+'%.3f'%pvalue4[j], xy=(0.05, 0.95), xycoords='axes fraction')

		#n1, bins, patches = ax1.hist(a2_ran,bins=15,normed=True,histtype='stepfilled')
		#n2, bins2, patches2 = ax2.hist(a4_ran,bins=15,normed=True,histtype='stepfilled')

		#ax1.axvline(a2[j],color='r',linestyle='dashed',label='a2')
		#ax2.axvline(a4[j],color='r',linestyle='dashed',label='a4')

		#ax1.legend()
		#ax2.legend()
		#plt.show()
		#------------------------------------------------------------------------

	col_Nt = Table.Column(Nt,name='Nt',dtype='int64',description='Total n of shell gxs')
	col_b1 = Table.Column(beta1,name='beta1',dtype='float64')
	col_b2 = Table.Column(beta2,name='beta2',dtype='float64')
	col_N  = Table.Column(N,name='N',dtype='int64',description='N of bin gxs')
	#col_a1 = Table.Column(a1,name='a1',dtype='float64')
	col_a2 = Table.Column(a2,name='a2',dtype='float64')
	#col_a3 = Table.Column(a3,name='a3',dtype='float64')
	col_a4 = Table.Column(a4,name='a4',dtype='float64')
	col_p2 = Table.Column(pvalue2,name='pvalue2',dtype='float64')
	col_p4 = Table.Column(pvalue4,name='pvalue4',dtype='float64')

	table.add_columns([col_Nt,col_b1,col_b2,col_N,col_a2,col_a4,col_p2,col_p4])

	paramList.append( table )
