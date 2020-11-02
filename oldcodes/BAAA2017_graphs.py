"""
Programa para hacer todos los graficos que van al BAAA2017
Excepto el grafico del ejemplo del metodo: posteraaa2017_ejemplos.py

"""
#---------------------------------

fs = 20
height = 6
figsize= (height*1.618,height)

rmin, rmax = .7, 1.1
#-------------------------------------------------------
"""
+ Lectura de Datos
+ Creacion de Tabla 'gxs' (restringida en Masa) con Pos., Spin, Masa, SFR
+ Creacion del periodicCKDTree
"""

print 'Reading data...'
#voidsFinal = ascii.read('voidsFinal_gxstracers.txt')#,format='commented_header')#,converters=converters)
profile = ascii.read('voidsProfileFinal.txt')
voidsFinal = ascii.read('voidsFinal_gxstracers.txt')

massfilter = (np.log10(subgroups['SubhaloMass'])>-1.)

gxs = Table(subgroups['SubhaloPos'],names=['x','y','z'],dtype=['float64','float64','float64'])
#N = len(gxs)
#N = len(subgroups['SubhaloPos'][(np.log10(subgroups['SubhaloMass'])>-1.)])
#rho_m = float(N) / (75000.**3)

col_id = Table.Column([i for i in range(len(gxs))],name='ID')
col_sp = Table.Column(subgroups['SubhaloSpin'],name='sp',dtype='float64')
col_sfr = Table.Column(subgroups['SubhaloSFR'],name='SFR')
col_m = Table.Column(subgroups['SubhaloMass'],name='mass',dtype='float64')
gxs.add_columns([col_sp,col_m,col_sfr,col_id])


col_spn = Table.Column(np.sum(np.abs(gxs['sp'])**2,axis=-1)**(1./2),name='sp_n',dtype='float64')
gxs.add_column(col_spn,index=5)

gxs = gxs[(np.log10(gxs['mass'])>-1.)]

g_m,g_b,g_rvalue,g_pvalue,g_std=scipy.stats.linregress(np.log10(gxs['mass']),np.log10(gxs['sp_n'])) # Este ajuste es de las galaxias en total (ajuste global "g_...")



#converters = {'x': [ascii.convert_numpy(np.float64)],
#	      'y': [ascii.convert_numpy(np.float64)],
#	      'z': [ascii.convert_numpy(np.float64)],
#	      'r': [ascii.convert_numpy(np.float64)]}

lbox = 75000. #kpc
bounds = [lbox,lbox,lbox]
tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)
nv=0

#---------------------------------------------------------
"""Plot Profiles"""

rho_m = float(len(gxs)) / (lbox**3)

plt.figure(3,dpi=100,figsize=figsize)
plt.tick_params(labelsize=fs)

dashList = [(5,2,20,2),(5,2),(2,5),(4,10),(3,3,2,2)] 

i=0
for void in profile.keys():
	if void=='r': break
	dash = dashList[i]
	plt.plot(profile['r']/1000.,(profile[void]/rho_m)-1,
		ls='--',
		dashes=dash,
		label='{} - Radius: {:.3}Mpc'.format(void,voidsFinal['r'][i]/1000.))
	i+=1
plt.xlabel(r'$r\,\, [Mpc]$',size=fs)
plt.ylabel(r'$\rho(<r)/\bar{\rho}-1$',size=fs)
plt.legend(loc=4,fontsize=16)
plt.tight_layout()
plt.savefig('BAAA_profiles.pdf')
plt.show()

#---------------------------------------------------------


#for _ in range(1):
for _ in range(len(voidsFinal)): 
	
	plt.figure(1,figsize=(9,4))#,dpi=500)
	#if nv==4:
	#	plt.close()
	#	plt.figure(1,figsize=(9,5))
	plt.tick_params(labelsize=fs)
	plt.xticks([])
	#plt.box('on')

	plt.title('Void {}'.format(nv),size=fs,position=(plt.xlim()[1]/2.,plt.ylim()[1]*.9))
	#if nv==4:
	plt.xlabel(r'$cos(\theta)$',fontsize=fs)
	plt.xticks([-1.,-.5,0.,.5,1.])

	plt.ylabel('Residues',fontsize=fs)
	#---------------------------------------------------------
	"""
	+ 'gals' son las galaxias en el shell del void
	"""
	
	x = voidsFinal[nv]['x']
	y = voidsFinal[nv]['y']
	z = voidsFinal[nv]['z']
	r = voidsFinal[nv]['r']


	idx1 = tree.query_ball_point([x,y,z],r*rmax)
	idx2 = tree.query_ball_point([x,y,z],r*rmin)
	shell = [g for g in idx1 if g not in idx2]

	gals = gxs[shell]
	print 'N of galaxies:',len(gals)

	#-------------------------------------------------------
	"""
	+ Seleccion de la muestra segun la relacion J vs Mass
	+ Histogramas de Cosenos para las dos poblaciones (#)
	"""
	m,b,rvalue,pvalue,std=scipy.stats.linregress(np.log10(gals['mass']),np.log10(gals['sp_n'])) # El ajuste tiene que ser con las 'gxs' (no con las 'gals')
	m1 = -1./m
	b1 = -.7*(m-m1)+b
	b2 = -.3*(m-m1)+b

	b_s = b
	b_i = b


	M = np.log10(gals['mass'])
	S = np.log10(gals['sp_n'])

	#-----------------------------------------------------------
	for sec in [123,456]:
		x = voidsFinal[nv]['x']
		y = voidsFinal[nv]['y']
		z = voidsFinal[nv]['z']
		r = voidsFinal[nv]['r']

		if sec == 123: 
			gals_h = gals[(S > M*m+b_s)]
			n_high = len(gals_h)
		else: 
			gals_h = gals[(S < M*m+b_i)]
			n_low = len(gals_h)
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
		y1 = ecdf(cos_h)-(cos_h+1.)/2. # Para cuando tomamos el valor verdadero de cos (no el absoluto) 
									   # Cambia tmb la 'func' de abajo

		def func(x,a,b,c,d,e):
			return a + b*np.sin( np.pi*(x+1.)/2. ) + c*np.sin( 2.*np.pi*(x+1.)/2. ) + d*np.sin( 3.*np.pi*(x+1.)/2. ) + e*np.sin( 4.*np.pi*(x+1.)/2. )
		def dfunc(x,b,c,d,e):
			return np.pi/2.*b*np.cos( np.pi*(x+1.)/2. ) + np.pi*c*np.cos( 2.*np.pi*(x+1.)/2. ) + np.pi*3./2.*d*np.cos( 3.*np.pi*(x+1.)/2. ) + np.pi*2.*e*np.cos( 4.*np.pi*(x+1.)/2. )

		x1 = np.array(cos_h)

		coeffs, cov = curve_fit(func, x1, y1)
		yfit = func(x1,coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4])
		d_yfit = dfunc(x1,coeffs[1],coeffs[2],coeffs[3],coeffs[4])

		a1 = coeffs[1]
		a2 = coeffs[2]
		a3 = coeffs[3]
		a4 = coeffs[4]
		print 'a1=',a1,'; a2=',a2,'; a3=',a3,'; a4=',a4
	
		if sec==123: 
			ls,label=['b--','High Spin']
			nhigh = len(gals_h)
		else: 
			ls,label=['g-.','Low Spin']
			nlow = len(gals_h)
		plt.step(x1,y1)

		plt.plot(x1,yfit,ls,label=label)
	
	#----------------------------------------------------
	"""
	+ PARA RANDOMS
	+ Calculo y Grafico de la ECDF
	+ Fitteado de la resta ECDF-Recta Lineal Random
	+ Coeffs a1 a2

	"""
	print 'Calculating randoms'
	#nran = np.min((n_high,n_low))
	nran = n_high
	
#	y1_slice1 = []
#	y1_slice2 = []
#	y1_slice3 = []
	y1_slices = Table()
	for i in range(100):
	#	print '1'
		cos_ = np.random.uniform(0.,1.,nran)
		cos__ = cos_*-1
		
		cos1 = np.concatenate((cos_,cos__))
	#	print '2'
		
		cos_h = np.sort(cos1)

		ecdf = ECDF(cos_h) # Empirical cumulated distribution function
		y1 = ecdf(cos_h)-(cos_h+1.)/2. 
		
		x1 = np.array(cos_h)

		coeffs, cov = curve_fit(func, x1, y1)
		yfit = func(x1,coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4])
		d_yfit = dfunc(x1,coeffs[1],coeffs[2],coeffs[3],coeffs[4])

		a1_ran = coeffs[1]
		a2_ran = coeffs[2]
		a3_ran = coeffs[3]
		a4_ran = coeffs[4]

		#plot.step(x1,y1,c='#7f7f7f',alpha=0.2,lw=0.5)
		
		# Para calcular errores
	
		#p1, p2, p3 = np.percentile(x1,[25,50,75],interpolation='nearest')
		#y1_slice1.append( y1[np.where(x1==p1)[0][0]] )
		#y1_slice2.append( y1[np.where(x1==p2)[0][0]] )
		#y1_slice3.append( y1[np.where(x1==p3)[0][0]] )

		y1_col = Table.Column(y1,'yran_{}'.format(i))
		y1_slices.add_column(y1_col)
		
	#------------------------------------------------------------------------
	if nv==0:
		plt.legend(loc=0,fontsize=fs)
				
	plt.ylim(-.06,.05)

	#-----------------------------------------------------------------------
	"""Plot de errores"""
	
	nerr = 100 # Number of error bars
	percs = np.percentile(x1,np.linspace(0.,100.,nerr),interpolation='nearest')
	
	uperr = []
	loerr = []
	#errors = []
	for i in range(nerr):
		idx = np.where(x1==percs[i])[0][0]
		y1_slice = []
		for col in y1_slices.colnames:
			y1_slice.append( y1_slices[col][idx] )
		uperr.append( np.percentile(y1_slice,84,interpolation='nearest') )
		loerr.append( np.percentile(y1_slice,16,interpolation='nearest') )
		#errors.append( np.percentile(y1_slice,[16,84],interpolation='nearest') )
	
	xerr = 	np.linspace(-1.,1.,nerr)
	plt.fill_between(xerr,loerr,uperr,alpha=.5)

	asigma = []
	for i in range(nerr):
		asigma.append( np.mean((uperr[i],abs(loerr[i]))) )
	sigma = np.array(asigma)

	sigmaf = scipy.interpolate.interp1d(xerr,sigma,kind='cubic')
	#sigmafit = sigmaf(xerr) 
	#plot.plot(xerr,sigma,'--')

	#------------------------------------------------------------------------
	"""Tabla"""	
	row_labels=['N of low spin Gxs','N of high spin Gxs','N of total Gxs','Void radius']
	#row_labels=['row1']
	table_vals=[[nlow],[nhigh],[len(gals)],['{:.3}Mpc'.format(r/1000.)]]
	# the rectangle is where I want to place the table
	the_table = plt.table(cellText=table_vals,
				      colWidths = [.19],
				      rowLabels=row_labels,
				      loc=4)
	the_table.set_zorder(0.5)
	the_table.auto_set_font_size(False)
	the_table.set_fontsize(fs-1)
	the_table.scale(1,2)
	#------------------------------------------------------------------------
	
	plt.tight_layout()
	plt.savefig('BAAA_Void{}.pdf'.format(nv))
	plt.savefig('BAAA_Void{}.png'.format(nv))
	#-------------
	if nv==0:
		"""Plot JvsMass"""

		plt.figure(2,dpi=100,figsize=figsize)
		
		plt.tick_params(labelsize=fs)
		#plt.box('on')


		plt.xlabel('$log_{10}(M/M_{\odot})$',fontsize=fs)
		plt.ylabel(r'$log_{10}(\vec{J})\,\, [(kpc/h)(km/s)]$',fontsize=fs)

		plt.scatter(M,S,c=seaborn.color_palette()[0],label='High Spin Galaxies')
		plt.scatter(np.log10(gals_h['mass']),np.log10(gals_h['sp_n']),c=seaborn.color_palette()[1],label= 'Low Spin Galaxies')

		plt.plot(M,M*m+b,ls=':',label='Shell Galaxies Spin-Mass Linear Regression')
		plt.plot(M,M*g_m+g_b,c=seaborn.color_palette()[2],ls='--',label = 'Box Galaxies Spin-Mass Linear Regression')
		plt.legend(loc=2,fontsize=16)
		plt.tight_layout()
		plt.savefig('BAAA_LRPlot.pdf')
		plt.show()

	nv+=1
	
	plt.show()
plt.close()

