"""
Programa para mostrar graficos de como el metodo utilizado refleja las orientaciones
Para la AAA2017

31/1/2018 "posteraaa2017_ejemplos_1.py" - Modificacion para re-añadir un tercer gráfico intermedio
"""

#---------------------------------
# Genero datos artificiales para plot del poster AAA2017
gauss = np.random.normal(0.,.6,100000)
gauss = gauss[abs(gauss)<1.]
gauss1 = gauss[(gauss>0.)]*-1.+1.
gauss2 = gauss[(gauss<0.)]*-1.-1.
gaussinv = np.concatenate((gauss1,gauss2)) 

unif = np.random.uniform(-1.,1.,100000)


#cos1 = gaussinv
#----
gs = gridspec.GridSpec(2, 1)

fig = plt.figure(figsize=(6*1.618,6))#,dpi=100)
#plt.box('on')

ax1 = fig.add_subplot(gs[0,0])
ax3 = fig.add_subplot(gs[1,0])
count=0
for cos1 in [gauss,gaussinv,unif]:
	count+=1
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

	#------------------------------------------------------------------------
	if count==1: c,alpha,ls=[seaborn.color_palette()[0],1.,'--']
	if count==2: c,alpha,ls=[seaborn.color_palette()[1],1.,'-.']
	if count==3: c,alpha,ls=[seaborn.color_palette()[2],1.,':']

	ax1.plot(x,d_yfit+.5,ls=ls,c=c)
	ax1.tick_params(labelsize=fs)
	#ax1.box('on')
	errN = ax1.hist(cos1,bins=40,normed=True,histtype='step',cumulative=False,linewidth=1.,alpha=.7)
	#ax1.set_xlabel(r'$cos(\theta)$',size=fs)
	#ax1.legend()

	ax3.plot(x,yfit,ls=ls,c=c)
	ax3.step(x,y,alpha=.7)
	ax3.set_xlabel(r'$cos(\theta)$',size=fs)
	ax3.set_ylabel('Residues',size=fs)
	ax3.tick_params(labelsize=fs)
	#ax3.box('on')

#ax1.set_title('Cosine',size=fs)
#ax2.set_title('ECDF',size=12)
#ax3.set_title('Residues',size=fs)

plt.tight_layout()
plt.savefig('BAAA_Examples.pdf')
plt.show()
