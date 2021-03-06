"""
El nombre 'spin.py' tiene muy poco que ver con lo que hace el programa
"""

#converters = {'x': [ascii.convert_numpy(np.float64)],
#	      'y': [ascii.convert_numpy(np.float64)],
#	      'z': [ascii.convert_numpy(np.float64)],
#	      'r': [ascii.convert_numpy(np.float64)]}


#print 'Reading Voids...'
#voidsFinal = rdr.read('voidsFinal.txt')#,converters=converters)
#print 'Reading Data...'
#subgroups = il.groupcat.loadSubhalos(basePath,135)

#---------------------------------------------------------

r = voidsFinal[nv]['r']
a1_ran = []
a2_ran = []
#rmin = .7
#rmax = 1.2
#perc = 75
#---------------------------------------------------------


gxs = Table(subgroups['SubhaloPos'],names=['x','y','z'],dtype=['float64','float64','float64'])
col_id = Table.Column([i for i in range(len(gxs))],name='ID')
col_sp = Table.Column(subgroups['SubhaloSpin'],name='sp',dtype='float64')
col_m = Table.Column(subgroups['SubhaloMass'],name='mass',dtype='float64')

gxs.add_columns([col_sp,col_m,col_id])
gxs = gxs[(np.log10(gxs['mass'])<3.)&(np.log10(gxs['mass'])>-1.)]

col_spn = Table.Column(np.sum(np.abs(gxs['sp'])**2,axis=-1)**(1./2),name='sp_n',dtype='float64')
gxs.add_column(col_spn,4)

print 'Selecting gxs with high spin...'

# Nueva forma de definir gxs (cortadas en masa) con spin alto o bajo
m,b = 0.65202943253484447, 1.840361622864827 # Sacados con stats.linregress() de todas las gxs con 0.<log10(Mass)<3. en un grafico J vs Mass

g_highsp = gxs[np.log10(gxs['sp_n'])>(np.log10(gxs['mass'])*m+b)] # Gxs de alto spin siguiendo la relacion J 

#---------------------------------------------------------

#lbox = 75000. #kpc
#bounds = [lbox,lbox,lbox]
tree = spatial.cKDTree(data=np.column_stack((g_highsp['x'],g_highsp['y'],g_highsp['z'])),boxsize=lbox)
#tree = spatial.cKDTree(data=np.column_stack((gals_h['x'],gals_h['y'],gals_h['z'])),boxsize=lbox)
#---------------------------------------------------------

nv = 0 # Void ID
#print 'Calculating Coefficients [a1, a2] of void',nv,'...'
#execfile('orientation_fit.py')

#a1 = coeffs[1]
#a2 = coeffs[2]

#print a1,a2
#---------------------------------------------------------
def func(x,a,b,c,d,e):
	return a + b*np.sin( np.pi*(x+1.)/2. ) + c*np.sin( 2.*np.pi*(x+1.)/2. ) + d*np.sin( 3.*np.pi*(x+1.)/2. ) + e*np.sin( 4.*np.pi*(x+1.)/2. )
# Esta definida en orientation_fit.py
#---------------------------------------------------------
#a1 = 0.0179067226379
#a2 = -0.0122519325717

for j in range(500):

	print j+1
		
	p = (75000.)*np.random.rand(1,3)
	idx1 = tree.query_ball_point(p,r*rmax)
	idx2 = tree.query_ball_point(p,r*rmin)
	shell = [g for g in idx1 if g not in idx2]

	gals=g_highsp[shell]
	#gals = gals_h[shell]
	cos_h=[]
	
	for i in range(len(gals)):
		u = [gals[i]['x']-p[0][0],gals[i]['y']-p[0][1],gals[i]['z']-p[0][2]]
		v = [gals[i]['sp'][0],gals[i]['sp'][1],gals[i]['sp'][2]]
		cos_h.append( abs(np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v))) )
	
	cos_h.sort()

	ecdf = ECDF(cos_h) # Empirical cumulated distribution function

	y = ecdf(cos_h)-np.array(cos_h) # Diferencia entre la acumulada y la recta 

	#plt.step(cos_h,y)
	#plt.show()


	x = np.array(cos_h)

	coeffs, cov = curve_fit(func, x, y)
	yfit=func(x,coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4])

	a1_ran.append(coeffs[1])
	a2_ran.append(coeffs[2])


#pvalue1 =  float(len([x for x in a1_ran if x>a1])) / float(len(a1_ran)) 
#pvalue2 =  float(len([x for x in a2_ran if x<a2])) / float(len(a2_ran)) 

#print pvalue1, pvalue2

#fig = plt.figure()
#ax1 = fig.add_subplot(1, 2, 1)
#ax2 = fig.add_subplot(1, 2, 2)

#fig.suptitle('Rmin, Rmax: '+str(rmin)+', '+str(rmax)+'. Percentile: '+str(perc) )

#ax1.set_title('$a1$ Coefficient',size=18)
#ax2.set_title('$a2$ Coefficient',size=18)

#ax1.text(left,top,'pvalue:'+str(pvalue1))#,ha='center', va='center',transform=ax.transAxes)

#ax1.annotate('P value: '+str(round(pvalue1,2)), xy=(0.05, 0.95), xycoords='axes fraction')
#ax2.annotate('P value: '+str(round(pvalue2,2)), xy=(0.05, 0.95), xycoords='axes fraction')

#n1, bins, patches = ax1.hist(a1_ran,normed=True,histtype='stepfilled')
#n2, bins2, patches2 = ax2.hist(a2_ran,normed=True,histtype='stepfilled')

#ax1.axvline(a1,color='r',linestyle='dashed',label='a1')
#ax2.axvline(a2,color='r',linestyle='dashed',label='a2')

#ax1.legend()
#ax2.legend()
#plt.show()


	
