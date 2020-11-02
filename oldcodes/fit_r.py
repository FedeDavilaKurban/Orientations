"""
Plots de como varian el ajuste (funcion con parametros a1-4) con el radio r
"""
#rmin, rmax = .7, .9
perc = 75
nv = 14
#-------------------

print 'Reading data...'
#voidsProfile = rdr.read('voidsProfile.txt')#,format='commented_header')
#voidsFinal = rdr.read('voidsFinal.txt')
#subgroups = il.groupcat.loadSubhalos(basePath,135)

rho_m = 1.8284091259259258e-08

lbox = 75000. #kpc
bounds = [lbox,lbox,lbox]
#tree = PeriodicCKDTree(bounds,subgroups['SubhaloPos'])

#gxs = ascii.read('gxs.txt',format='commented_header')

#gxs = Table(subgroups['SubhaloPos'],names=['x','y','z'],dtype=['float64','float64','float64'])
#col_id = Table.Column([i for i in range(len(gxs))],name='ID')
#col_spx = Table.Column(subgroups['SubhaloSpin'][:,0],name='sp_x',dtype='float64')
#col_spy = Table.Column(subgroups['SubhaloSpin'][:,1],name='sp_y',dtype='float64')
#col_spz = Table.Column(subgroups['SubhaloSpin'][:,2],name='sp_z',dtype='float64')
#gxs.add_columns([col_spx,col_spy,col_spz,col_id])


yfit = []
xfit = []

def func(x,a,b,c,d,e):
	return a + b*np.sin( np.pi*(x+1.)/2. ) + c*np.sin( 2.*np.pi*(x+1.)/2. ) + d*np.sin( 3.*np.pi*(x+1.)/2. ) + e*np.sin( 4.*np.pi*(x+1.)/2. )

	
for rmin,rmax in [(.7,.9),(.8,1.),(.9,1.1),(1.,1.2)]:

	x,y,z,r = voidsFinal[nv]
	print rmin
	idx1 = tree.query_ball_point([x,y,z],r*rmax)
	idx2 = tree.query_ball_point([x,y,z],r*rmin)
	shell = [g for g in idx1 if g not in idx2]

	gals=gxs[shell]

	#----------------------------------------------------
	n=[]
	for i in range(len(gals)):
		n.append( np.linalg.norm( [gals[i]['sp_x'],gals[i]['sp_y'],gals[i]['sp_z']] ) )
	
	p1 = np.percentile(np.log(n),perc)
	
	g_highsp = gals[np.where(np.log(n)>p1)]
	
	cos_h=[]
	for i in range(len(g_highsp)):
		u = [g_highsp[i]['x']-x,g_highsp[i]['y']-y,g_highsp[i]['z']-z]
		v = [g_highsp[i]['sp_x'],g_highsp[i]['sp_y'],g_highsp[i]['sp_z']]
		cos_h.append( np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)) )
	
	cos_h.sort()

	ecdf = ECDF(cos_h) # Empirical cumulated distribution function

	y1 = ecdf(cos_h)-(np.array(cos_h)+1.)/2. # Diferencia entre la acumulada y la recta 

	xfit.append( np.array(cos_h) )

	coeffs, cov = curve_fit(func, xfit[-1], y1)
	yfit.append( func(xfit[-1],coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4]) )

	#a1 = coeffs[1]
	#a2 = coeffs[2]
	#print a1,a2
#------------------------------------------------------------------------
gs = gridspec.GridSpec(2, 1)

fig = plt.figure()
fig.suptitle('Void {} - VoidRadius = {}'.format(nv,int(voidsFinal[nv]['r'])),size=10)

ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[1,0])

ax1.plot(xfit[0],yfit[0],color='MidnightBlue',ls='-',label='(0.7-0.9)rVoid')
ax1.plot(xfit[1],yfit[1],color='DarkBlue',ls='--',label='(0.8-1.0)rVoid')
ax1.plot(xfit[2],yfit[2],color='Blue',ls='-.',label='(0.9-1.1)rVoid')
ax1.plot(xfit[3],yfit[3],color='DarkCyan',ls=':',label='(1.0-1.2)rVoid')
ax1.legend()

ax2.plot(voidsProfile['r'],(voidsProfile['Void{}'.format(nv)])/rho_m-1.,color='DarkRed')

ax1.set_title('Fits',size=14)
ax2.set_title('Density Profile',size=14)

plt.show()
#-------------------------------------------------

