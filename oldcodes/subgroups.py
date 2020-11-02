
rmin, rmax = .7, 1.2
perc = 75
nv = 0 

#----------------------------------------------------
subgroups = il.groupcat.loadSubhalos(basePath,135)

gxs = Table(subgroups['SubhaloPos'],names=['x','y','z'],dtype=['float64','float64','float64'])
col_id = Table.Column([i for i in range(len(gxs))],name='ID')
col_spx = Table.Column(subgroups['SubhaloSpin'][:,0],name='sp_x',dtype='float64')
col_spy = Table.Column(subgroups['SubhaloSpin'][:,1],name='sp_y',dtype='float64')
col_spz = Table.Column(subgroups['SubhaloSpin'][:,2],name='sp_z',dtype='float64')
gxs.add_columns([col_spx,col_spy,col_spz,col_id])


converters = {'x': [ascii.convert_numpy(np.float64)],
	      'y': [ascii.convert_numpy(np.float64)],
	      'z': [ascii.convert_numpy(np.float64)],
	      'r': [ascii.convert_numpy(np.float64)]}


print 'Reading data...'
voidsFinal = ascii.read('voidsFinal.txt',converters=converters)

lbox = 75000. #kpc
bounds = [lbox,lbox,lbox]
tree = PeriodicCKDTree(bounds,subgroups['SubhaloPos'])
#tree = PeriodicCKDTree(bounds,np.column_stack((gxs['x'],gxs['y'],gxs['z'])))


#---------------------------------------------------------
#for nv in range(len(voidsFinal)):

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
#cos=[]
#for i in range(len(gals)):
#	u = [gals[i]['x']-x,gals[i]['y']-y,gals[i]['z']-z]
#	v = [gals[i]['sp_x'],gals[i]['sp_y'],gals[i]['sp_z']]
#	cos.append( np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)) )
	
#plt.close()
#plt.hist(cos,bins=3000,histtype='stepfilled',cumulative=True)
#plt.xlabel('$cos \\theta $')
#plt.show()



#----------------------------------------------------
n=[]
for i in range(len(gals)):
	n.append( np.linalg.norm( [gals[i]['sp_x'],gals[i]['sp_y'],gals[i]['sp_z']] ) )
	
p1 = np.percentile(np.log(n),perc)
#p2 = np.percentile(np.log(n),10)
	
g_highsp = gals[np.where(np.log(n)>p1)]
#g_lowsp = gals[np.where(np.log(n)<p2)]
	
#print len(g_highsp),len(g_lowsp)
cos_h=[]
for i in range(len(g_highsp)):
	u = [g_highsp[i]['x']-x,g_highsp[i]['y']-y,g_highsp[i]['z']-z]
	v = [g_highsp[i]['sp_x'],g_highsp[i]['sp_y'],g_highsp[i]['sp_z']]
	cos_h.append( np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)) )
	
#cos_l=[]
#for i in range(len(g_lowsp)):
#	u = [g_lowsp[i]['x']-x,g_lowsp[i]['y']-y,g_lowsp[i]['z']-z]
#	v = [g_lowsp[i]['sp_x'],g_lowsp[i]['sp_y'],g_lowsp[i]['sp_z']]
#	cos_l.append( np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)) )
	
#plt.close()

#fig = plt.figure()
#ax1 = fig.add_subplot(1, 2, 1)
#ax2 = fig.add_subplot(1, 2, 2)

#n1, bins, patches = ax1.hist(cos_h,bins=50,normed=True,histtype='stepfilled',cumulative=True,alpha=30.)
#ax1.set_xlabel('$cos \\theta (High Spin)$',size=18)


#n2, bins1, patches1 = ax2.hist(cos_l,bins=50,normed=True,histtype='stepfilled',cumulative=True,alpha=30.)
#ax2.set_xlabel('$cos \\theta (Low Spin)$',size=18)

#plt.show()

#-------------------------------------------------------------
#ran = [(1.-(-1.))*np.random.random()+(-1.) for i in range(len(g_highsp))] #(b - a) * random_sample() + a

#ran = np.linspace(-1,1,num=len(g_highsp))


#plt.close()

#fig = plt.figure()
#ax1 = fig.add_subplot(3, 1, 1)
#ax2 = fig.add_subplot(3, 1, 2)
#ax3 = fig.add_subplot(3, 1, 3)

#n1, bins, patches = ax1.hist(cos_h,bins=100,normed=True,histtype='stepfilled',cumulative=True,alpha=30.)
#ax1.set_title('$cos \\theta (High Spin)$',size=10)


#n2, bins2, patches2 = ax2.hist(ran,bins=100,normed=True,histtype='stepfilled',cumulative=True,alpha=30.)

#ax2.plot(ran,'r-')
#ax2.set_title('$cos \\theta (Random)$',size=10)

#ax3.plot((np.array(cos_h)-np.array(ran)))
#ax3.plot((n1-n2),'m:')
#ax3.set_title('$HighSpin - Random$',size=10)

#plt.show()

#plt.step(cos_h,(ecdf(cos_h)-(cos_h+1.)/2.))

#-----------------------------------------------------------






#-----------------------------------------------------------
cos_h.sort()

ecdf = ECDF(cos_h) # Empirical cumulated distribution function

y = ecdf(cos_h)-(np.array(cos_h)+1.)/2. # Diferencia entre la acumulada y la recta 

#plt.step(cos_h,y)
#plt.show()

def func(x,a,b,c,d,e):
	return a + b*np.sin( np.pi*(x+1.)/2. ) + c*np.sin( 2.*np.pi*(x+1.)/2. ) + d*np.sin( 3.*np.pi*(x+1.)/2. ) + e*np.sin( 4.*np.pi*(x+1.)/2. )

x = np.array(cos_h)

coeffs, cov = curve_fit(func, x, y)
yfit=func(x,coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4])

a1 = coeffs[1]
a2 = coeffs[2]
print a1,a2

plt.close()

plt.plot(x,y)
plt.plot(x,yfit,'r--')
plt.show()

#------------------------------------------------------------------------
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

n1, bins, patches = ax1.hist(a1_ran,normed=True,histtype='stepfilled')
n2, bins2, patches2 = ax2.hist(a2_ran,normed=True,histtype='stepfilled')

ax1.axvline(a1,color='r',linestyle='dashed',label='a1')
ax2.axvline(a2,color='r',linestyle='dashed',label='a2')

ax1.legend()
ax2.legend()
plt.show()




