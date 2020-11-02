x,y,z,r = voidsFinal[nv]
print tree
idx1 = tree.query_ball_point([x,y,z],r*1.2)
idx2 = tree.query_ball_point([x,y,z],r*.7)
shell = [g for g in idx1 if g not in idx2]

g_highsp=g_highsp[shell]
#----------------------------------------------------

#g_highsp = gals

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

#-----------------------------------------------------------
cos_h.sort()

ecdf = ECDF(cos_h) # Empirical cumulated distribution function

y = ecdf(cos_h)-(np.array(cos_h)+1.)/2. # Diferencia entre la acumulada y la recta 

plt.step(cos_h,y)
plt.show()

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

	
