a2_lowMass = []
a2_midMass = []
a2_highMass = []

for i in range (22):
	a2_lowMass.append( paramList[i]['a2'][0] )
	a2_midMass.append( paramList[i]['a2'][1] )
	a2_highMass.append( paramList[i]['a2'][2] )
     
s1,p1 = scipy.stats.ttest_ind(a2_ran,a2_lowMass,equal_var=False)
s2,p2 = scipy.stats.ttest_ind(a2_ran,a2_midMass,equal_var=False)
s3,p3 = scipy.stats.ttest_ind(a2_ran,a2_highMass,equal_var=False)

fig = plt.figure()
ax1 = fig.add_subplot(3, 1, 1)
ax2 = fig.add_subplot(3, 1, 2)
ax3 = fig.add_subplot(3, 1, 3)

ax1.set_title('Low Mass')
ax2.set_title('Mid Mass')
ax3.set_title('High Mass')

ax1.annotate('P value: '+'%.3f'%p1, xy=(0.05, 0.95), xycoords='axes fraction')
ax2.annotate('P value: '+'%.3f'%p2, xy=(0.05, 0.95), xycoords='axes fraction')
ax3.annotate('P value: '+'%.3f'%p3, xy=(0.05, 0.95), xycoords='axes fraction')

n1, bins, patches = ax1.hist(a2_ran,bins=15,normed=True,histtype='stepfilled')
n2, bins2, patches2 = ax2.hist(a2_ran,bins=15,normed=True,histtype='stepfilled')
n3, bins3, patches3 = ax3.hist(a2_ran,bins=15,normed=True,histtype='stepfilled')

for i in range(22):
	ax1.axvline(a2_lowMass[i],color='r',linestyle='dotted')
	ax2.axvline(a2_midMass[i],color='r',linestyle='dotted')
	ax3.axvline(a2_highMass[i],color='r',linestyle='dotted')

ax1.axvline(np.mean(a2_lowMass),color='r')
ax2.axvline(np.mean(a2_midMass),color='r')
ax3.axvline(np.mean(a2_highMass),color='r')
ax1.axvline(np.mean(a2_ran),color='b')
ax2.axvline(np.mean(a2_ran),color='b')
ax3.axvline(np.mean(a2_ran),color='b')



ax1.legend()
ax2.legend()
ax3.legend()
plt.show()


#plt.boxplot([a2_ran,a2_lowMass,a2_midMass,a2_highMass])
#plt.show()
