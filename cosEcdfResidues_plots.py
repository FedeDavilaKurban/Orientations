#Plots
fig, (ax1,ax2,ax3) = plt.subplots(3, 1,figsize=(5,6))
ax1.set_title('Sec={}'.format(sec))

ax1.hist(cos,bins=40,density=True,histtype='stepfilled',cumulative=False,alpha=1)#,alpha=30)
ax1.plot(cos,d_yfit+.5,'r--')

ax2.plot(cos,ecdf(cos))
ax2.plot(cos,(cos+1.)/2.,'b--')
ax2.plot(cos,yfit+(cos+1.)/2.,'r--')

ax3.step(cos,y)
ax3.plot(cos,yfit,'r--')
ax3.fill_between(xmean,ymean-ystd,ymean+ystd,color='k',alpha=.6,label=r'$1\sigma$')
ax3.fill_between(xmean,ymean-2*ystd,ymean+2*ystd,color='k',alpha=.4,label=r'$2\sigma$')
ax3.legend()

#ax3.text(-.75,.0,'a2=%.3f'%(a2))

# ax1.set_title('Cos (Mirrored)',size=12)
# ax2.set_title('ECDF',size=12)
# ax3.set_title('Residues',size=12)

plt.tight_layout()
plt.savefig('../plots/cosines_{}.png'.format(sec),dpi=300)
plt.show()
