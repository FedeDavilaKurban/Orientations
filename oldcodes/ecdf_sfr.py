# Correr primero coeffs_sfr.py

sfr_hsfo = np.sort(gals_fo_hs['SFR'])
sfr_hseo = np.sort(gals_eo_hs['SFR'])
sfr_lsfo = np.sort(gals_fo_ls['SFR'])
sfr_lseo = np.sort(gals_eo_ls['SFR'])

sfr_hsfo = sfr_hsfo[(sfr_hsfo)>10**-3.5]
sfr_hseo = sfr_hseo[(sfr_hseo)>10**-3.5]
sfr_lsfo = sfr_lsfo[(sfr_lsfo)>10**-3.5]
sfr_lseo = sfr_lseo[(sfr_lseo)>10**-3.5]

y1=ECDF(sfr_hsfo)
y2=ECDF(sfr_hseo)
y3=ECDF(sfr_lsfo)
y4=ECDF(sfr_lseo)

plt.step(sfr_hsfo,y1(sfr_hsfo),label='HighSpin FaceOn')
plt.step(sfr_hseo,y2(sfr_hseo),ls='--',label='HighSpin EdgeOn')
plt.step(sfr_lsfo,y3(sfr_lsfo),label='LowSpin FaceOn')
plt.step(sfr_lseo,y4(sfr_lseo),ls='--',label='LowSpin EdgeOn')

plt.xscale('log')
plt.xlabel('SFR')
plt.ylabel('ECDF')
plt.legend(loc=4)
plt.show()

