# Importing needed libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance

# dat2 = pd.read_csv(filename, sep = ',', header=2)

# Trying to read in all files at once
import glob
path = '/Users/andrewbowen/BO_EB-LSST/code/Opsim'
# path ='atb0836@quest.it.northwestern.edu:/projects/p30137/ageller/testing/EBLSST/fix_omega/output_files' # use your path
allFiles = glob.glob(path + "/*.csv")
frame = pd.DataFrame()
list_1 = []
list_2 = []
for file_ in allFiles:
    dat1 = pd.read_csv(file_, sep = ',', header=0, nrows = 1)
    dat2 = pd.read_csv(file_, sep = ',', header=2)
    dat2['RA'] = dat1['OpSimRA'][0]
    dat2['Dec'] = dat1['OpSimDec'][0]
    list_1.append(dat1)
    list_2.append(dat2)
frame1 = pd.concat(list_1) #RA, Dec and id stuff
frame2 = pd.concat(list_2) #Actual data

#print(frame2)

PeriodIn = frame2['p'] # input period -- 'p' in data file
PeriodOut = frame2['LSM_PERIOD'] #LSM_PERIOD in data file
Mass1 = frame2['m1']
Mass2 = frame2['m2']
radius1 = frame2['r1']
radius2 = frame2['r2']

RA = frame1['OpSimRA']
Dec = frame1['OpSimDec']
OpSimCoords = SkyCoord(RA,Dec,unit=(u.degree, u.degree),frame='icrs')

# u filter mollweide
plt.figure()
plt.subplot( projection = 'mollweide' )
plt.grid(True)
plt.scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
            cmap = 'Blues', c = frame1['NOpSimObs_u'], vmin = 0, vmax = 200)
clb = plt.colorbar()
clb.set_label(r'Number of Observations', rotation = 270)
plt.xlabel(r'$\alpha$', fontsize = 16)
plt.ylabel(r'$\delta$', fontsize = 16)


plt.title('u Coordinates', fontsize = 16)
plt.savefig('u_mollweide.pdf')

# g filter mollweide
plt.figure()
plt.subplot( projection = 'mollweide' )
plt.grid(True)

plt.scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian,s = 5,\
            cmap = 'Blues', c = frame1['NOpSimObs_g'], vmin = 0, vmax = 200)
clb = plt.colorbar()
clb.set_label(r'Number of Observations', rotation = 270)
plt.xlabel(r'$\alpha$', fontsize = 16)
plt.ylabel(r'$\delta$', fontsize = 16)
plt.title('g Coordinates', fontsize = 16)

# r filter mollweide
plt.figure()
plt.subplot( projection = 'mollweide' )
plt.grid(True)
plt.scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
            cmap = 'Blues', c = frame1['NOpSimObs_r'], vmin = 0, vmax = 200)
clb = plt.colorbar()
clb.set_label(r'Number of Observations', rotation = 270)
plt.xlabel(r'$\alpha$', fontsize = 16)
plt.ylabel(r'$\delta$', fontsize = 16)

plt.title('r Coordinates', fontsize = 16)

# i filter mollweide
plt.figure()
plt.subplot( projection = 'mollweide' )
plt.grid(True)

plt.scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
            cmap = 'Blues', c = frame1['NOpSimObs_i'], vmin = 0, vmax = 200)
clb = plt.colorbar()
clb.set_label(r'Number of Observations', rotation = 270)
plt.xlabel(r'$\alpha$', fontsize = 16)
plt.ylabel(r'$\delta$', fontsize = 16)

plt.title('i Coordinates', fontsize = 16)

# z filter mollweide
plt.figure()
plt.subplot( projection = 'mollweide' )
plt.grid(True)

plt.scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian,s = 5,\
            cmap = 'Blues', c = frame1['NOpSimObs_z'], vmin = 0, vmax = 200)
clb = plt.colorbar()
clb.set_label(r'Number of Observations', rotation = 270)
plt.xlabel(r'$\alpha$', fontsize = 16)
plt.ylabel(r'$\delta$', fontsize = 16)

plt.title('z Coordinates', fontsize = 16)

# y filter mollweide
plt.figure()
plt.subplot( projection = 'mollweide' )
plt.grid(True)

plt.scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s = 5,\
            cmap = 'Blues', c = frame1['NOpSimObs_y'], vmin = 0, vmax = 200)
clb = plt.colorbar()
clb.set_label(r'Number of Observations', rotation = 270)
plt.xlabel(r'$\alpha$', fontsize = 16)
plt.ylabel(r'$\delta$', fontsize = 16)

plt.title('y Coordinates', fontsize = 16)

# Period Recovery 
Perr1 = PeriodIn - (0.1*PeriodIn)
Perr2 = PeriodIn + (0.1*PeriodIn)
pdiff = abs(PeriodIn - PeriodOut)/PeriodIn
P_Recover = frame2.loc[pdiff < 0.1]

# ID/All Binaries

Total_Recovered = P_Recover['p'].count
Total_Periods = frame2['p'].count
# print(Total_Recovered)
ID_rate = P_Recover.shape[0]/frame2.shape[0]

# Half-Period correction
Half_Period = 0.5 * PeriodIn
LSM_Half_Period = 0.5 * PeriodIn
hpdiff = abs(Half_Period - LSM_Half_Period)/Half_Period
# HP_Recover = frame2.loc[hpdiff < 0.1]
pdiff = abs(PeriodIn.values - PeriodOut.values)/PeriodIn.values
P_Recover = frame2.loc[(pdiff < 0.1) | (hpdiff < 0.1)]

# Detected Binaries
detected_bins = frame2.loc[PeriodOut != -999]
# print(detected_bins.shape[0])
detection_rate = detected_bins.shape[0]/frame2.shape[0]
# print(detection_rate)

# ID/Detected
finder = P_Recover.shape[0]/detected_bins.shape[0]

### Defining Variables for grid plot  ###
# Period
p = frame2['p']
d_p = detected_bins['p']
i_p = P_Recover['p']

# Eccentricity
ecc = frame2['e']
d_ecc = detected_bins['e']
i_ecc = P_Recover['e']

# inclination
inc = frame2['i']
d_inc = detected_bins['i']
i_inc = P_Recover['i']

# Mass
mass1= frame2['m1']
mass2 = frame2['m2']
mass_ratio = mass2/mass1

i_mass1 = P_Recover['m1']
i_mass2 = P_Recover['m2']
i_mass_ratio = i_mass2/i_mass1
d_mass1 = detected_bins['m1']
d_mass2 = detected_bins['m2']
d_mass_ratio = d_mass2/d_mass1

# Radius
radius1 = frame2['r1']
radius2 = frame2['r2']
radius_ratio = radius2/radius1
d_radius = detected_bins['r2']/detected_bins['r1']
i_radius = P_Recover['r2']/P_Recover['r1']

################################################################################################################

IDdict = {'RA':np.array([]), 'Dec':np.array([]), 'frac':np.array([])}

for file in allFiles:
    dat1 = pd.read_csv(file, sep = ',', header=0, nrows = 1)
    dat2 = pd.read_csv(file, sep = ',', header=2)
#     appending RA/Dec to dictionary
    if dat2['p'].values[0] != -1:
        IDdict['RA'] = np.append(IDdict['RA'],dat1['OpSimRA'])
        IDdict['Dec'] = np.append(IDdict['Dec'], dat1['OpSimDec'])

    #     Criteria for period recovery
        PeriodIn = dat2['p']
        PeriodOut = dat2['LSM_PERIOD']
        pdiff = abs(PeriodIn.values - PeriodOut.values)/PeriodIn.values
        f = 0

        recovered_p = dat2['LSM_PERIOD'].loc[pdiff<0.1]
        all_p = dat2['LSM_PERIOD']
    #     Putting recovered period into 'frac' lists
        f = recovered_p.shape[0]/all_p.shape[0]

    
        IDdict['frac'] = np.append(IDdict['frac'],f)


# Mollweide plot of percent recovered at each RA/Dec
RA = IDdict['RA']
Dec = IDdict['Dec']
Coords = SkyCoord(RA,Dec,unit=(u.degree, u.degree),frame='icrs')

plt.figure()
plt.subplot( projection = 'mollweide' )
plt.grid(True)

plt.scatter(Coords.ra.wrap_at(180.*u.degree).radian,Coords.dec.radian, s = 20,\
           cmap = 'Blues', c = IDdict['frac']*100, vmin = 0, vmax = 0.2)

plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

clb = plt.colorbar()
clb.set_label(r'% Periods Recovered', fontsize = 18, rotation = 270, labelpad = 20)
plt.xlabel(r'$\alpha$', fontsize = 16)
plt.ylabel(r'$\delta$', fontsize = 16) 
plt.savefig('mollweide_percent_rec.pdf')

##################################################################################################################

# Big Grid Plot

f, axarr = plt.subplots(4,5, figsize = (18,12), sharex = 'col')
f.subplots_adjust(wspace=0.3, hspace = 0)

# column titles
axarr[0,0].set_title('Period', fontsize = 24)
axarr[0,1].set_title('$M_1$/$M_2$', fontsize = 24)
axarr[0,2].set_title('Eccentricity', fontsize = 24)
axarr[0,3].set_title('$R_1$/$R_2$', fontsize = 24)
axarr[0,4].set_title('Inclination', fontsize = 24)
axarr[0,0].set_ylabel('All Binaries', fontsize = 18)
axarr[1,0].set_ylabel('Observable EBs', fontsize = 18)
axarr[2,0].set_ylabel('Recoverable EBs', fontsize = 18)
axarr[3,0].set_ylabel('Cumulative', fontsize = 18)

# Period
axarr[0,0].hist(np.log10(p), bins = 50, range = (0,10), color = '#13294B')
axarr[1,0].hist(np.log10(d_p), bins = 50, range = (0,10), color = '#356897')
axarr[2,0].hist(np.log10(i_p), bins = 50, range = (0,10), color = '#99badd')
axarr[3,0].hist(np.log10(p), bins = 1000, range = (0,10), color = '#13294B',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,0].hist(np.log10(d_p), bins = 1000, range = (0,10), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,0].hist(np.log10(i_p), bins = 1000, range = (0,10), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,0].set_xlabel('Days')

# Mass Ratio
axarr[0,1].hist(mass_ratio, bins = 50, range = (0,2), color = '#13294B')
axarr[1,1].hist(d_mass_ratio, bins = 50, range = (0,2), color = '#356897')
axarr[2,1].hist(i_mass_ratio, bins = 50, range = (0,2), color = '#99badd')
axarr[3,1].hist(mass_ratio, bins = 1000, range = (0,2), color = '#13294B',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,1].hist(d_mass_ratio, bins = 1000, range = (0,2), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,1].hist(i_mass_ratio, bins = 1000, range = (0,2), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)


# Eccentricity
axarr[0,2].hist(ecc, bins = 50, range = (0,1), color = '#13294B')
axarr[1,2].hist(d_ecc, bins = 50, range = (0,1), color = '#356897')
axarr[2,2].hist(i_ecc, bins = 50, range = (0,1), color = '#99badd')
axarr[3,2].hist(ecc, bins = 1000, range = (0,1), color = '#13294B',density = True,\
             histtype = 'step', fill = False,cumulative = True)
axarr[3,2].hist(d_ecc, bins = 1000, range = (0,1), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,2].hist(i_ecc, bins = 1000, range = (0,1), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)

# Radius Ratio
axarr[0,3].hist(radius_ratio, bins = 50, range = (0,2), color = '#13294B')
axarr[1,3].hist(d_radius, bins = 50, range = (0,2), color = '#356897')
axarr[2,3].hist(i_radius, bins = 50, range = (0,2), color = '#99badd')
axarr[3,3].hist(radius_ratio, bins = 1000, range = (0,2), color = '#13294B',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,3].hist(d_radius, bins = 1000, range = (0,2), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,3].hist(i_radius, bins = 1000, range = (0,2), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)


#Inclination 
axarr[0,4].hist(inc, bins = 50, range = (0,90), color = '#13294B')
axarr[1,4].hist(d_inc, bins = 50, range = (0,90), color = '#356897')
axarr[2,4].hist(i_inc, bins = 50, range = (0,90), color = '#99badd')
axarr[3,4].hist(inc, bins = 1000, range = (0,90), color = '#13294B',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,4].hist(d_inc, bins = 1000, range = (0,90), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,4].hist(i_inc, bins = 1000, range = (0,90), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,4].set_xlabel('Degrees')

f.savefig('Poster_grid_plot.pdf')
