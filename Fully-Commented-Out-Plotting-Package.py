### Making plots of data from ellcgatspyOpSim simulation -> output_file.csv ###
###In main code: ofile = 'output_file.csv' Line 684 ###

#importing needed libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pylab
from scipy.optimize import curve_fit


#pulling data file output_file.csv from running OpSim and inputting data into program
dat = pd.read_csv('output_file1.csv', sep=",", header=0)

#print(dat)

#These variables are for the histograms of Period/eccentricity vs. number
#f,ax = plt.subplots(figsize =(25,15))
PeriodIn = dat["p"]
PeriodOut = dat["LSM_PERIOD"]
ecc = dat['e']
#Need these for a mass-radius scatter plot
radius_1 = dat["r1"]
radius_2 = dat["r2"]
mass_1 = dat["m1"]
mass_2 = dat["m2"]

# #Making first plot

# ### PERIOD ###

f,ax = plt.subplots(figsize =(15,15))
#Sets up format of histogram
ax.hist(np.log10(PeriodIn), bins = 200, range = (0,10))
#Give hsitogram title
ax.set_title('Histogram of period')
#Setting axes labels
ax.set_xlabel('log10(Input Period)')
ax.set_ylabel('Number of binaries Identified, N')
#Displays histogram
pylab.show()
#Saves Plot as a pdf
f.savefig('periodin_histogram.pdf')

# f,ax = plt.subplots(figsize =(15,15))
# #Sets up format of histogram
ax.hist(np.log10(PeriodOut), bins = 200, range = (0,10))
# #Give hsitogram title
ax.set_title('Histogram of period')
# #Setting axes labels
ax.set_xlabel('log10(Output Period)')
ax.set_ylabel('Number of binaries Identified, N')
# #Displays histogram
# #pylab.show()
# #Saves Plot as a pdf
f.savefig('periodout_histogram.pdf')



# ##Ok, now onto the next thing

# ### ECCENTRICITY ###
f,ax = plt.subplots(figsize = (15,15))
pylab.hist(ecc, bins = 100, range =(0,1))
pylab.title('Distribution of eccentricities')
ax.set_xlabel('Eccentricity')
ax.set_ylabel('Number of binaries, N')
# #pylab.show()
f.savefig('eccentricity_histogram.pdf')

# ##### Scatter plot of mass_1 to radius_1 #####

# #Since this is a scatterplot, can just use .plot() function instead of .hist()
# f,ax = plt.subplots(figsize = (15,15))
# ax.plot(radius_1, mass_1, 'o', markersize = 5)
# pylab.title('Mass 1 versus Radius 1')
# ax.set_xlabel('Radius 1')
# ax.set_ylabel('Mass 1')
# # #pylab.show()
# f.savefig('radius_1-mass_1-scatter-plot.pdf')

# #Contour map for mass1 and radius 1:
# # from numpy import linspace, meshgrid
# # from matplotlib.mlab import griddata

# # def grid(radius_1, mass_1, resX=100, resY=100):
# #     xi = linspace(min(radius_1), max(radius_1), resX)
# #     yi = linspace(min(mass_1), max(mass_1), resY)

# #     Z = griddata(radius_1, mass_1, xi, yi)
# #     X, Y = meshgrid(xi, yi)
# #     return X, Y, Z
# # X, Y, Z = grid(radius_1, mass_1)
# # plt.contour(X, Y, Z)
# # pylab.show()

# # ### Same thing, but with mass_2 and radius_2 ###
# f,ax = plt.subplots(figsize = (15,15))
# ax.plot(radius_2, mass_2, 'o', markersize = 5)
# pylab.title('Mass 2 versus Radius 2')
# ax.set_xlabel('Radius 2')
# ax.set_ylabel('Mass 2')
# # # #pylab.show()
# f.savefig('radius_2-mass_2-scatter-plot.pdf')

####remember to save plots for future use

h2D = plt.subplots(figsize = (10,10))
plt.hist2d(radius_1, mass_1, bins = [50,50], range = [[0,1],[0,1]], cmap = 'binary')
plt.show()
h2D.savefig('radius_1-mass_1-heatmap.pdf')

#mass2 and radius 2 heat map
h2D = plt.subplots(figsize = (10,10))
plt.hist2d(radius_2, mass_2, bins = [50,50], range = [[0,1],[0,1]], cmap = 'binary')
plt.show()
h2D.savefig('radius_2-mass_2-heatmap.pdf')