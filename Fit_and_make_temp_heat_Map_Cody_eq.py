#Non-normlaized fit
# Fits 2 lorenzians to data and claculates temp with Cody equation, then makes a seaborn and custom heatmap. 





import numpy as np
import seaborn as sns
import os
import glob
import scipy as sp
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit
from scipy.integrate import simps
import sys

import csv

#declare heat map sort arrays


name = "Pre_damage_Raman_map_19"

cords_x = np.array([])
cords_y = np.array([])
temp_l = np.array([])
temp_i = np.array([])
carbon_amount = np.array([])

def get_x_y_from_filename(fname):
  fname_split = [
    item for item in fname.split('_')
    if item != ''
  ]
  x = fname_split[fname_split.index('X') + 1]
  y = fname_split[fname_split.index('Y') + 1]
  map_i = fname_split[fname_split.index('map')+1]
  return x, y, map_i


# make obj_array
datadir = '/Users/amylebleu-debartola/Desktop/Raman_Damage/Jbilet_Winselwan/focus/Map_19_undamaged/' #string to find position
AllFiles = np.array([])                       #declare empty array to be filled
for file in glob.glob(datadir+'*.txt'):       #start loop to go to the directory above, get all files with .txt
    AllFiles = np.append(AllFiles, file)      #add each file that has .txt to AllFiles array
#test to see if things got fed in properly
#print(AllFiles)

#declare golbal variables
k = len(AllFiles)  #size of AllFiles 
redFlag = 0        #error marker
t2=0               #time counter

#stuff here is to faciliate wrting out into a csv file
#file = open('Murchisan_test_map_python.txt', 'w')
data_fit_all = None

def plot_fitted_data(x, y, fitted, peak1, peak2,baseline,extra):
        """Plot data following curve fitting"""
        #plt.style.use('dark_background')
        plt.title("2 Lorentzian fit"+extra)
        plt.xlabel('Raman shift in cm-1')
        plt.ylabel('Normalizied counts')
        # Captures the most meaningful range of X to visualize plt.xlim((950, 2000))
               # Plot the original values
        
        plt.plot(x, y, '#ffcfea')  #prink
        # Plot the fitted curve
        plt.plot(x, fitted,'#6EB5FF' )  # blue
        plt.plot(x, peak1, '#b19cd9')  #lavandar
        plt.plot(x, peak2, '#afe9ff')  #light blue
        #plt.plot(x, baseline, '#bfb9ff')  # purple?
        plt.gcf()
        plt.savefig(extra+"_2lor_fit_non_norm.png") 
        plt.clf()
        #plt.show()

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

#for j in tqdm(range(k), desc="Processing data files...")::  #the range of 0 to the number of files
for j in range(k):
    t2=t2+1
    
    #read in txt file, convert 2d array to 2 1d arrays, X Y
    CurFileName =AllFiles[j]        # get name of 1 file and path 
    parse_one = CurFileName.split('_')
   # print(parse_one)
    x, y, map_i = get_x_y_from_filename(CurFileName)
    x_cords_list, y_cords_list = float(x),float(y)
    map_i = float(map_i)
    #print(x_cords_list)
    
    CurFileX = np.array([])         # make empty array for X
    CurFileY = np.array([])         # make empty array for Y
    f = open(CurFileName, 'r')      # open up file, read only(r) 
    lines = np.loadtxt(f,)          # load all data into array "lines"
    #print(lines)                   # check

    # Pull out section of data we are intersted in
    truncated_data = lines[170:399]  #range in X axis of 1000-2000: similar to range on origin_plot 
    #truncated_data = lines
    CurFileX, CurFileY = truncated_data[:,0], truncated_data[:,1]

    #find highest y, divide rest of the data by this to normilze(between 0-1)
    high_y = np.amax(CurFileY)
    

    #translate x and y into more easly workable variables, normlizae y
    x = CurFileX
    y = CurFileY

    #y = smooth(y,)
    #y[0:15] = y[16]
    #y[:205220] = y[203]

 
    ##playing with gaussian from here
    def lorinzian(x, amp, ctr, wid):
        """Gaussian for fitting"""
        # Utilize elementwise functions from numpy
        return (amp / (1+((x-ctr)/wid)**2))

    def double_lorinzian_extended(x, *params):
        """Combined double gaussian w/ polynomial etc"""
      # Destructure and assign parameters passed in
        const, y1, y2, amp_1, ctr_1, wid_1, amp_2, ctr_2, wid_2 = params
      # Calculate our double gaussian with additive extras
        return (lorinzian(x, amp_1, ctr_1, wid_1)
                + lorinzian(x, amp_2, ctr_2, wid_2)
                + const + (y1 * x) + (y2 * np.power(x, 2.)))



    # Set initial parameters and bounds
    # -> All values in this order:
    #                 const, y1,   y2,  amp_d, ctr_d, wid_d, amp_g, ctr_g, wid_g 
    initial_params = [0.,    0.45, 0.2, 1e4,   1370., 150.,   1e4,  1595.,  30.] 
    lower_bounds = [-1000.,  -1.5, 0.,  0.,    1340., 50.,   0.,   1550.,  10.] 
    upper_bounds = [1500.,   1.8,  1.4, 5.6e4, 1410., 250.,  5e4,  1610.,  75.]

    # Fit the curve, taking  and tossing away the covariance
    # each parma defined below
    fitted_params, _ = curve_fit(

    # Function to fit curve with
    f = double_lorinzian_extended,

    # Training data
    xdata = x,
    ydata = y,
    
    # Initial parameters and bounds
    p0 = initial_params,
    bounds = (lower_bounds, upper_bounds),
    
    # trf = Trust-Region
    method = 'trf',
    
    # Maximum iterations
    max_nfev = 5e10,
    
    # Tolerance for termination wrt change in loss function 
    ftol = 1e-15,
    # Tolerance for termination wrt change in xdata
    xtol = 1e-15
    )

    
    # Apply fitted parameters to calculate fitted y values
    fitted_y = double_lorinzian_extended(x, *fitted_params)

    #calculate one peak
    #first speratre out fitted amp,center, wid
    par_S1 = fitted_params[3:5]
    par_S2 = fitted_params[6:8]
    
    #fitted values named for legablity                                                                                                                                          
    const_2p = fitted_params[0]
    y1_2p =fitted_params[1]
    y2_2p= fitted_params[2]
    amp_d_2p =fitted_params[3]
    ctr_d_2p =fitted_params[4]
    wid_d_2p =fitted_params[5]
    amp_g_2p =fitted_params[6]
    ctr_g_2p =fitted_params[7]
    wid_g_2p =fitted_params[8]

    #calc y for each individual peak
    peak1 = lorinzian(x,amp_d_2p, ctr_d_2p, wid_d_2p )
    area_d_sims =simps(peak1, dx=5)
    #print(area_d_sims)

    peak2 = lorinzian(x,amp_g_2p, ctr_g_2p, wid_g_2p )
    area_g_sims =simps(peak2, dx=5)

    #calcualte the baseline y
    baseline = const_2p+(y1_2p*x)+(y2_2p*(x*x))
    y_no_base = y - baseline
    fitted_no_base = fitted_y - baseline
    # Plot our fitted curve with the prepared function                     
    extra = np.str(x_cords_list) + np.str(y_cords_list)        
    plot_fitted_data(x, y_no_base, fitted_no_base, peak1, peak2,baseline, extra ) 
    plt.clf()

     

    #better bad data poitns fitting 
    #print(amp_d_2p ,' = amp_d_2p')
    #print(ctr_d_2p, ' = ctr_d_2p')
    #print(wid_d_2p ,' = wid_d_2p')
    #print(amp_g_2p ,' = amp_g_2p')
    #print(ctr_g_2p ,' = ctr_g_2p')
    #print(wid_g_2p ,' = wid_g_2p')

    #print('ratio maker g/d1 =', amp_g_2p/ amp_d_2p)
    ratio = amp_g_2p/ amp_d_2p

    fit_params = (fitted_params, high_y)
   # print("number of fitted params")
    #print(len(fitted_params))

   #print("norm factor", high_y)

    #calcalute temp based on the two peak fit Chan paper                                                              
    fwhm = wid_d_2p*2
    temp2= 899.9-(3*fwhm)+(1.4*0.001*(fwhm*fwhm))

    ctotal = ((amp_d_2p*3.14*wid_d_2p)/2) +((amp_g_2p*3.14*wid_g_2p)/2)

    Fit_params = (const_2p,y1_2p,y2_2p, amp_d_2p,ctr_d_2p,wid_d_2p,amp_g_2p,ctr_g_2p,wid_g_2p, high_y, temp2, x_cords_list, y_cords_list, map_i, area_d_sims, area_g_sims)

    #append fitted_params to fitted_data_all
    if data_fit_all is None:
        data_fit_all = np.array([Fit_params]) 
    else:
        data_fit_all = np.concatenate((data_fit_all, np.array([Fit_params])))

    #print("this is temp2_2lor", temp2)
    #print('amp_d_2p = ', amp_d_2p, 'amp_g_2p =', amp_g_2p)
    if np.float(amp_d_2p) < 0.01 and np.float(amp_g_2p) < 0.01 and ratio >2.5 :
      temp_i = np.append(temp_i, np.NaN)
    else:
      temp_i = np.append(temp_i, temp2)

    if np.float(amp_d_2p) < 0.01 and np.float(amp_g_2p) < 0.01 and ratio >2.5 :
      carbon_amount = np.append(carbon_amount, np.NaN)
    else:
      carbon_amount = np.append(carbon_amount, ctotal)

    cords_x = np.append(cords_x, x_cords_list)
    cords_y = np.append(cords_y, y_cords_list)
    #print(cords_x)

np.savetxt(name+"_2lor_Cody.csv", data_fit_all, delimiter=",")
    
#sys.exit(0)
min_x = cords_x.min()
min_y = cords_y.min()

max_x = cords_x.max()
max_y = cords_y.max()

print(min_x , ' = Min_x')
print(min_y , ' = Min_y')

print(max_x , ' = Max_x')
print(max_y , ' = Max_y')

sorted_x = sorted(np.unique(cords_x))
sorted_y = sorted(np.unique(cords_y))

arranged_data = np.zeros((len(np.unique(cords_x)), len(np.unique(cords_y))))

for idx, temp in enumerate(temp_i):
  x_idx = sorted_x.index(cords_x[idx])
  y_idx = sorted_y.index(cords_y[idx])
  arranged_data[(y_idx, x_idx)] = temp

fig, ax = plt.subplots()

import matplotlib.image as mpimg 
#map_img = mpimg.imread('Tag_overplot.png')

plt.title('Temperature with Cody Equation Pre-Heating')
plt.xlabel('x-axis')
plt.ylabel('y-axis')

sns.heatmap(arranged_data, vmin=0, vmax=400, alpha = 1, cmap=sns.color_palette('Spectral_r', 100), ax=ax, zorder = 1)

#hmax.imshow(map_img,
#          aspect = hmax.get_aspect(),
#          extent = hmax.get_xlim() + hmax.get_ylim(), zorder = 0)


print(sorted_y)
def format_x_ticks(_, idx):
    if idx is None:
        return None
    return sorted_x[idx]

def format_y_ticks(_, idx):
    if idx is None:
        return None
    return sorted_y[idx]


ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_x_ticks))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(format_y_ticks))

# ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
# ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
plt.xticks(rotation=45)
plt.yticks(rotation=315)
plt.savefig(name+"_2lor_cody.png")
plt.show()




plt.clf
