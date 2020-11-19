##converted matlab to python program 
# Fits 2 lorenzians to data and claculates temp with Cody equation,allows comparision of damages
#takes all center points, outupts code with fitted params, normalization factor, temp, the error in the temp, map number
#also pulls the hotest and coldest points.
#makes several scatter plots with key facctors that varry 
#makes several spectra plots, both normaolized and non nomolarized. 

#SPECIAL NOTE: THIS VERSION DOES NOT INCLUDE X AND Y CORIDANTES AS INPUT DATA DOSNT HAVE THEM! MAKE SURE TO ADD CORDS ON YOUR OWN. 



#decalre libarires, some may not be used 
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

from collections import defaultdict

import csv


#declare variables needed later in code

#this is the name of the outupt/goes in the title of plots or CSV files
name = "Jbilet_Winselwan"
name_plots = 'Jbilet_Winselwan'
#Varibles for finding and fitting and printing the hotest and coldes fits
hottest_x = np.array([])
hottest_y = np.array([])
coldest_x = np.array([])
coldest_y = np.array([])
hotest_temp = 300
coldest_temp =250

#dummy varibles and cordinates variables
cords_x = np.array([])
cords_y = np.array([])
temp_l = np.array([])
temp_i = np.array([])
carbon_amount = np.array([])
map = np.array([])
map_for_csv = np.array([])
x_y_and_power_full = {}
power_float = float

#function to extract the X and Y cordinate from the file name. Additionaly for Raman damage fitter it also extracts the map number, and hopefuly also energy at some point. 
def get_x_y_from_filename(fname):
  fname_split = [
    item for item in fname.split('_')
    if item != ''
  ]
#  x = fname_split[fname_split.index('X') + 1]
#  y = fname_split[fname_split.index('Y') + 1]
  map = fname_split[fname_split.index('Y--90') + 1]
  power = fname_split[fname_split.index('static1125') + 2]
  laser_time_str = fname_split[fname_split.index('static1125') + 1]
  return (  map, power, laser_time_str)


#locate correct data and feed it in. 
datadir = '/Users/amylebleu-debartola/Desktop/Raman_Damage/Jbilet_Winselwan/focus/focus/' #directory where data is 
AllFiles = np.array([])                       #declare empty array to be filled
for file in glob.glob(datadir+'Jbilet_Winselwan1'+'*.txt'):       #start loop to go to the directory above, get all files with .txt
    AllFiles = np.append(AllFiles, file)      #add each file that has .txt to AllFiles array

#test to see if things got fed in properly

#declare golbal variables for loop 
k = len(AllFiles)  #size of AllFiles 
redFlag = 0        #error marker
t2=0               #time counter

#stuff here is to faciliate wrting out into a csv file
#file = open('Murchisan_test_map_python.txt', 'w')
data_fit_all = None


#funciton to plot fitted data, takes x, y, the fitted y, two peaks info, the baseline, and extra is a name string 
def plot_fitted_data(x, y, fitted, peak1, peak2,baseline,extra):
        """Plot data following curve fitting"""
        #plt.style.use('dark_background')
        plt.title("2 Lorentzian fit"+extra)
        plt.xlabel('Ramen shift in cm-1')
        plt.ylabel('Normalizied counts')
        # Captures the most meaningful range of X to visualize plt.xlim((950, 2000))
               # Plot the original values
        
        plt.plot(x, y, '#ffcfea')  #prink
        # Plot the fitted curve
        plt.plot(x, fitted,'#6EB5FF' )  # blue
        plt.plot(x, peak1, '#b19cd9')  #lavandar
        plt.plot(x, peak2, '#afe9ff')  #light blue
        #plt.plot(x, baseline, '#bfb9ff')  # purple?
        plt.show()

#currently unused smoth fuction, doesn't help fitting, was explored 
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

#decalre smaller arrays for the purposre of making plots eaiser post fitting 
all_fitted_params = []
x_y_and_power = {}

#loop to read in each file, process it(find best fit, trasnlate that to temp), and write to a CSV
for j in range (k):  #the range of 0 to the number of files
    t2=t2+1
    
    #read in txt file, convert 2d array to 2 1d arrays, X Y
    CurFileName =AllFiles[j]        # get name of 1 file and path 
    parse_one = CurFileName.split('_')
   # print(parse_one)
    map, power, laser_time_str = get_x_y_from_filename(CurFileName)
    #print(power)
    #x_cords_list, y_cords_list = float(x),float(y)
    map_for_csv = float(map)
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
    y = CurFileY/high_y

  
    #translate str laser_time to an int 
    if laser_time_str == '5s4x':
      laser_time = 20
    elif laser_time_str == '8s9x':
      laser_time = 72
 
    ##define the type of fit and how they slot together. Adding is fine for Raman damage as chekcing carbon changee
    def lorinzian(x, amp, ctr, wid):
        """Gaussian for fitting"""
        # Utilize elementwise functions from numpy
        return (amp / (1+((x-ctr)/wid)**2))

    #This is the whole fomula, not just lorinzians but also baseline 
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
    initial_params = [0.,    0.45, 0.2, 1e4,   1370., 70.,   1e4,  1595.,  30.] 
    lower_bounds = [-1000.,  -1.5, 0.,  0.,    1340., 50.,   0.,   1550.,  10.] 
    upper_bounds = [1500.,   1.8,  1.4, 5.6e4, 1410., 250.,  5e4,  1610.,  75.]

    # Fit the curve, taking parameters and tossing away the covariance
    # each parma defined below
    fitted_params, pcov = curve_fit(

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

    #calculate error array
    error_cords = np.sqrt(np.diag(pcov))
    wid_d_error = error_cords[5]
    wid_g_error = error_cords[8]
    ctr_d_error = error_cords[4]
    ctr_g_error= error_cords[7]    


# print(error)
#    print(fitted_y.shape)
#    print(error.shape)
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
    area_d_sims =simps(peak1, dx=1)
    true_area_d = area_d_sims*high_y

    peak2 = lorinzian(x,amp_g_2p, ctr_g_2p, wid_g_2p )
    area_g_sims =simps(peak2, dx=1)
    true_area_g = area_g_sims*high_y

    #calcualte the baseline y
    baseline = const_2p+(y1_2p*x)+(y2_2p*(x*x))
    y_no_base = y - baseline
    fitted_no_base = fitted_y - baseline
    # Plot our fitted curve with the prepared function                     
    extra = np.str(map_for_csv)
   # plot_fitted_data(x, y_no_base, fitted_no_base, peak1, peak2,baseline, extra )    


    #print('ratio maker g/d1 =', amp_g_2p/ amp_d_2p)
    ratio = amp_g_2p/ amp_d_2p

   
    #calcalute temp based on the two peak fit Cody paper 
    fwhm = wid_d_2p*2
    temp2= 899.9-(3*fwhm)+(1.4*0.001*(fwhm*fwhm))

    #calaclute carbon baseed on theory that the area of the D and G bands together has a carbon inmplication, not proved, likley not 100% correct 
    #more for corlation than actual content(IE lots or not lots of carbon indicator)
    ctotal = ((amp_d_2p*3.14*wid_d_2p)/2) +((amp_g_2p*3.14*wid_g_2p)/2)


    #this is wrong!
    #likley something wrong with inital wid_d_error taking, as it is VERY high. 
    temp_error =  (((wid_d_error)/(wid_d_2p*2)*(3*fwhm))**2 + ((2*(1.4*0.001*(fwhm*fwhm)))*((wid_d_error)/(wid_d_2p*2)))**2)**0.5
    #print(wid_d_error)
    #print(wid_d_2p*2)
    #print(temp_error)

    # makes a new variable so that legend later on can have nice labeles 
    #print(power)
    if power == '0-1pct':
      power_float = 0.1
    elif power == '0-5pct':
      power_float = 0.5
    elif power == '1pct':
      power_float = 1.0
    elif power == '5pct':
      power_float = 5.0
    elif power == '10pct':
      power_float = 10.0

    #print(power_float)
    #shove all the the things we care about into a variable to prep for print to csv
    Fit_params = (const_2p,y1_2p,y2_2p, amp_d_2p,ctr_d_2p,wid_d_2p,amp_g_2p,ctr_g_2p,wid_g_2p, map_for_csv,power, high_y, temp2,temp_error, wid_d_error, ctr_d_error,wid_g_error, ctr_g_error, true_area_d, true_area_g, laser_time)

    #declare an array with the abrivated x,y, power, and highy*y 
    x_y_and_power[Fit_params[9]] = (x, y, Fit_params[10], y*high_y, power_float, laser_time)
    
    #declare an array with the non-abviated x,y, power, and high_y*y, add variables to arrays 
    full_x, full_y = lines[:,0], lines[:,1]
    full_y_high_y= np.amax(full_y)
    full_y_n = full_y/full_y_high_y
    x_y_and_power_full[Fit_params[9]] = (full_x, full_y_n, Fit_params[10], full_y, power_float, laser_time)
    all_fitted_params.append(Fit_params)
    print(all_fitted_params)

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
    plt.clf

###############Plotting section###################### 

#decalre global var for sorting, sort infomration such that things are correctly sorted temporaly(Ie first trial is first)
map_number_index = 9
power_index = 10
print(all_fitted_params)
sorted_fitted_params = sorted(all_fitted_params, key=lambda f: f[map_number_index])

params_bucketed_by_power = defaultdict(list)
for fp in sorted_fitted_params:
  params_bucketed_by_power[fp[power_index]].append(fp)


#start of scatter plot section(temp vs dose, temp vs isntance, high_y(florecince) vs instance, ctr of D and G peaks vs instance )

#decalre variables that are mostly for admin and arrays to fill with sorted stuff 
Fit_param_size = len(sorted_fitted_params)
i = 0
tot_power = 0
Dose=[]
laser_photons_per_sec=[]
laser_photons = 0
ordered_laser_time = []
temp_sort = []
temp_error_sort =[]
k = 0

#binning var, so that each power can have a color in plots 
infocus_temp = []
infocus_temp_error= []
infocus_high_y =[]
infocus_map_num =[]
infocus_ctr_d=[]
infocus_crt_g=[]
infocus_ctr_d_error=[]
infocus_ctr_g_error=[]
infocus_area_d = []
infocus_area_g = []


outfocus_clockwise_1_temp= []
outfocus_clockwise_1_temp_error= []
outfocus_clockwise_1_high_y =[]
outfocus_clockwise_1_map_num =[]
outfocus_clockwise_1_ctr_d=[]
outfocus_clockwise_1_crt_g=[]
outfocus_clockwise_1_ctr_d_error=[]
outfocus_clockwise_1_ctr_g_error=[]
outfocus_clockwise_1_area_d = []
outfocus_clockwise_1_area_g = []


outfocus_clockwise_2_temp = []
outfocus_clockwise_2_temp_error= []
outfocus_clockwise_2_high_y =[]
outfocus_clockwise_2_map_num =[]
outfocus_clockwise_2_ctr_d=[]
outfocus_clockwise_2_crt_g=[]
outfocus_clockwise_2_ctr_d_error=[]
outfocus_clockwise_2_ctr_g_error=[]
outfocus_clockwise_2_area_d = []
outfocus_clockwise_2_area_g = []


outfocus_clockwise_3_temp= []
outfocus_clockwise_3_temp_error= []
outfocus_clockwise_3_high_y =[]
outfocus_clockwise_3_map_num =[]
outfocus_clockwise_3_ctr_d=[]
outfocus_clockwise_3_crt_g=[]
outfocus_clockwise_3_ctr_d_error=[]
outfocus_clockwise_3_ctr_g_error=[]
outfocus_clockwise_3_area_d = []
outfocus_clockwise_3_area_g = []


outfocus_counterclockwise_1_temp = []
outfocus_counterclockwise_1_temp_error= []
outfocus_counterclockwise_1_high_y =[]
outfocus_counterclockwise_1_map_num =[]
outfocus_counterclockwise_1_ctr_d=[]
outfocus_counterclockwise_1_crt_g=[]
outfocus_counterclockwise_1_ctr_d_error=[]
outfocus_counterclockwise_1_ctr_g_error=[]
outfocus_counterclockwise_1_area_d = []
outfocus_counterclockwise_1_area_g = []


outfocus_counterclockwise_2_temp = []
outfocus_counterclockwise_2_temp_error= []
outfocus_counterclockwise_2_high_y =[]
outfocus_counterclockwise_2_map_num =[]
outfocus_counterclockwise_2_ctr_d=[]
outfocus_counterclockwise_2_crt_g=[]
outfocus_counterclockwise_2_ctr_d_error=[]
outfocus_counterclockwise_2_ctr_g_error=[]
outfocus_counterclockwise_2_area_d = []
outfocus_counterclockwise_2_area_g = []


outfocus_counterclockwise_3_temp = []
outfocus_counterclockwise_3_temp_error= []
outfocus_counterclockwise_3_high_y =[]
outfocus_counterclockwise_3_map_num =[]
outfocus_counterclockwise_3_ctr_d=[]
outfocus_counterclockwise_3_crt_g=[]
outfocus_counterclockwise_3_ctr_d_error=[]
outfocus_counterclockwise_3_ctr_g_error=[]
outfocus_counterclockwise_3_area_d = []
outfocus_counterclockwise_3_area_g = []

map_sort =[]
high_y_full =[]
ctr_g_full = []
ctr_d_full =[]
area_d_full =[]
area_g_full =[]
fix_count = 0 

legend_list =[]

# For every Fit_param tuple, grab the corresponding x, y ,power, high_y, map_num, Ctr_d, ctr_derror, ctr_g_error, temp, and temp error
while k < Fit_param_size:
  var_set = sorted_fitted_params[k]
  this_power = var_set[10]
  this_temp_sort = var_set[12]
  this_t_error = var_set[13]
  time = var_set[18]
  this_high_y =var_set[11]
  this_map_num =var_set[9]
  this_area_d = var_set[18]
  this_area_g = var_set[19]
  print(this_map_num)
  #print(type(this_map_num))
  this_ctr_d=var_set[4]
  this_ctr_g=var_set[7]
  this_ctr_d_err=var_set[15]
  this_ctr_g_err=var_set[17]
  if this_map_num == 3.0 or this_map_num == 7.0 or this_map_num == 10.0:
    power_float = 0.1
    if len(infocus_temp) == 0:
      legend_list.append('infocus')
    infocus_temp.append(this_temp_sort)
    infocus_temp_error.append(this_t_error)
    infocus_high_y.append(this_high_y)
    infocus_map_num.append(this_map_num)
    infocus_ctr_d.append(this_ctr_d)
    infocus_crt_g.append(this_ctr_g)
    infocus_ctr_d_error.append(this_ctr_d_err)
    infocus_ctr_g_error.append(this_ctr_g_err)
    infocus_area_d.append(this_area_d)
    infocus_area_g.append(this_area_g)    
    fix_count = fix_count+1

  elif this_map_num == 4.0:
    power_float = 0.1 
    if len(outfocus_clockwise_1_temp) == 0:
      legend_list.append('outfocus_clockwise_1')  
    outfocus_clockwise_1_temp.append(this_temp_sort)
    outfocus_clockwise_1_temp_error.append(this_t_error)
    outfocus_clockwise_1_high_y.append(this_high_y)
    outfocus_clockwise_1_map_num.append(this_map_num)
    outfocus_clockwise_1_ctr_d.append(this_ctr_d)
    outfocus_clockwise_1_crt_g.append(this_ctr_g)
    outfocus_clockwise_1_ctr_d_error.append(this_ctr_d_err)
    outfocus_clockwise_1_ctr_g_error.append(this_ctr_g_err)
    outfocus_clockwise_1_area_d.append(this_area_d)
    outfocus_clockwise_1_area_g.append(this_area_g)

  elif this_map_num == 5:
    power_float = 0.1
    if len(outfocus_clockwise_2_temp) == 0:
      legend_list.append('outfocus_clockwise_2')
    power_float =  0.1
    outfocus_clockwise_2_temp.append(this_temp_sort)
    outfocus_clockwise_2_temp_error.append(this_t_error)
    outfocus_clockwise_2_high_y.append(this_high_y)
    outfocus_clockwise_2_map_num.append(this_map_num)
    outfocus_clockwise_2_ctr_d.append(this_ctr_d)
    outfocus_clockwise_2_crt_g.append(this_ctr_g)
    outfocus_clockwise_2_ctr_d_error.append(this_ctr_d_err)
    outfocus_clockwise_2_ctr_g_error.append(this_ctr_g_err)
    outfocus_clockwise_2_area_d.append(this_area_d)
    outfocus_clockwise_2_area_g.append(this_area_g)   


  elif this_map_num == 6:
    power_float = 0.1 
    if len(outfocus_clockwise_3_temp) == 0:
      legend_list.append('outfocus_clockwise_3')
    power_float =  0.1
    outfocus_clockwise_3_temp.append(this_temp_sort)
    outfocus_clockwise_3_temp_error.append(this_t_error)
    outfocus_clockwise_3_high_y.append(this_high_y)
    outfocus_clockwise_3_map_num.append(this_map_num)
    outfocus_clockwise_3_ctr_d.append(this_ctr_d)
    outfocus_clockwise_3_crt_g.append(this_ctr_g)
    outfocus_clockwise_3_ctr_d_error.append(this_ctr_d_err)
    outfocus_clockwise_3_ctr_g_error.append(this_ctr_g_err)
    outfocus_clockwise_3_area_d.append(this_area_d)
    outfocus_clockwise_3_area_g.append(this_area_g)

  elif this_map_num == 8:
    power_float = 0.1
    if len(outfocus_counterclockwise_1_temp) == 0:
      legend_list.append('outfocus_counterclockwise_1')
    power_float =  0.1
    outfocus_counterclockwise_1_temp.append(this_temp_sort)
    outfocus_counterclockwise_1_temp_error.append(this_t_error)
    outfocus_counterclockwise_1_high_y.append(this_high_y)
    outfocus_counterclockwise_1_map_num.append(this_map_num)
    outfocus_counterclockwise_1_ctr_d.append(this_ctr_d)
    outfocus_counterclockwise_1_crt_g.append(this_ctr_g)
    outfocus_counterclockwise_1_ctr_d_error.append(this_ctr_d_err)
    outfocus_counterclockwise_1_ctr_g_error.append(this_ctr_g_err)
    outfocus_counterclockwise_1_area_d.append(this_area_d)
    outfocus_counterclockwise_1_area_g.append(this_area_g)

  elif this_map_num == 9:
    power_float = 0.1
    if len(outfocus_counterclockwise_2_temp) == 0:
      legend_list.append('outfocus_counterclockwise_2')
    power_float = 0.1
    outfocus_counterclockwise_2_temp.append(this_temp_sort)
    outfocus_counterclockwise_2_temp_error.append(this_t_error)
    outfocus_counterclockwise_2_high_y.append(this_high_y)
    outfocus_counterclockwise_2_map_num.append(this_map_num)
    outfocus_counterclockwise_2_ctr_d.append(this_ctr_d)
    outfocus_counterclockwise_2_crt_g.append(this_ctr_g)
    outfocus_counterclockwise_2_ctr_d_error.append(this_ctr_d_err)
    outfocus_counterclockwise_2_ctr_g_error.append(this_ctr_g_err)
    outfocus_counterclockwise_2_area_d.append(this_area_d)
    outfocus_counterclockwise_2_area_g.append(this_area_g)

  elif this_map_num == 9.5:
    power_float = 0.1
    if len(outfocus_counterclockwise_3_temp) == 0:
      legend_list.append('outfocus_counterclockwise_3')
    power_float = 0.1
    outfocus_counterclockwise_3_temp.append(this_temp_sort)
    outfocus_counterclockwise_3_temp_error.append(this_t_error)
    outfocus_counterclockwise_3_high_y.append(this_high_y)
    outfocus_counterclockwise_3_map_num.append(this_map_num)
    outfocus_counterclockwise_3_ctr_d.append(this_ctr_d)
    outfocus_counterclockwise_3_crt_g.append(this_ctr_g)
    outfocus_counterclockwise_3_ctr_d_error.append(this_ctr_d_err)
    outfocus_counterclockwise_3_ctr_g_error.append(this_ctr_g_err)
    outfocus_counterclockwise_3_area_d.append(this_area_d)
    outfocus_counterclockwise_3_area_g.append(this_area_g)

  if power_float == 0.1:
    laser_photons =4.43E16
    #print(this_power)
    #print(power_float)
  elif power_float == 0.5:
    laser_photons =4.07E17
  elif power_float == 1:
    laser_photons =7.61E17
  elif power_float == 5:
    laser_photons =5.28E18
  elif power_float == 10:
    laser_photons=1.08E19
    print(this_power)
    print(power_float)

  laser_photons_per_sec.append(laser_photons)
  ordered_laser_time.append(time)
  temp_sort.append(this_temp_sort)
  temp_error_sort.append(this_t_error)
  map_sort.append(this_map_num)
  high_y_full.append(this_high_y)
  ctr_g_full.append(this_ctr_g)
  ctr_d_full.append(this_ctr_d)
  area_d_full.append(this_area_d)
  area_g_full.append(this_area_g)
  k = k+1
      
#decalre smaller dose arrays for binning 


#area_vs_focus
plt.clf
plt.figure()
ax = plt.gca()
plt.plot(infocus_map_num, infocus_area_d ,'o',c='#e6194B')
plt.plot(outfocus_clockwise_1_map_num, outfocus_clockwise_1_area_d,'o', c='#f58231')
plt.plot(outfocus_clockwise_2_map_num,outfocus_clockwise_2_area_d,'o',c='#ffe119')
plt.plot(outfocus_clockwise_3_map_num,outfocus_clockwise_3_area_d,'o',c='#3cb44b')
plt.plot(outfocus_counterclockwise_1_map_num, outfocus_counterclockwise_1_area_d,'o',c='#42d4f4')
plt.plot(outfocus_counterclockwise_2_map_num, outfocus_counterclockwise_2_area_d,'o',c='#4363d8')
plt.plot(outfocus_counterclockwise_3_map_num, outfocus_counterclockwise_3_area_d,'o',c='#911eb4')
plt.gca().legend(legend_list)
plt.plot(map_sort,area_d_full, c = 'k')
plt.title('Area_D vs Trial of '+name_plots)
plt.xlabel('Trial number')
plt.ylabel('Area')
plt.savefig('Area_d_vs_Trial_num_'+name_plots)
#area_G_vs_focus
plt.clf
plt.figure()
ax = plt.gca()
plt.plot(infocus_map_num, infocus_area_g ,'o',c='#e6194B')
plt.plot(outfocus_clockwise_1_map_num, outfocus_clockwise_1_area_g,'o', c='#f58231')
plt.plot(outfocus_clockwise_2_map_num,outfocus_clockwise_2_area_g,'o',c='#ffe119')
plt.plot(outfocus_clockwise_3_map_num,outfocus_clockwise_3_area_g,'o',c='#3cb44b')
plt.plot(outfocus_counterclockwise_1_map_num, outfocus_counterclockwise_1_area_g,'o',c='#42d4f4')
plt.plot(outfocus_counterclockwise_2_map_num, outfocus_counterclockwise_2_area_g,'o',c='#4363d8')
plt.plot(outfocus_counterclockwise_3_map_num, outfocus_counterclockwise_3_area_g,'o',c='#911eb4')
plt.gca().legend(legend_list)
plt.plot(map_sort,area_g_full, c = 'k')
plt.title('Area_G vs Trial of '+name_plots)
plt.xlabel('Trial number')
plt.ylabel('Area')
plt.savefig('Area_G_vs_Trial_num_'+name_plots)


#high_y vs map_number, colors coresponding to dose power 
plt.clf
plt.figure()
ax = plt.gca()
plt.errorbar(infocus_map_num, infocus_high_y,c='#e6194B', fmt='o')
plt.errorbar(outfocus_clockwise_1_map_num, outfocus_clockwise_1_high_y , c='#f58231',fmt='o')
plt.errorbar(outfocus_clockwise_2_map_num, outfocus_clockwise_2_high_y, c='#ffe119',fmt='o')
plt.errorbar(outfocus_clockwise_3_map_num, outfocus_clockwise_3_high_y,c='#3cb44b', fmt='o')
plt.errorbar(outfocus_counterclockwise_1_map_num, outfocus_counterclockwise_1_high_y, c='#42d4f4', fmt='o')
plt.errorbar(outfocus_counterclockwise_2_map_num, outfocus_counterclockwise_2_high_y,c='#4363d8', fmt='o')
plt.errorbar(outfocus_counterclockwise_3_map_num, outfocus_counterclockwise_3_high_y,c='#911eb4', fmt='o')
plt.gca().legend(legend_list)
plt.plot(map_sort,high_y_full, c = 'k')
plt.title('Fluorescence vs Trial of '+name_plots)
plt.xlabel('Trial number')
plt.ylabel('Fluorescence')
plt.savefig('Flor_vs_Trial_num_'+name_plots)

#ctr_d vs map number with error bars colors coresponding to dose power 
plt.clf
plt.figure()
ax = plt.gca()
plt.errorbar(infocus_map_num, infocus_ctr_d, yerr= infocus_ctr_d_error ,c='#e6194B', fmt='o')
plt.errorbar(outfocus_clockwise_1_map_num, outfocus_clockwise_1_ctr_d ,yerr = outfocus_clockwise_1_ctr_d_error  , c='#f58231',fmt='o')
plt.errorbar(outfocus_clockwise_2_map_num, outfocus_clockwise_2_ctr_d, yerr = outfocus_clockwise_2_ctr_d_error  ,c='#ffe119',fmt='o')
plt.errorbar(outfocus_clockwise_3_map_num, outfocus_clockwise_3_ctr_d, yerr = outfocus_clockwise_3_ctr_d_error  ,c='#3cb44b', fmt='o')
plt.errorbar(outfocus_counterclockwise_1_map_num, outfocus_counterclockwise_1_ctr_d, yerr = outfocus_counterclockwise_1_ctr_d_error  ,c='#42d4f4', fmt='o')
plt.errorbar(outfocus_counterclockwise_2_map_num, outfocus_counterclockwise_2_ctr_d, yerr = outfocus_counterclockwise_2_ctr_d_error ,c='#4363d8', fmt='o')
plt.errorbar(outfocus_counterclockwise_3_map_num, outfocus_counterclockwise_3_ctr_d, yerr = outfocus_counterclockwise_3_ctr_d_error ,c='#911eb4', fmt='o')
plt.gca().legend(legend_list)
plt.plot(map_sort,ctr_d_full, c = 'k')
plt.title('Ctr_d vs Trial of '+name_plots)
plt.xlabel('Trial number')
plt.ylabel('Center of D peak')
plt.savefig('Ctr_d_vs_Trial_num_'+name_plots)

#ctr_g vs map number with error bars colors coresponding to dose power 
plt.clf
plt.figure()
ax = plt.gca()
plt.errorbar(infocus_map_num, infocus_crt_g, yerr= infocus_ctr_g_error ,c='#e6194B', fmt='o')
plt.errorbar(outfocus_clockwise_1_map_num, outfocus_clockwise_1_crt_g ,yerr = outfocus_clockwise_1_ctr_g_error  , c='#f58231',fmt='o')
plt.errorbar(outfocus_clockwise_2_map_num, outfocus_clockwise_2_crt_g, yerr = outfocus_clockwise_2_ctr_g_error  ,c='#ffe119',fmt='o')
plt.errorbar(outfocus_clockwise_3_map_num, outfocus_clockwise_3_crt_g, yerr = outfocus_clockwise_3_ctr_g_error  ,c='#3cb44b', fmt='o')
plt.errorbar(outfocus_counterclockwise_1_map_num, outfocus_counterclockwise_1_crt_g, yerr = outfocus_counterclockwise_1_ctr_g_error  ,c='#42d4f4', fmt='o')
plt.errorbar(outfocus_counterclockwise_2_map_num, outfocus_counterclockwise_2_crt_g, yerr = outfocus_counterclockwise_2_ctr_g_error ,c='#4363d8', fmt='o')
plt.errorbar(outfocus_counterclockwise_3_map_num, outfocus_counterclockwise_3_crt_g, yerr = outfocus_counterclockwise_3_ctr_g_error ,c='#911eb4', fmt='o')
plt.gca().legend(legend_list)
plt.plot(map_sort,ctr_g_full, c = 'k')
plt.title('Ctr_g vs Trial of '+name_plots)
plt.xlabel('Trial number')
plt.ylabel('Center of g peak')
plt.savefig('Ctr_g_vs_Trial_num_'+name_plots)


#section for abrivtaed ragne spectra plots; plots first and last,(norm/notnorm), last-first(norm/notnorm), subtraction series, power set(norm/not norm, first of each set(norm/not norm)

#plot first and last plot, normolized
first_set_params = sorted_fitted_params[0]
last_set_params= sorted_fitted_params[-1]

#print(x_y_and_power[first_set_params[9]])
first_set = x_y_and_power[first_set_params[9]]
last_set = x_y_and_power[last_set_params[9]]

first_x = first_set[0]
first_y_n = first_set[1]
first_y = first_set[3]
last_x = last_set[0]
last_y_n = last_set[1]
last_y = last_set[3]

#quick translation
ax = plt.gca()

# plot first and last sets on the same graph, normalized 
plt.figure()
ax = plt.gca()
ax.plot(first_x, first_y_n, label='First trial')
ax.plot(last_x, last_y_n, label = 'Last trial')
ax.legend()
plt.title('First and Last trial of '+name_plots+', normalized')
plt.xlabel('Ramen shift in cm-1')
plt.ylabel('Normalizied counts')
plt.savefig('First_and_Last_trial_of_'+name_plots+'_normalized_mini')
plt.show()
plt.clf()
plt.cla()
plt.close('all')
plt.savefig('First_and_Last_trial_of_'+name_plots+'_normalized_mini')


#plot first and last sets on the same graph, not norm
plt.clf()
plt.figure()
ax = plt.gca()
ax.plot(first_x, first_y, label='First trial')
ax.plot(last_x, last_y, label = 'Last trial')
ax.legend()
plt.title('First and Last trial of '+name_plots)
plt.xlabel('Ramen shift in cm-1')
plt.ylabel('Counts')
#ax.set(xlim=(1000,2000), ylim=(first_y[0], first_y[-1]))
plt.savefig('First_and_Last_trial_of_'+name_plots+'_mini')
plt.show()





#last minus first subtraction plot 
f_l_n_x = first_x

#talk to chris?? do subtraction set after fixing.
f_l_n_y = last_y-first_y
#print(f_l_n_y)
plt.figure()
ax = plt.gca()
ax.plot(first_x, f_l_n_y, label = 'Last-First')
ax.legend()
plt.title('Change in spectra from first to last trial of '+name_plots+', normalized')
plt.xlabel('Ramen shift in cm-1')
plt.ylabel('Normalizied counts (last trial - first trial)')
plt.savefig('Last_min_First_'+name_plots+'_normalized_mini')
plt.show()

#do an automatic sort by power to assign each plot a new set of colros, plots all spectra of a cetin power on same plot 
def sort_together_by_power_norm(param_sets_by_power, reverse_palette=False):
  # Each power has a corresponding set of Fit_param sets
  for power, param_set in param_sets_by_power.items():
    # Generate a sequential color palette with a hue for every set of Fit_params in this power's set
    base_palette = sns.color_palette("GnBu_d", len(param_set))
    color_palette = list(reversed(base_palette)) if reverse_palette else base_palette

    ordered_map_numbers = []

    # For every Fit_param tuple, grab the corresponding x, y , and power by map number
    # then, grab the color for that particular set of parameters
    # then, plot the x and y with that color
    for idx, params in enumerate(param_set):
      map_number = params[map_number_index]

      this_x, this_y_n, this_map_num, this_y, power_num, laser_time = x_y_and_power[map_number]
      this_color = color_palette[idx]

      ordered_map_numbers.append(map_number)

      plt.plot(this_x, this_y_n, color=this_color)

    plt.gca().legend(ordered_map_numbers)
    plt.title('Power set: ' + this_map_num +' normalized spectra of '+name_plots)
    plt.xlabel('Ramen shift in cm-1')
    plt.ylabel('Normalizied counts')
    plt.savefig('Power_set:' + this_map_num +'_normalized_spectra_of_'+name_plots+'mini')
    plt.show()

# plot just the frist instance of each power on a plot to show large-scale cahnges due to that power 
def sort_all_first_param_sets_of_powers_norm(param_sets_by_power, reverse_palette=False):
  # Generate a color palette with a color for each power
  base_palette = sns.color_palette("husl", len(param_sets_by_power))
  color_palette = list(reversed(base_palette)) if reverse_palette else base_palette

  ordered_map_numbers = []
  ordered_power_num = []

  for idx, (power, param_set) in enumerate(param_sets_by_power.items()):
    first_param_set_for_power = param_set[0]
    map_number = first_param_set_for_power[map_number_index]

    this_x, this_y_n, this_power, this_y, power_num, laser_time = x_y_and_power[map_number]
    this_color = color_palette[idx]

    ordered_map_numbers.append(map_number)
    ordered_power_num.append(power_num)
    plt.plot(this_x, this_y_n, color=this_color)
 
  plt.gca().legend(ordered_power_num)
  plt.title('Inital normalized spectra of each power of '+name_plots)
  plt.xlabel('Ramen shift in cm-1')
  plt.ylabel('Normalizied counts')
  plt.savefig('Inital_normalized_spectra_of_each_power_of_'+name_plots+'_mini')
  plt.show()

#these are the same as above but not normalized. 
def sort_together_by_power(param_sets_by_power, reverse_palette=False):
  # Each power has a corresponding set of Fit_param sets
  for power, param_set in param_sets_by_power.items():
    # Generate a sequential color palette with a hue for every set of Fit_params in this power's set
    base_palette = sns.color_palette("GnBu_d", len(param_set))
    color_palette = list(reversed(base_palette)) if reverse_palette else base_palette

    ordered_map_numbers = []

    # For every Fit_param tuple, grab the corresponding x, y , and power by map number
    # then, grab the color for that particular set of parameters
    # then, plot the x and y with that color
    for idx, params in enumerate(param_set):
      map_number = params[map_number_index]

      this_x, this_y_n, this_power, this_y, power_num, laser_time = x_y_and_power[map_number]
      this_color = color_palette[idx]

      ordered_map_numbers.append(map_number)
      
      plt.plot(this_x, this_y, color=this_color)

    plt.gca().legend(ordered_map_numbers)
    
    plt.title('Power set: ' + this_power +' spectra of '+name_plots)
    plt.xlabel('Ramen shift in cm-1')
    plt.ylabel('Counts')
    plt.savefig('Power_set:' + this_power +'_spectra_of_'+name_plots+'mini')
    plt.show()

def sort_all_first_param_sets_of_powers(param_sets_by_power, reverse_palette=False):
  # Generate a color palette with a color for each power
  base_palette = sns.color_palette("husl", len(param_sets_by_power))
  color_palette = list(reversed(base_palette)) if reverse_palette else base_palette

  ordered_map_numbers = []
  ordered_power_num = []
  for idx, (power, param_set) in enumerate(param_sets_by_power.items()):
    first_param_set_for_power = param_set[0]
    map_number = first_param_set_for_power[map_number_index]

    this_x, this_y_n, this_power, this_y, power_num, laser_time = x_y_and_power[map_number]
    this_color = color_palette[idx]

    ordered_map_numbers.append(map_number)
    ordered_power_num.append(power_num)
    plt.plot(this_x, this_y, color=this_color)

  plt.gca().legend(ordered_power_num)
  plt.title('Inital spectra of each power of '+name_plots)
  plt.xlabel('Ramen shift in cm-1')
  plt.ylabel('Counts')
  plt.savefig('Inital_spectra_of_each_power_of_'+name_plots+'_mini')
  plt.show()
#makes a graph of each power set, 0.1-10, with color coded in acending colors for differnet maps/dose
sort_together_by_power_norm(params_bucketed_by_power, reverse_palette=True)
sort_together_by_power(params_bucketed_by_power, reverse_palette=True)

#Makes one graph, with the first of each set sorted by dose. 
sort_all_first_param_sets_of_powers_norm(params_bucketed_by_power, reverse_palette=False)
sort_all_first_param_sets_of_powers(params_bucketed_by_power, reverse_palette=False)




########## section for non abrivated plots. makes same plots. ##################
#x_y_and_power_full_full- series to draw x/y/ynorm form 
first_set = x_y_and_power_full[first_set_params[9]]
last_set = x_y_and_power_full[last_set_params[9]]

first_x_full = first_set[0]
first_y_n_full = first_set[1]
first_y_full  = first_set[3]
last_x_full  = last_set[0]
last_y_n_full  = last_set[1]
last_y_full  = last_set[3]

#quick translation
ax = plt.gca()

# plot first and last sets on the same graph, normalized 
plt.figure()
ax = plt.gca()
ax.plot(first_x_full , first_y_n_full , label='First trial')
ax.plot(last_x_full , last_y_n_full , label = 'Last trial')
ax.legend()
plt.title('First and Last trial of '+name_plots+', normalized, full spectra')
plt.xlabel('Ramen shift in cm-1')
plt.ylabel('Normalizied counts')
plt.savefig('First_and_Last_trial_of_'+ name_plots +'_normalized_full')
plt.show()
plt.clf()
plt.cla()
plt.close('all')

#plot first and last sets on the same graph, not norm
plt.figure()
ax = plt.gca()
ax.plot(first_x_full , first_y_full , label='First trial')
ax.plot(last_x_full , last_y_full , label = 'Last trial')
ax.legend()
plt.title('First and Last trial of '+name_plots+' , full spectra')
plt.xlabel('Ramen shift in cm-1')
plt.ylabel('Counts')
plt.savefig('First_and_Last_trial_of_'+name_plots+'_full')
plt.show()




#last minus first subtraction plot 
f_l_n_x = first_x_full


f_l_n_y = last_y_full-first_y_full
#print(f_l_n_y)
plt.figure()
ax = plt.gca()
ax.plot(first_x_full, f_l_n_y, label = 'Last-First')
ax.legend()
plt.title('Change in spectra from first to last trial of '+name_plots+', normalized, full spectra')
plt.xlabel('Ramen shift in cm-1')
plt.ylabel('Normalizied counts (last trial - first trial)')
plt.savefig('Last_min_First_'+name_plots+'_normalized_full')
plt.show()


def sort_together_by_power_norm(param_sets_by_power, reverse_palette=False):
  # Each power has a corresponding set of Fit_param sets
  for power, param_set in param_sets_by_power.items():
    # Generate a sequential color palette with a hue for every set of Fit_params in this power's set
    base_palette = sns.color_palette("GnBu_d", len(param_set))
    color_palette = list(reversed(base_palette)) if reverse_palette else base_palette

    ordered_map_numbers = []

    # For every Fit_param tuple, grab the corresponding x, y , and power by map number
    # then, grab the color for that particular set of parameters
    # then, plot the x and y with that color
    for idx, params in enumerate(param_set):
      map_number = params[map_number_index]

      this_x, this_y_n, this_power, this_y, power_num, laser_time = x_y_and_power_full[map_number]
      this_color = color_palette[idx]

      ordered_map_numbers.append(map_number)

      plt.plot(this_x, this_y_n, color=this_color)

    plt.gca().legend(ordered_map_numbers)
    plt.title('Power set: ' + this_power +' normalized spectra of '+name_plots+', full spectra')
    plt.xlabel('Ramen shift in cm-1')
    plt.ylabel('Normalizied counts')
    plt.savefig('Power_set:' +this_power+'_normalized_spectra_of_'+name_plots+'full')
    plt.show()

def sort_all_first_param_sets_of_powers_norm(param_sets_by_power, reverse_palette=False):
  # Generate a color palette with a color for each power
  base_palette = sns.color_palette("husl", len(param_sets_by_power))
  color_palette = list(reversed(base_palette)) if reverse_palette else base_palette

  ordered_map_numbers = []
  ordered_power_num = []

  for idx, (power, param_set) in enumerate(param_sets_by_power.items()):
    first_param_set_for_power = param_set[0]
    map_number = first_param_set_for_power[map_number_index]

    this_x, this_y_n, this_power, this_y, power_num, laser_time = x_y_and_power_full[map_number]
    this_color = color_palette[idx]

    ordered_map_numbers.append(map_number)
    ordered_power_num.append(power_num)
    plt.plot(this_x, this_y_n, color=this_color)
 
  plt.gca().legend(ordered_power_num)
  plt.title('Inital normalized spectra of each power of '+name_plots+', full spectra')
  plt.xlabel('Ramen shift in cm-1')
  plt.ylabel('Normalizied counts')
  plt.savefig('Inital_normalized_spectra_of_each_power_of_'+name_plots+'_full')
  plt.show()

def sort_together_by_power(param_sets_by_power, reverse_palette=False):
  # Each power has a corresponding set of Fit_param sets
  for power, param_set in param_sets_by_power.items():
    # Generate a sequential color palette with a hue for every set of Fit_params in this power's set
    base_palette = sns.color_palette("GnBu_d", len(param_set))
    color_palette = list(reversed(base_palette)) if reverse_palette else base_palette

    ordered_map_numbers = []

    # For every Fit_param tuple, grab the corresponding x, y , and power by map number
    # then, grab the color for that particular set of parameters
    # then, plot the x and y with that color
    for idx, params in enumerate(param_set):
      map_number = params[map_number_index]

      this_x, this_y_n, this_power, this_y, power_num, laser_time = x_y_and_power_full[map_number]
      this_color = color_palette[idx]

      ordered_map_numbers.append(map_number)
      
      plt.plot(this_x, this_y, color=this_color)

    plt.gca().legend(ordered_map_numbers)
    
    plt.title('Power set: ' + this_power +' spectra of '+name_plots+', full spectra')
    plt.xlabel('Ramen shift in cm-1')
    plt.ylabel('Counts')
    plt.savefig('Power_set:' + this_power +'_spectra_of_'+name_plots+'full')
    plt.show()

def sort_all_first_param_sets_of_powers(param_sets_by_power, reverse_palette=False):
  # Generate a color palette with a color for each power
  base_palette = sns.color_palette("husl", len(param_sets_by_power))
  color_palette = list(reversed(base_palette)) if reverse_palette else base_palette

  ordered_map_numbers = []
  ordered_power_num = []
  for idx, (power, param_set) in enumerate(param_sets_by_power.items()):
    first_param_set_for_power = param_set[0]
    map_number = first_param_set_for_power[map_number_index]

    this_x, this_y_n, this_power, this_y, power_num, laser_time = x_y_and_power_full[map_number]
    this_color = color_palette[idx]

    ordered_map_numbers.append(map_number)
    ordered_power_num.append(power_num)
    plt.plot(this_x, this_y, color=this_color)

  plt.gca().legend(ordered_power_num)
  plt.title('Inital spectra of each power of '+name_plots+', full spectra')
  plt.xlabel('Ramen shift in cm-1')
  plt.ylabel('Counts')
  plt.savefig('Inital_spectra_of_each_power_of_'+name_plots+'_full')
  plt.show()
#makes a graph of each power set, 0.1-10, with color coded in acending colors for differnet maps/dose
sort_together_by_power_norm(params_bucketed_by_power, reverse_palette=True)
sort_together_by_power(params_bucketed_by_power, reverse_palette=True)

#Makes one graph, with the first of each set sorted by dose. 
sort_all_first_param_sets_of_powers_norm(params_bucketed_by_power, reverse_palette=False)
sort_all_first_param_sets_of_powers(params_bucketed_by_power, reverse_palette=False)




#save/backup infomation in csv format 
np.savetxt(name+"_2lor.csv", data_fit_all, delimiter=",", fmt='%s')






