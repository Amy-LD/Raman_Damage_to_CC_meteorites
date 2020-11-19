#fits the D and G band with two lorentzians
#uses the short heating event modle to fit temperature 
#makes maps 





#decalre libarires, some may not be used 
import numpy as np
import seaborn as sns
import os
import glob
import scipy as sp
import scipy.optimize
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit
from scipy.integrate import simps
import math 

from collections import defaultdict

import csv


#declare variables needed later in code

#this is the name of the outupt/goes in the title of plots or CSV files
name = "Allende_Raman_damage"
name_plots = 'Allende'
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

#function to extract the X and Y cordinate from the file name. Additionaly for Raman damage fitter it also extracts the map number, and hopefuly also energy at some point. 
def get_x_y_from_filename(fname):
  fname_split = [
    item for item in fname.split('_')
    if item != ''
  ]
#  x = fname_split[fname_split.index('X') + 1]
#  y = fname_split[fname_split.index('Y') + 1]
  map = fname_split[fname_split.index('Y--362-7') + 1]
  power = fname_split[fname_split.index('514nm') - 1]
  laser_time_str = fname_split[fname_split.index('M20') + 1]
  #print(laser_time_str)
  return (  map, power, laser_time_str)


#locate correct data and feed it in. 
datadir = '/Users/amylebleu-debartola/Desktop/Raman_Damage/Allende1/center/' #directory where data is 
AllFiles = np.array([])                   #declare empty array to be filled
for file in glob.glob(os.path.join(datadir, 'Allende1*.txt')):       #start loop to go to the directory above, get all files with .txt
    #print(file)
    AllFiles = np.append(AllFiles, file)      #add each file that has .txt to AllFiles array
#test to see if things got fed in properly
#print(AllFiles)

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
    #print(f"Retrieved data for map/trial: {map}")
    
    #print(map)
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
    #print(x)

  
    #translate str laser_time to an int 
    if laser_time_str == '5s4x':
      laser_time = 20
    elif laser_time_str == 's8x9':
      laser_time = 72
    elif laser_time_str == 's10x4':
      laser_time = 40
 
    ##define the type of fit and how they slot together. Adding is fine for Raman damage as chekcing carbon changee
    def lorentzian(x, amp, ctr, wid):
        """Lorentzian for fitting"""
        # Utilize elementwise functions from numpy
        return (amp / (1+((x-ctr)/wid)**2))

    #This is the whole fomula, not just Lorentzians but also baseline 
    def double_lorentzian_extended(x, *params):
        """Combined double Lorentzianw/ polynomial etc"""
      # Destructure and assign parameters passed in
        const, y1, y2, amp_d, ctr_d, wid_d, amp_g, ctr_g, wid_g = params
      # Calculate our double Lorentzian with the baseline added
        return (lorentzian(x, amp_d, ctr_d, wid_d)
                + lorentzian(x, amp_g, ctr_g, wid_g)
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
    f = double_lorentzian_extended,

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
    fitted_y = double_lorentzian_extended(x, *fitted_params)

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
    const = fitted_params[0]
    y1 =fitted_params[1]
    y2= fitted_params[2]
    amp_d =fitted_params[3]
    ctr_d =fitted_params[4]
    wid_d =fitted_params[5]
    amp_g =fitted_params[6]
    ctr_g =fitted_params[7]
    wid_g =fitted_params[8]

    #calculate the fitted y and area for each peak
    interval_for_simps = (x[1]-x[0])

    peak1 = lorentzian(x,amp_d, ctr_d, wid_d )
    area_d_sims =simps(peak1, dx=interval_for_simps)
    true_area_d = area_d_sims*high_y

    peak2 = lorentzian(x,amp_g, ctr_g, wid_g )
    area_g_sims =simps(peak2, dx=interval_for_simps)
    true_area_g = area_g_sims*high_y


    #calcualte the baseline y
    baseline = const+(y1*x)+(y2*(x*x))
    y_no_base = y - baseline
    fitted_no_base = fitted_y - baseline
    # Plot our fitted curve with the prepared function                     
    extra = np.str(map_for_csv)
   # plot_fitted_data(x, y_no_base, fitted_no_base, peak1, peak2,baseline, extra )    


    #print('ratio maker g/d1 =', amp_g/ amp_d)

   
    #calcalute temp based on the two peak fit Schmidt paper 
    fwhm = wid_g*2
    temp2= 1426.6*math.exp(-fwhm/55.55)
    temp_error =  -((temp2*2)/55.55)*wid_g_error

    print(temp2)





    # makes a new variable so that legend later on can have nice labeles 
    #print(power)
    if power == '0-1pct':
      power_float = 0.1
    elif power == '0-5pct':
      power_float = 0.5
    elif power == '1-0pct':
      power_float = 1.0
    elif power == '1pct':
      power_float = 1.0
    elif power == '5-0pct':
      power_float = 5.0
    elif power == '10pct':
      power_float = 10.0
    elif power == '50pct':
      power_float = 50.0
    else:
      raise ValueError("Label for power on file name not valid!")

    #un-normalize const, y1, y2
    const = const *high_y
    y1 = y1*high_y
    y2 = y2 *high_y

    #shove all the the things we care about into a variable to prep for #print to csv
    Fit_params = (const,y1,y2, amp_d,ctr_d,wid_d,amp_g,ctr_g,wid_g, map_for_csv,power, high_y, temp2,temp_error, wid_d_error, ctr_d_error,wid_g_error, ctr_g_error, laser_time, true_area_d, true_area_g)


    #declare an array with the abrivated x,y, power, and highy*y 
    x_y_and_power[Fit_params[9]] = (x, y, Fit_params[10], y*high_y, power_float, laser_time)
    
    #declare an array with the non-abviated x,y, power, and high_y*y, add variables to arrays 
    full_x, full_y = lines[:,0], lines[:,1]
    full_y_high_y= np.amax(full_y)
    full_y_n = full_y/full_y_high_y
    x_y_and_power_full[Fit_params[9]] = (full_x, full_y_n, Fit_params[10], full_y, power_float, laser_time)
    all_fitted_params.append(Fit_params)

    if data_fit_all is None:
        data_fit_all = np.array([Fit_params]) 
    else:
        data_fit_all = np.concatenate((data_fit_all, np.array([Fit_params])))


   

###############Plotting section###################### 

#decalre global var for sorting, sort infomration such that things are correctly sorted temporaly(Ie first trail is first)
map_number_index = 9
power_index = 10
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
area_d_full =[]
area_g_full = []

c_full = []

#binning var, so that each power can have a color in plots 
#binning var, so that each power can have a color in plots 
p_0_1_temp = []
p_0_1_temp_error= []
p_0_1_high_y =[]
p_0_1_map_num =[]
p_0_1_ctr_d=[]
p_0_1_crt_g=[]
p_0_1_ctr_d_error=[]
p_0_1_ctr_g_error=[]
p_0_1_area_d = []
p_0_1_area_g = []
p_0_1_c =[]


p_0_5_temp= []
p_0_5_temp_error= []
p_0_5_high_y =[]
p_0_5_map_num =[]
p_0_5_ctr_d=[]
p_0_5_crt_g=[]
p_0_5_ctr_d_error=[]
p_0_5_ctr_g_error=[]
p_0_5_area_g = []
p_0_5_area_d = []
p_0_5_c =[]

p_1_temp = []
p_1_temp_error= []
p_1_high_y =[]
p_1_map_num =[]
p_1_ctr_d=[]
p_1_crt_g=[]
p_1_ctr_d_error=[]
p_1_ctr_g_error=[]
p_1_area_d = []
p_1_area_g = []
p_1_c =[]

p_5_temp= []
p_5_temp_error= []
p_5_high_y =[]
p_5_map_num =[]
p_5_ctr_d=[]
p_5_crt_g=[]
p_5_ctr_d_error=[]
p_5_ctr_g_error=[]
p_5_area_d = []
p_5_area_g = []
p_5_c =[]

p_10_temp = []
p_10_temp_error= []
p_10_high_y =[]
p_10_map_num =[]
p_10_ctr_d=[]
p_10_crt_g=[]
p_10_ctr_d_error=[]
p_10_ctr_g_error=[]
p_10_area_d = []
p_10_area_g = []
p_10_c =[]

p_50_temp= []
p_50_temp_error= []
p_50_high_y =[]
p_50_map_num =[]
p_50_ctr_d=[]
p_50_ctr_g=[]
p_50_ctr_d_error=[]
p_50_ctr_g_error=[]
p_50_area_d = []
p_50_area_g = []
p_50_c =[]

map_sort =[]
high_y_full =[]
ctr_g_full = []
ctr_d_full =[]

fix_count = 0 

# For every Fit_param tuple, grab the corresponding x, y ,power, high_y, map_num, ctr_d, ctr_derror, ctr_g_error, temp, and temp error
while k < Fit_param_size:
      var_set = sorted_fitted_params[k]
      this_power = var_set[10]
      this_temp_sort = var_set[12]
      this_t_error = var_set[13]
      time = var_set[18]
      this_high_y =var_set[11]
      this_map_num =var_set[9]
      this_ctr_d=var_set[4]
      this_ctr_g=var_set[7]
      this_ctr_d_err=var_set[15]
      this_ctr_g_err=var_set[17]
      this_area_d = var_set[19]
      this_area_g = var_set[20]
      this_c = var_set[0]
      if this_power == '0-1pct':
        power_float = 0.1
        p_0_1_c.append(this_c)
        p_0_1_temp.append(this_temp_sort)
        p_0_1_temp_error.append(this_t_error)
        p_0_1_high_y.append(this_high_y)
        p_0_1_map_num.append(this_map_num)
        p_0_1_ctr_d.append(this_ctr_d)
        p_0_1_crt_g.append(this_ctr_g)
        p_0_1_ctr_d_error.append(this_ctr_d_err)
        p_0_1_ctr_g_error.append(this_ctr_g_err)
        p_0_1_area_d.append(this_area_d)
        p_0_1_area_g.append(this_area_g)
      elif this_power == '0-5pct':
        power_float = 0.5
        p_0_5_c.append(this_c)
        p_0_5_temp.append(this_temp_sort)
        p_0_5_temp_error.append(this_t_error)
        p_0_5_high_y.append(this_high_y)
        p_0_5_map_num.append(this_map_num)
        p_0_5_ctr_d.append(this_ctr_d)
        p_0_5_crt_g.append(this_ctr_g)
        p_0_5_ctr_d_error.append(this_ctr_d_err)
        p_0_5_ctr_g_error.append(this_ctr_g_err)
        p_0_5_area_d.append(this_area_d)
        p_0_5_area_g.append(this_area_g)
      elif this_power == '1pct' or this_power == '1-0pct':
        power_float = 1.0
        p_1_c.append(this_c)
        p_1_temp.append(this_temp_sort)
        p_1_temp_error.append(this_t_error)
        p_1_high_y.append(this_high_y)
        p_1_map_num.append(this_map_num)
        p_1_ctr_d.append(this_ctr_d)
        p_1_crt_g.append(this_ctr_g)
        p_1_ctr_d_error.append(this_ctr_d_err)
        p_1_ctr_g_error.append(this_ctr_g_err)
        p_1_area_d.append(this_area_d)
        p_1_area_g.append(this_area_g)
      elif this_power == '5-0pct':
        power_float = 5.0
        p_5_c.append(this_c)
        p_5_temp.append(this_temp_sort)
        p_5_temp_error.append(this_t_error)
        p_5_high_y.append(this_high_y)
        p_5_map_num.append(this_map_num)
        p_5_ctr_d.append(this_ctr_d)
        p_5_crt_g.append(this_ctr_g)
        p_5_ctr_d_error.append(this_ctr_d_err)
        p_5_ctr_g_error.append(this_ctr_g_err)
        p_5_area_d.append(this_area_d)
        p_5_area_g.append(this_area_g)
      elif this_power == '10pct':        
        power_float = 10.0
        p_10_c.append(this_c)
        p_10_temp.append(this_temp_sort)
        p_10_temp_error.append(this_t_error)
        p_10_high_y.append(this_high_y)
        p_10_map_num.append(this_map_num)
        p_10_ctr_d.append(this_ctr_d)
        p_10_crt_g.append(this_ctr_g)
        p_10_ctr_d_error.append(this_ctr_d_err)
        p_10_ctr_g_error.append(this_ctr_g_err)
        p_10_area_d.append(this_area_d)
        p_10_area_g.append(this_area_g)
      elif this_power == '50pct':
        power_float = 50.0
        p_50_c.append(this_c)
        p_50_temp.append(this_temp_sort)
        p_50_temp_error.append(this_t_error)
        p_50_high_y.append(this_high_y)
        p_50_map_num.append(this_map_num)
        p_50_ctr_d.append(this_ctr_d)
        p_50_ctr_g.append(this_ctr_g)
        p_50_ctr_d_error.append(this_ctr_d_err)
        p_50_ctr_g_error.append(this_ctr_g_err)
        p_50_area_d.append(this_area_d)
        p_50_area_g.append(this_area_g)

      if power_float == 0.1:
        laser_photons =4.43E16
      elif power_float == 0.5:
        laser_photons =4.07E17
      elif power_float == 1:
        laser_photons =7.61E17
      elif power_float == 5:
        laser_photons =5.28E18
      elif power_float == 10:
        laser_photons=1.08E19
      elif power_float == 50:
        laser_photons =5.45E19

      area_d_full.append(this_area_d)
      area_g_full.append(this_area_g)
      laser_photons_per_sec.append(laser_photons)
      ordered_laser_time.append(time)
      temp_sort.append(this_temp_sort)
      temp_error_sort.append(this_t_error)
      map_sort.append(this_map_num)
      high_y_full.append(this_high_y)
      ctr_g_full.append(this_ctr_g)
      ctr_d_full.append(this_ctr_d)
      c_full.append(this_c)
      k = k+1
      
#decalre smaller dose arrays for binning 
#print(fix_count)

Dose_c = 0 
Dose =[]
Dose_0_1 = []
Dose_0_5 = []
Dose_1 = []
Dose_5 =[]
Dose_10 =[]
Dose_50 =[]
i = 0 

#calaclute each individual iradation instance dose, based on read in time and power, then apprend uniqure dose instance to total dose for a running totoal 
while i < Fit_param_size:
  if i == 0:
    Dose_c = laser_photons_per_sec[0] * ordered_laser_time[0]
    tot_power = tot_power + laser_photons_per_sec[0] * ordered_laser_time[0]
    Dose.append(Dose_c)
    if laser_photons_per_sec[i] == 4.43E16:
      Dose_0_1.append(tot_power)
      #print(len(Dose_0_1))
    elif laser_photons_per_sec[i] == 4.07E17:
      Dose_0_5.append(tot_power)
    elif laser_photons_per_sec[i] == 7.61E17:
      Dose_1.append(tot_power)
    elif laser_photons_per_sec[i] == 5.28E18:
      Dose_5.append(tot_power)
    elif laser_photons_per_sec[i] == 1.08E19:
      Dose_10.append(tot_power)
    elif laser_photons_per_sec[i] == 5.45E19:
      Dose_50.append(tot_power)  
    i = i+1
   
  elif tot_power > 0:
    Dose_c = (laser_photons_per_sec[i] * ordered_laser_time[i]) + tot_power
    tot_power = tot_power + laser_photons_per_sec[i] * ordered_laser_time[i]
    Dose.append(Dose_c)
    if laser_photons_per_sec[i] == 4.43E16:
      Dose_0_1.append(tot_power)
      #print(len(Dose_0_1))
    elif laser_photons_per_sec[i] == 4.07E17:
      Dose_0_5.append(tot_power)
    elif laser_photons_per_sec[i] == 7.61E17:
      Dose_1.append(tot_power)
    elif laser_photons_per_sec[i] == 5.28E18:
      Dose_5.append(tot_power)
    elif laser_photons_per_sec[i] == 1.08E19:
      Dose_10.append(tot_power)
    elif laser_photons_per_sec[i] == 5.45E19:
      Dose_50.append(tot_power)  
    i = i+1
    
power_legend = ['all','0.1%', '0.5%', '1%','5%','10%','50%']

#area_d vs dose with colors coresponding to dose power 
plt.figure()
ax = plt.gca()
plt.plot(Dose_0_1, p_0_1_area_d, 'o', c='#AA00FF')
plt.plot(Dose_0_5, p_0_5_area_d, 'o', c='#0050EF')
plt.plot(Dose_1, p_1_area_d, 'o', c='#60A917')
plt.plot(Dose_5, p_5_area_d , 'o', c='#E3C800')
plt.plot(Dose_10, p_10_area_d, 'o', c='#E51400')
plt.plot(Dose_50, p_50_area_d, 'o', c='#8b0000')
plt.plot(Dose, area_d_full, c = 'k')
plt.gca().legend(power_legend)
plt.title('Area_D vs Dose of '+name_plots)
plt.xlabel('Total photons')
plt.ylabel('Area of D Lorenzian')
plt.savefig('Area_d_vs_dose_'+name_plots)
##plt.show()
plt.clf()

    #area_g vs dose with colors coresponding to dose power. 
plt.figure()
ax = plt.gca()
plt.plot(Dose_0_1, p_0_1_area_g, 'o', c='#AA00FF')
plt.plot(Dose_0_5, p_0_5_area_g, 'o', c='#0050EF')
plt.plot(Dose_1, p_1_area_g, 'o', c='#60A917')
plt.plot(Dose_5, p_5_area_g , 'o', c='#E3C800')
plt.plot(Dose_10, p_10_area_g, 'o', c='#E51400')
plt.plot(Dose_50, p_50_area_g, 'o', c='#8b0000')
plt.plot(Dose, area_g_full, c = 'k')
plt.gca().legend(power_legend)
plt.title('Area_G vs Dose of '+name_plots)
plt.xlabel('Total photons')
plt.ylabel('Area of G Lorenzian')
plt.savefig('Area_g_vs_dose_'+name_plots)
##plt.show()
plt.clf()


#temp vs dose with error bars and colors coresponding to dose power 
plt.figure()
ax = plt.gca()
plt.errorbar(Dose_0_1, p_0_1_temp, yerr= p_0_1_temp_error ,c='#AA00FF', fmt='o')
plt.errorbar(Dose_0_5, p_0_5_temp ,yerr = p_0_5_temp_error , c='#0050EF',fmt='o')
plt.errorbar(Dose_1, p_1_temp, yerr = p_1_temp_error ,c='#60A917',fmt='o')
plt.errorbar(Dose_5, p_5_temp, yerr = p_5_temp_error ,c='#E3C800', fmt='o')
plt.errorbar(Dose_10, p_10_temp, yerr = p_10_temp_error ,c='#E51400', fmt='o')
plt.errorbar(Dose_50, p_50_temp, yerr = p_50_temp_error , c='#8b0000',fmt='o' )
plt.plot(Dose,temp_sort, c = 'k')
plt.gca().legend(power_legend)
plt.ylim(-400, 900)
plt.title('Temp vs Dose of '+name_plots)
plt.xlabel('Total photons')
plt.ylabel('Temperature')
plt.savefig('Mini_Temp_vs_dose_'+name_plots)

#temp vs dose with error bars and colors coresponding to dose power, logorythmic
plt.figure()
ax = plt.gca()
plt.errorbar(Dose_0_1, p_0_1_temp, yerr= p_0_1_temp_error ,c='#AA00FF', fmt='o')
plt.errorbar(Dose_0_5, p_0_5_temp ,yerr = p_0_5_temp_error , c='#0050EF',fmt='o')
plt.errorbar(Dose_1, p_1_temp, yerr = p_1_temp_error ,c='#60A917',fmt='o')
plt.errorbar(Dose_5, p_5_temp, yerr = p_5_temp_error ,c='#E3C800', fmt='o')
plt.errorbar(Dose_10, p_10_temp, yerr = p_10_temp_error ,c='#E51400', fmt='o')
plt.errorbar(Dose_50, p_50_temp, yerr = p_50_temp_error , c='#8b0000',fmt='o' )
plt.plot(Dose,temp_sort, c = 'k')
ax.set_xscale('log')
plt.ylim(-400, 900)
plt.gca().legend(power_legend)
plt.title('Temp vs Dose of '+name_plots)
plt.xlabel('Total photons')
plt.ylabel('Temperature')
plt.savefig('Mini_Temp_vs_dose_log_'+name_plots)

#temp vs map number with error bars colors coresponding to dose power 
plt.clf
plt.figure()
ax = plt.gca()
plt.errorbar(p_0_1_map_num, p_0_1_temp, yerr= p_0_1_temp_error ,c='#AA00FF', fmt='o')
plt.errorbar(p_0_5_map_num, p_0_5_temp ,yerr = p_0_5_temp_error , c='#0050EF',fmt='o')
plt.errorbar(p_1_map_num, p_1_temp, yerr = p_1_temp_error ,c='#60A917',fmt='o')
plt.errorbar(p_5_map_num, p_5_temp, yerr = p_5_temp_error ,c='#E3C800', fmt='o')
plt.errorbar(p_10_map_num, p_10_temp, yerr = p_10_temp_error ,c='#E51400', fmt='o')
plt.errorbar(p_50_map_num, p_50_temp, yerr = p_50_temp_error , c='#8b0000',fmt='o' )
plt.plot(map_sort,temp_sort, c = 'k')
plt.gca().legend(power_legend)
plt.title('Temp vs Trial of '+name_plots)
plt.xlabel('Trial number')
plt.ylabel('Temperature')
plt.savefig('Temp_vs_Trial_num_'+name_plots)

#high_y vs map_number, colors coresponding to dose power 
plt.clf
plt.figure()
ax = plt.gca()
plt.errorbar(p_0_1_map_num, p_0_1_high_y,c='#AA00FF', fmt='o')
plt.errorbar(p_0_5_map_num, p_0_5_high_y , c='#0050EF',fmt='o')
plt.errorbar(p_1_map_num, p_1_high_y, c='#60A917',fmt='o')
plt.errorbar(p_5_map_num, p_5_high_y,c='#E3C800', fmt='o')
plt.errorbar(p_10_map_num, p_10_high_y, c='#E51400', fmt='o')
plt.plot(p_50_map_num, p_50_high_y, 'o', c='#8b0000')
plt.plot(map_sort,high_y_full, c = 'k')
plt.gca().legend(power_legend)
plt.title('Fluorescence vs Trial of '+name_plots)
plt.xlabel('Trial number')
plt.ylabel('Fluorescence')
plt.savefig('Flor_vs_Trial_num_'+name_plots)

#ctr_d vs map number with error bars colors coresponding to dose power 
plt.clf
plt.figure()
ax = plt.gca()
plt.errorbar(p_0_1_map_num, p_0_1_ctr_d, yerr= p_0_1_ctr_d_error ,c='#AA00FF', fmt='o')
plt.errorbar(p_0_5_map_num, p_0_5_ctr_d ,yerr = p_0_5_ctr_d_error  , c='#0050EF',fmt='o')
plt.errorbar(p_1_map_num, p_1_ctr_d, yerr = p_1_ctr_d_error  ,c='#60A917',fmt='o')
plt.errorbar(p_5_map_num, p_5_ctr_d, yerr = p_5_ctr_d_error  ,c='#E3C800', fmt='o')
plt.errorbar(p_10_map_num, p_10_ctr_d, yerr = p_10_ctr_d_error  ,c='#E51400', fmt='o')
plt.errorbar(p_50_map_num, p_50_ctr_d, yerr = p_50_ctr_d_error , c='#8b0000',fmt='o' )
plt.plot(map_sort,ctr_d_full, c = 'k')
plt.gca().legend(power_legend)
plt.title('ctr_d vs Trial of '+name_plots)
plt.xlabel('Trial number')
plt.ylabel('Center of D peak')
plt.savefig('ctr_d_vs_Trial_num_'+name_plots)

#ctr_g vs map number with error bars colors coresponding to dose power 
plt.clf
plt.figure()
ax = plt.gca()
plt.errorbar(p_0_1_map_num, p_0_1_crt_g, yerr= p_0_1_ctr_g_error ,c='#AA00FF', fmt='o')
plt.errorbar(p_0_5_map_num, p_0_5_crt_g ,yerr = p_0_5_ctr_g_error  , c='#0050EF',fmt='o')
plt.errorbar(p_1_map_num, p_1_crt_g, yerr = p_1_ctr_g_error  ,c='#60A917',fmt='o')
plt.errorbar(p_5_map_num, p_5_crt_g, yerr = p_5_ctr_g_error  ,c='#E3C800', fmt='o')
plt.errorbar(p_10_map_num, p_10_crt_g, yerr = p_10_ctr_g_error  ,c='#E51400', fmt='o')
plt.errorbar(p_50_map_num, p_50_ctr_g, yerr = p_50_ctr_g_error , c='#8b0000',fmt='o' )
plt.plot(map_sort,ctr_g_full, c = 'k')
plt.gca().legend(power_legend)
plt.title('ctr_g vs Trial of '+name_plots)
plt.xlabel('Trial number')
plt.ylabel('Center of g peak')
plt.savefig('ctr_g_vs_Trial_num_'+name_plots)

plt.clf
plt.figure()
ax = plt.gca()
plt.errorbar(p_0_1_map_num, p_0_1_c,c='#AA00FF', fmt='o')
plt.errorbar(p_0_5_map_num, p_0_5_c, c='#0050EF',fmt='o')
plt.errorbar(p_1_map_num, p_1_c, c='#60A917',fmt='o')
plt.errorbar(p_5_map_num, p_5_c,c='#E3C800', fmt='o')
plt.errorbar(p_10_map_num, p_10_c, c='#E51400', fmt='o')
plt.errorbar(p_50_map_num, p_50_c,  c='#8b0000', fmt='o')
plt.plot(map_sort,c_full, c = 'k')
plt.gca().legend(power_legend)
plt.title('Baseline constant vs trial '+name_plots)
plt.xlabel('Trial number')
plt.ylabel('Baseline Constant')
plt.savefig('Const_vs_Trial_num_'+name_plots)


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
plt.xlabel('Raman shift in cm-1')
plt.ylabel('Normalizied counts')
plt.savefig('First_and_Last_trial_of_'+name_plots+'_normalized_mini')
#plt.show()
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
plt.xlabel('Raman shift in cm-1')
plt.ylabel('Counts')
#ax.set(xlim=(1000,2000), ylim=(first_y[0], first_y[-1]))
plt.savefig('First_and_Last_trial_of_'+name_plots+'_mini')
#plt.show()





#last minus first subtraction plot 
f_l_n_x = first_x[0:200]

last_y = last_y[0:200]
first_y = first_y[0:200]

#talk to chris?? do subtraction set after fixing.
f_l_n_y = last_y-first_y
#print(f_l_n_y)
plt.figure()
ax = plt.gca()
ax.plot(f_l_n_x, f_l_n_y, label = 'Last-First')
ax.legend()
plt.title('Change in spectra from first to last trial of '+name_plots+', normalized')
plt.xlabel('Raman shift in cm-1')
plt.ylabel('Normalizied counts (last trial - first trial)')
plt.savefig('Last_min_First_'+name_plots+'_normalized_mini')
#plt.show()

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

      this_x, this_y_n, this_power, this_y, power_num, laser_time = x_y_and_power[map_number]
      this_color = color_palette[idx]

      ordered_map_numbers.append(map_number)

      plt.plot(this_x, this_y_n, color=this_color)

    plt.gca().legend(ordered_map_numbers)
    plt.title('Power set: ' + this_power +' normalized spectra of '+name_plots)
    plt.xlabel('Raman shift in cm-1')
    plt.ylabel('Normalizied counts')
    plt.savefig('Power_set:' + this_power +'_normalized_spectra_of_'+name_plots+'mini')
    #plt.show()

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
  plt.xlabel('Raman shift in cm-1')
  plt.ylabel('Normalizied counts')
  plt.savefig('Inital_normalized_spectra_of_each_power_of_'+name_plots+'_mini')
  #plt.show()

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
    plt.xlabel('Raman shift in cm-1')
    plt.ylabel('Counts')
    plt.savefig('Power_set:' + this_power +'_spectra_of_'+name_plots+'mini')
    #plt.show()

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
  plt.xlabel('Raman shift in cm-1')
  plt.ylabel('Counts')
  plt.savefig('Inital_spectra_of_each_power_of_'+name_plots+'_mini')
  #plt.show()
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
plt.xlim(200,2000)
ax.legend()
plt.title('First and Last trial of '+name_plots+', normalized, full spectra')
plt.xlabel('Raman shift in cm-1')
plt.ylabel('Normalizied counts')
plt.savefig('First_and_Last_trial_of_'+ name_plots +'_normalized_full')
#plt.show()
plt.clf()
plt.cla()
plt.close('all')

#plot first and last sets on the same graph, not norm
plt.figure()
ax = plt.gca()
ax.plot(first_x_full , first_y_full , label='First trial')
ax.plot(last_x_full , last_y_full , label = 'Last trial')
plt.xlim(200,2000)
ax.legend()
plt.title('First and Last trial of '+name_plots+' , full spectra')
plt.xlabel('Raman shift in cm-1')
plt.ylabel('Counts')
plt.savefig('First_and_Last_trial_of_'+name_plots+'_full')
#plt.show()




#last minus first subtraction plot 
f_l_n_x = first_x_full[0:390]


f_l_n_y = last_y_full[0:390]-first_y_full[0:390]
#print(f_l_n_y)
plt.figure()
ax = plt.gca()
ax.plot(f_l_n_x, f_l_n_y, label = 'Last-First')
ax.legend()
plt.title('Change in spectra from first to last trial of '+name_plots+', normalized, full spectra')
plt.xlim(200,2000)
plt.xlabel('Raman shift in cm-1')
plt.ylabel('Normalizied counts (last trial - first trial)')
plt.savefig('Last_min_First_'+name_plots+'_normalized_full')
#plt.show()


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
      plt.xlim(200,2000)

    plt.gca().legend(ordered_map_numbers)
    plt.title('Power set: ' + this_power +' normalized spectra of '+name_plots+', full spectra')
    plt.xlabel('Raman shift in cm-1')
    plt.ylabel('Normalizied counts')
    plt.savefig('Power_set:' +this_power+'_normalized_spectra_of_'+name_plots+'full')
    #plt.show()

def sort_all_first_param_sets_of_powers_norm(param_sets_by_power, reverse_palette=False):
  # Generate a color palette with a color for each power
  base_palette = sns.color_palette("husl", len(param_sets_by_power))
  color_palette = list(reversed(base_palette)) if reverse_palette else base_palette

  ordered_map_numbers = []
  ordered_power_num = []

  #print(param_sets_by_power.keys())
  for power, param_set in param_sets_by_power.items():
    print(f"{power}: {param_set[0]}")

  for idx, (power, param_set) in enumerate(param_sets_by_power.items()):
    first_param_set_for_power = param_set[0]
    map_number = first_param_set_for_power[map_number_index]

    this_x, this_y_n, this_power, this_y, power_num, laser_time = x_y_and_power_full[map_number]
    this_color = color_palette[idx]
    #print(power_num)

    ordered_map_numbers.append(map_number)
    ordered_power_num.append(power_num)
    plt.plot(this_x, this_y_n, color=this_color)
    plt.xlim(200,2000)
 
  plt.gca().legend(ordered_power_num)
  plt.title('Inital normalized spectra of each power of '+name_plots+', full spectra')
  plt.xlabel('Raman shift in cm-1')
  plt.ylabel('Normalizied counts')
  plt.savefig('Inital_normalized_spectra_of_each_power_of_'+name_plots+'_full')
  #plt.show()

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
      plt.xlim(200,2000)

    plt.gca().legend(ordered_map_numbers)
    
    plt.title('Power set: ' + this_power +' spectra of '+name_plots+', full spectra')
    plt.xlabel('Raman shift in cm-1')
    plt.ylabel('Counts')
    plt.savefig('Power_set:' + this_power +'_spectra_of_'+name_plots+'full')
    #plt.show()

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
    #print(map_number)
    plt.xlim(200,2000)

  plt.gca().legend(ordered_power_num)
  plt.title('Inital spectra of each power of '+name_plots+', full spectra')
  plt.xlabel('Raman shift in cm-1')
  plt.ylabel('Counts')
  plt.savefig('Inital_spectra_of_each_power_of_'+name_plots+'_full')
  #plt.show()
#makes a graph of each power set, 0.1-10, with color coded in acending colors for differnet maps/dose
sort_together_by_power_norm(params_bucketed_by_power, reverse_palette=True)
sort_together_by_power(params_bucketed_by_power, reverse_palette=True)

#Makes one graph, with the first of each set sorted by dose. 
sort_all_first_param_sets_of_powers_norm(params_bucketed_by_power, reverse_palette=False)
sort_all_first_param_sets_of_powers(params_bucketed_by_power, reverse_palette=False)





#save/backup infomation in csv format 
np.savetxt(name+"_2lor_Schmidt.csv", data_fit_all, delimiter=",", fmt='%s')






