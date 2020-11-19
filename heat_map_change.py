import numpy as np
import seaborn as sns
import os
import glob
import scipy as sp
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit

import csv

#program reads in data in format change in temp, x, y
#procededs to use this data to produce a heatmap showing the change in temp between two maps

temp_i = np.array([])
cords_x= np.array([])
cords_y= np.array([])

with open('temp_change_btw_62_19_cody.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
            temp_i = np.append(temp_i, row[0])
            cords_x = np.append(cords_x, row[1])
            cords_y = np.append(cords_y, row[2])
            line_count += 1
    

print(cords_x.dtype)
cords_x = cords_x.astype(np.float)
cords_y = cords_y.astype(np.float)
temp_i = temp_i.astype(np.float)

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

plt.title('Temperature change between map19 and map62 Cody Fit')
plt.xlabel('x-axis')
plt.ylabel('y-axis')

sns.heatmap(arranged_data, vmin=-250, vmax=350, alpha = 1, cmap=sns.color_palette('RdBu_r', 100), ax=ax, zorder = 1)

plt.savefig('JBW_Cody_Change_in_temp')
plt.show()


##Schmidt fit ##############################
temp_i = np.array([])
cords_x= np.array([])
cords_y= np.array([])


with open('temp_change_btw_62_19_schmidt.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
            temp_i = np.append(temp_i, row[0])
            cords_x = np.append(cords_x, row[1])
            cords_y = np.append(cords_y, row[2])
            line_count += 1
    

print(cords_x.dtype)
cords_x = cords_x.astype(np.float)
cords_y = cords_y.astype(np.float)
temp_i = temp_i.astype(np.float)

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

plt.title('Temperature change between map19 and map62 Schmidt Fit')
plt.xlabel('x-axis')
plt.ylabel('y-axis')

sns.heatmap(arranged_data, vmin=-250, vmax=350, alpha = 1, cmap=sns.color_palette('RdBu_r', 100), ax=ax, zorder = 1)
plt.savefig('JBW_Schmidt_Change_in_temp')
plt.show()
