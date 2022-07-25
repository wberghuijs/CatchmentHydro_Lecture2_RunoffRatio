#!/usr/bin/env python
# coding: utf-8

# # Assignment: large-sample hydrology
# ### Wouter Berghuijs (w.r.berghuijs@vu.nl) & Ype van der Velde
# 
# In this Jupyter notebook you will work with a large hydrogical dataset called CAMELS (Newman et al. 2015). The dataset contains information for  671 small- to medium- sized catchments across the contiguous United States (median basin size of 336 km$^2$) that span a  wide range of hydroclimatic conditions.  Daily streamflow data was compiled from the United States Geological Survey National Water Information System. The catchments are minimally impacted by human activities. The catchments have 20 years of daily streamflow data (expressed in mm/d). From these data several **streamflow signatures** are calculate (see Table below). 
# 
# ![Screenshot%202022-07-24%20at%2012.59.00.png](attachment:Screenshot%202022-07-24%20at%2012.59.00.png)
# 
# These streamflow are quantitative metrics or indices that describe statistical or dynamical properties of streamflow. They can be useful in understanding how catchments behave, and how their behaviour various between different catchments. In addition, their applications include extracting biologically relevant attributes of streamflow data, monitoring hydrologic change, analysing runoff generation processes, defining similarity between watersheds, and calibrating and evaluating hydrologic models. Hydrologic signatures allow us to extract meaningful information about watershed processes from streamflow series, and are therefore seeing increasing use in emerging information-rich areas such as global-scale hydrologic modelling, machine learning and large-sample hydrology. 
# 
# The CAMELS dataset also contains information about the climate and landscape of these catchments. The catchment propreties that are available span several classes: topography, climate, land cover, soil, and geology (see Addor et al., 2017)
# 
# 
# ## References
# Addor, N., Newman, A. J., Mizukami, N., & Clark, M. P. (2017). The CAMELS data set: catchment attributes and meteorology for large-sample studies. Hydrology and Earth System Sciences, 21(10), 5293-5313.
# 
# Newman, A. J., Clark, M. P., Sampson, K., Wood, A., Hay, L. E., Bock, A., Viger, R. J., Blodgett, D., Brekke, L., Arnold, J. R., Hopson, T., and Duan, Q.: Development of a large-sample watershed-scale hydrometeorological data set for the contiguous USA: data set characteristics and assessment of regional variability in hydrologic model performance, Hydrol. Earth Syst. Sci., 19, 209–223, https://doi.org/10.5194/hess-19-209-2015, 2015.
# 
# 
# ## Assignment
# 
# Your task in this assignment is to figure out which climate and landscape characteristics control between-catchment differences in streamflow signatures.
# 
# You will pick one of the hydrological signatures that is provided in the CAMELS dataset (see table above).
# 
# Based on this signature you write a short (max 3-page) report that includes
# 1. Cover page, including your names
# 1. List which signature you picked and a short description the signature's relevance.
# 2. A description of how this signature varies across the USA.
# 3. Hypotheses on how 2-3 climate and/or landscape properties may affect your signature substantially
# 4. Graphs that show which factors (most strongly) seem to control your hydrological signatures across the USA
# 5. A discussion of these results.
# 
# We provide code that automatically makes graphics you can use for this assignment. You are free to create extra graphics yourself.

# In[15]:


#get_ipython().system('source activate frontiers_in.yml')
#conda info --envs


# In[8]:


# Import relevant packages and data

#make students do: conda env create -f frontiers_in.yml
# and activate conda env create -f environment.yml 

# Packages
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

import numpy as np
import statistics
from scipy import stats
import seaborn as sns
import geopandas
#import datetime 
import math
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
get_ipython().run_line_magic('matplotlib', 'inline')

# CAMELS data
hydro = pd.read_csv('data/camels_hydro.txt', sep=";", header=0)
geo = pd.read_csv('data/camels_geol.txt', sep=";", header=0)
clim = pd.read_csv('data/camels_clim.txt', sep=";", header=0)
soil = pd.read_csv('data/camels_soil.txt', sep=";", header=0)
veg = pd.read_csv('data/camels_vege.txt', sep=";", header=0)
topo = pd.read_csv('data/camels_topo.txt', sep=";", header=0)
name =  pd.read_csv('data/camels_name.txt', sep=";", header=0)
frames = [hydro,soil, topo, clim, veg]

merged = pd.concat([hydro,soil, topo, clim, veg], axis=1)

hydro.reset_index(drop=True, inplace=True)
soil.reset_index(drop=True, inplace=True)
topo.reset_index(drop=True, inplace=True)
clim.reset_index(drop=True, inplace=True)
veg.reset_index(drop=True, inplace=True)

del merged['gauge_id']

# load shapefile USA
states = geopandas.read_file('shapefiles/usa-states-census-2014.shp')

plt.rcParams["figure.figsize"] = (12, 7)


# In[9]:


## A histogram that allows you to pick which variable you want to display and the number of bins

def run_histogram(variable,bin_count):
    
    ## Python example
    X = merged[variable]
    fig = plt.figure(figsize =(5, 3))
     # Creating axes instance
    ax = fig.add_axes([0, 0, 1, 1])
    # Creating plot
    plt.hist(X, bins=bin_count, alpha=0.5)
    # adding labels
    ax.set_xlabel(variable)
    ax.set_ylabel('Count')
    print(X.describe())
    plt.show()
    
widgets.interact(run_histogram,variable = merged, bin_count=(1,50,1))


# In[10]:


## Scatterplots

def run_scatter(variable_x,variable_y,variable_z,log_x_on,log_y_on,polyfit):

    ## Python example
    X = merged[variable_x]
    Y = merged[variable_y]
    Z = merged[variable_z]

    fig = plt.figure(figsize =(5, 3))
    
    # Creating axes instance
    cm=plt.get_cmap('jet') # picks color based on these options: 
    #https://matplotlib.org/stable/tutorials/colors/colormaps.html 
    
    ax = fig.add_axes([0, 0, 1, 1])
    
    print(min(Z))
    # Creating plot
    sc = plt.scatter(X, Y,c=Z,  vmin=min(Z), vmax=max(Z),s=10, cmap=cm)
    
    # adding labels
    ax.set_xlabel(variable_x)
    ax.set_ylabel(variable_y)
    cbar = plt.colorbar(sc)
    cbar.set_label(variable_z, rotation=270)

    new_X = X[np.logical_not(np.isnan(X*Y))]
    new_Y = Y[np.logical_not(np.isnan(X*Y))]

    #show plot
    linear_model=np.polyfit(new_X,new_Y,polyfit)
    linear_model_fn=np.poly1d(linear_model)
    x_s=np.arange(min(X),max(X))
    plt.plot(x_s,linear_model_fn(x_s),color="green")
    print(linear_model_fn)
    if log_x_on =="Yes": 
        plt.xscale("log")
        
    if log_y_on =="Yes": 
        plt.yscale("log")
        
    plt.show()
    
    print('spearman')
    print(stats.spearmanr(new_X,new_Y))
    print('pearson')
    print(stats.pearsonr(new_X,new_Y))
    
widgets.interact(run_scatter,variable_x = merged,
                             variable_y = merged,
                             variable_z = merged,
                             log_x_on = ["No","Yes"],
                             log_y_on = ["No","Yes"],
                             polyfit = [1, 2, 3, 4, 5, 6, 7, 8])


# In[4]:


# Make a map of a variable and adjust the upper and lower boundary of the colorscale

def map_maker(variable,min_fac,max_fac):
    plt.close()
    
    X = merged[variable]
    states.boundary.plot(color='black')
    cm=plt.get_cmap('jet') # choses color based on these options: https://matplotlib.org/stable/tutorials/colors/colormaps.html 
    sc = plt.scatter(topo.gauge_lon, topo.gauge_lat, c=X,  vmin=min(X)*min_fac, vmax=max(X)*max_fac, s=10, cmap=cm)
    cbar = plt.colorbar(sc, shrink=0.5)
    
    cbar.set_label(variable, rotation=270)
    plt.xlim(-127, -67)
    plt.ylim(24, 50)
    plt.show()

widgets.interact(map_maker,variable = merged,
                min_fac = (0,2,0.1), max_fac = (0,2,0.05))


# In[5]:


# As a reference, we make a correlation matrix of everything against everything. 
# You can manually pick which correlation coefficient is presented in the graph

def correlation_overview(corr_type):
    corrmat = merged.corr(method= corr_type, min_periods = 1)
    f, ax = plt.subplots(figsize =(15, 13))
    sns.heatmap(corrmat, ax = ax, cmap ="BrBG", linewidths = 0.1)

widgets.interact(correlation_overview, corr_type = ['pearson', 'kendall', 'spearman'])


# In[6]:


# In[ ]:




