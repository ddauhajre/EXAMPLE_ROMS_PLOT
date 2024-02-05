'''
EXAMPLE SCRIPT FOR READING AND ACCESSING ROMS DATA
AND MAKING SOME SIMPLE PLOTS
'''
########################################
#IMPORT things
#Standard python libraries
import os
from netCDF4 import Dataset as netcdf
import numpy as np
import matplotlib.pyplot as plt
import cmocean as cmocean
#ROMS speciific ones (see .py files in directory)
import ROMS_depths as RD
import ROMS_tools as RT
#####################################

#Path of solution (directory)
path_out = '/mnt/kamiya/dauhajre/Simulations/MidCal/L3/NO_WEC/Output_30min/'
#Path of grid file
path_grd = '/mnt/kamiya/dauhajre/Simulations/MidCal/L3/NO_WEC/Input/'

#Grid name
grd_name = 'usw3_grd.nc'

#Open grid netcdf
grd_nc = netcdf(path_grd + grd_name, 'r')


# Grab list of output files (and sort them)
path_files, dirs_out,his_names = os.walk(path_out).next()
his_names.sort()
nfiles = len(his_names)
#So here, his_names[0] is the first output file 

######################################################################
#       EXAMPLE OF MAKING A TIME-SERIES OF TEMPERATURE AT A GRID POINT

#Empty list of temperature values and ocean time
temp_grab=[]
otime = []

#Loop through files 10-->15 (subset of full output)
for n in range(10,15):
    print 'Loading file: ' + his_names[n]
    #Open netcdf of output file his_name[n]
    roms_nc = netcdf(path_out + his_names[n], 'r')
    #Find out how many time-steps are in that file
    nt_nc = len(roms_nc.variables['ocean_time'][:])
    #Loop through all times in the file
    for t in range(nt_nc):
        #Grab temperature at k=-1 (surface) and j=100,i=200
        temp_grab.append(roms_nc.variables['temp'][t,-1,100,200])
        #Save time ('ocean_time')
        otime.append(roms_nc.variables['ocean_time'][t])
    #Close netcdf
    roms_nc.close()


#Convert lists to arrays in case you want to do operations on them
temp_arr = np.asarray(temp_grab)
otime_arr = np.asarray(otime)
#Make a time-series of temperature
plt.figure(figsize=[6,6])
plt.plot(otime_arr/86400.,temp_arr,linewidth=3,color='k')
plt.xlabel('Ocean time [days]', fontsize=14)
plt.ylabel('Temp [deg C]',fontsize=14)
#You can save a figure
plt.savefig('example_temptseries')

#######################################################################



######################################################################
# EXAMPLE OF MAKING A 2-D MAP OF VELOCITY MAGNITUDE AT A TIME

#Grab some file from the output
nlook = 10 #look at file number 10 in his_names
tlook = 3 #look at time-step=3 (starting at 0) in netcdf output file
roms_nc = netcdf(path_out + his_names[nlook], 'r')

#Grab u-velocity  at surface
u_surf = roms_nc.variables['u'][tlook,-1,:,:]
#Grab v-velocity  at surface
v_surf = roms_nc.variables['v'][tlook,-1,:,:]

#See that their shapes are different, indicating the u and v grids
u_surf.shape
v_surf.shape

#Convert u and v to 'rho' points using ROMS_tools functions 
u_rho = RT.u2rho(u_surf)
v_rho = RT.v2rho(v_surf)

#Comput velocity magnitude
mag_uv = np.sqrt(u_rho**2 + v_rho**2)

#Grab land mask from grid file
mask_rho = grd_nc.variables['mask_rho'][:,:]
#Nan-out zero values so mask shows up as white in plots
mask_rho[mask_rho==0.0] = np.nan

#Now make a map of it (and multiply mag_uv *mask_rho)
#Some plotting parameters
cmap_uv = plt.cm.seismic
min_uv = 0
max_uv = 0.5
#Get grid dimensions in km
#pm, pn are 1/dx, 1/dy in ROMS --> they are not always constant, but for this grid it is ok
dx =  1. / grd_nc.variables['pm'][0,0]
dy = 1. / grd_nc.variables['pn'][0,0]
#Get shape of mag_uv
[ny,nx] = mag_uv.shape
im_ext_km = [0,nx*dx*1e-3,0,ny*dy*1e-3]
plt.figure(figsize=[12,6])
plt.imshow(mag_uv*mask_rho,origin='lower',vmin=min_uv,vmax=max_uv,cmap=cmap_uv,extent=im_ext_km)
#Put a colorbar()
cbar = plt.colorbar()
cbar.set_label(r'$\sqrt{u^2+v^2}\;(m/s)$',fontsize=12)
plt.xlabel('km',fontsize=14)
plt.ylabel('km',fontsize=14)
#Save
plt.savefig('example_uvmag_map')

#Close netcdf
roms_nc.close()
#####################################################################


######################################################################
# EXAMPLE OF MAKING A VERTICAL PROFILE OF SOMETHING AT A TIME
#Grab some file from the output
nlook = 10 #look at file number 10 in his_names
tlook = 3 #look at time-step=3 (starting at 0) in netcdf output file
roms_nc = netcdf(path_out + his_names[nlook], 'r')

#Compute depths with ROMS_depths.py tools
z_r,z_w = RD.get_zr_zw_tind(roms_nc,grd_nc,tlook,[0,ny,0,nx])
#z_r are 'rho' point depths, z_w are 'w-point' depths 
#Grab temperature at some spot 
temp_prof = roms_nc.variables['temp'][tlook,:,100,200]
zr_prof   = z_r[:,100,200]

#Make plot
plt.figure(figsize=[5,4])
plt.plot(temp_prof,zr_prof,color='k',linewidth=3)
plt.xlabel('temp [deg C]',fontsize=14)
plt.ylabel('z[m]',fontsize=14)
plt.savefig('example_vert_prof',bbox_inches='tight')



################################################################
