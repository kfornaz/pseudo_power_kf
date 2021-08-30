## By Pablo Motta	
#This script works executing a sequence of functions of the the PseudoPower code 
#It writes the A fullsky or partial sky C_ell's table in the format required for uclcl fitting
#You can define the input and output paths in this script itself or by command line as well


##### DEFAULT SETTING ########################

#GENERAL
#The parent directory of all PseudoPower executables
PseudoPower = '/home/pablo/codes/PseudoPower'
#Indicate Nside
Nside = 256
#Select Lmax for computation 
Lmax = 250
#Select if you want to correct the beam
correct_beam = True
#Put the FWHM 
fwhm = 40.

#PARTIAL SKY
#Select if you want partial sky
partial = True
#Select mask
mask = '/home/pablo/codes/PseudoPower/files/mask_larissa_w.fits'
#Select if you need to remove the average ( i.e.: the maps are absolute temperatures )
remove_average = True
#Select if you need to apply the mask ( i.e.: the maps are still fullsky )
apply_mask = True
#Select if you want Jlm normalisation instead of fsky normalisation
want_jlm = False
#Put ilmjlm of the mask (needed only for jlm normalisation)
ilmjlm = '/home/pablo/codes/PseudoPower/files/ilmjlm_larissa_wL750.dat'


#INPUT FILES ( if you use command line it will be neglected )
#Number of maps
number_maps = 3
maps = [None]*number_maps
#Path of the maps
maps[0] = '/home/pablo/codes/data/FlaskSims/test_bin965/organised/sim0001/map0.fits'
maps[1] = '/home/pablo/codes/data/FlaskSims/test_bin965/organised/sim0001/map0.fits'
maps[2] = '/home/pablo/codes/data/FlaskSims/test_bin965/organised/sim0001/map0.fits'
#path[1] = 

#OUTPUT FILES ( if you use command line it will be neglected )
#Path for intermediate files
path_out = '/home/pablo/codes/data/FlaskSims/test_bin965/out/test1_'
#Path for the cl file
path_cl_table = '/home/pablo/codes/data/FlaskSims/test_bin965/out/test1_cltable.dat'

##############################################


import os
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import subprocess as subprocess
import argparse



# Pars section
parser = argparse.ArgumentParser()
parser.add_argument('-I' , nargs='+',  type=str, required=False)
parser.add_argument('-O' , type=str, required=False)
parser.add_argument('-X' , type=str, required=False)
args = parser.parse_args()
if( args.X != None ):
    path_out = args.X 
if( args.O != None ):
    path_cl_table = args.O
if( args.I != None ):
    path =  args.I
    number_maps = len( path )
else:
    path = maps
          

#Compute alm's
for i in range(number_maps):
    mapi = hp.read_map( path[i] , dtype = float , verbose = False )
    if( remove_average == True ):
        average = sum( mapi)/len(mapi)
        mapi = mapi - average
        path_delta = path_out + 'map' + str(i) + '.fits'
        if( os.path.exists(path_delta)  == True ):
            os.remove( path_delta )
        hp.write_map( path_delta , mapi , dtype = float  )
        path[i] = path_delta
    path_alm = path_out + 'alm' + str(i) + '.fits'
    if( os.path.exists( path_alm )  == True ):
        os.remove( path_alm )
    if ( apply_mask == True ):    
        command = PseudoPower + '/Map2Alm -I ' + path[i] + ' -O ' + path_alm +  ' -L ' + str(Lmax) + ' -m ' + mask   
    else:   
        command = PseudoPower + '/Map2Alm -I ' + path[i] + ' -O ' + path_alm +  ' -L ' + str(Lmax) 
    print( command ); os.system(command)


#Define clmatrix variable
ClMatrix = np.zeros(  (  Lmax ,  1 + int( 0.5*number_maps*( number_maps + 1 ) )  )  )
header = 'Author: Pablo Motta \nPartial sky Pseudo Cl estimator - Computed with PseudoPower code\n 1:l'

#Compute Cl's
index=0
for i in range( number_maps ):
    for j in range( number_maps ):
        if( i <= j ):
                path_alm_i = path_out + 'alm' + str(i) + '.fits'
                path_alm_j = path_out + 'alm' + str(j) + '.fits'
                cl_path = path_out + 'cl' + str(i) + '_' + str(j) + '.fits'
                if(want_jlm == True and partial == True ):
                    command = PseudoPower + '/Alm2Cl -I ' + path_alm_i + ' -O ' + cl_path  + ' -L ' + str(Lmax) + ' -N ' + str(Nside) + ' -m ' + mask + ' -R ' + ilmjlm  + ' -P -o -T ' + path[i] + ' -c ' + path_alm_j + ' ' + path[j]
                elif( partial == True ):
                    command = PseudoPower + '/Alm2Cl -I ' + path_alm_i + ' -O ' + cl_path  + ' -L ' + str(Lmax) + ' -N ' + str(Nside) + ' -m ' + mask + ' -c ' + path_alm_j + ' ' + path[j]
                else:
                    command = PseudoPower + '/Alm2Cl -I ' + path_alm_i + ' -O ' + cl_path  + ' -L ' + str(Lmax) + ' -N ' + str(Nside)  + ' -c ' + path_alm_j + ' ' + path[j]    
                print( command ); os.system(command)


                #Read the Cl
                Cldata = np.genfromtxt( cl_path  )
                ClMatrix[:,index+1] = Cldata[2:,1]
                if( correct_beam == True   ):
                    b_ell = hp.gauss_beam ( fwhm=fwhm * np.pi/10800. , lmax = Lmax-1   )
                    ClMatrix[:,index+1] /= b_ell**2
                header += '      ' + str(index+2) + ':temp[' + str(i) + ']-temp[' + str(j) + ']'
                index +=1

                          

#Write Cl table
ClMatrix[:,0] = Cldata[2:,0]
np.savetxt(  path_cl_table  ,  ClMatrix , header = header  )


#Remove intermediate files
for i in range( number_maps ):
    for j in range( number_maps ):
        if( i <= j ):
                cl_path = path_out + 'cl' + str(i) + '_' + str(j) + '.fits'
                command = 'rm ' + cl_path 
                print( command ); os.system(command) 
                if( i==j ): 
                    path_alm = path_out + 'alm' + str(i) + '.fits'
                    path_delta = path_out + 'map' + str(i) + '.fits'
                    if( os.path.exists( path_alm )  == True ):
                        command = 'rm ' + path_alm  
                    print( command ); os.system(command)   
                    if( os.path.exists( path_delta )  == True ):
                        command = 'rm ' + path_delta  
                    print( command ); os.system(command)
                    
                     

