#by: Pablo Motta
#run:
#python binClmatrix.py -I inpputfile.dat -O outputfile.dat - DL delta_ell
#You dont need to edit anything

import numpy as np
import scipy.interpolate as interpolate
import os
import matplotlib.pyplot as plt
import argparse

# Pars section
parser = argparse.ArgumentParser()
parser.add_argument('-I' ,  type=str, required=True)
parser.add_argument('-O' , type=str, required=True)
parser.add_argument('-DL' , type=int, required=True)
args = parser.parse_args()
SL_file = args.I
output_file = args.O
delta_ell = args.DL
ell_init = 2

#print( SL_file )
#print( output_file )
#print( delta_ell )


#Get the Cl Matrix
data = np.genfromtxt(  SL_file  )
ell_data = data[:,0]

#Compute number of Cl and number of bins
SLnum = len(data[0,:]) -1
num_of_bins = int( - 0.5 + 0.5*np.sqrt(  1 + 8.*float( SLnum   )  ) )
print( SLnum )

#Ell handling 
ell_min = max( [ ell_init , ell_data[0]  ]  )
ell_max_temp =  ell_data[-1]
ellbin_max = int ( (ell_max_temp - ell_min +1.)/delta_ell -1) 
num_of_ellbin = int( ellbin_max + 1  )
print ( num_of_ellbin  )
ell_max = (ellbin_max + 1.)*delta_ell + ell_min -1.
ell = np.arange( ell_min, ell_max_temp + 1  )
print( 'ell min: ' + str(int(ell_min))   )
print( 'ell max: ' + str(int(ell_max))  )


#Bin the data
if( os.path.exists( output_file )  == True ):
        os.remove( output_file )
file_object = open( output_file, 'a'  )        
for a in range( num_of_ellbin  ): 
    ell_start = int(   a*delta_ell  + ell_min )
    file_object.write(  str(ell_start) + ' '  )
    for k in range(  SLnum   ):
        SL = data[ :, k+1]
        norm_a = 0.
        SLbin = 0.
        for m in range( ell_start , int( ell_start + delta_ell)   ):
             l = float(m)
             norm_a += 2.*l +1.
             SLbin += SL[m] * (  2.*l +1.   )
        SLbin /= norm_a
        file_object.write(  str(SLbin) + ' '  )
    file_object.write(  '\n'  )
        




