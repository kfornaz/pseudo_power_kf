#By Pablo Motta
#usage:
#python -I uclclfile.dat -O outputfile.dat -M mixmatrix.dat


import numpy as np
import os
import argparse

factor = True #( let it true if the cl file has an extra term of l(l+1)/2pi  )

# Pars section
parser = argparse.ArgumentParser()
parser.add_argument('-I' ,  type=str, required=True)
parser.add_argument('-O' , type=str, required=True)
parser.add_argument('-M' , type=str, required=True)
args = parser.parse_args()
UCLCL_file_path = args.I
mixing_matrix_path = args.M
output_path = args.O


#Read uclcl file
input_data = np.genfromtxt(  UCLCL_file_path   )
ell_uclcl = input_data[:,0]
CLnum = int(len( input_data[0,:] ) - 1 )
binnum = int( - 0.5 + 0.5*np.sqrt(  1 + 8.*float( CLnum   )  ) )
print(  CLnum  )
if( factor == True  ):
    for i in range( CLnum  ):
        input_data[:,int(i+1)] /= ell_uclcl*( ell_uclcl + np.ones_like( ell_uclcl )  )/( 2.*np.pi )


#Read mixing matrix file
mixing_matrix_data = np.genfromtxt( mixing_matrix_path )
ell_mm = mixing_matrix_data[:,0]

#ell handling
ell_max_uclcl = max( ell_uclcl  )
ell_max_mm = max( ell_mm )
ell_max = min( [ell_max_uclcl, ell_max_mm]   )
ellnum = int( ell_max + 1 )
ell = np.arange( ellnum )


#Reshape mixing matrix
mixing_matrix = np.reshape( mixing_matrix_data[:,2] ,  (  int(ell_max_mm + 1) , int(ell_max_mm + 1 ) )   )

#Define output matrix variable
data_output = np.zeros( (  ellnum , int( CLnum + 1 )   )  )
data_output[:,0] = ell

#Convolve
for k in range( CLnum  ):
    CL_temp = input_data[:, k+1]
    CL = np.insert( CL_temp, obj = 0, values = [0.,0.]   )
    SL = np.zeros( ellnum  )
    for i in range( ellnum  ):
        for j in range( ellnum  ):
            SL[i] += mixing_matrix[i,j]*CL[j]
    data_output[:,k+1] = SL[:]


#Header
header = 'Author: Pablo Motta \nAPS convolved with mixing matrix \n 1:l'
index = 0
for i in range(binnum):
    for j in range (binnum):
        if( i<= j  ):
            header += '      ' + str(index+2) + ':temp[' + str(i) + ']-temp[' + str(j) + ']'
            index +=1


#Write output
np.savetxt( output_path , data_output , header = header )            

 

	

