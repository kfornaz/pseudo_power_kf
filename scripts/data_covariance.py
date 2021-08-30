#by: Pablo Motta

### edit here ###
num_sim = 500
sim_list = []
for i in range( num_sim ):
    sim_list.append( '/home/pablo/codes/data/FlaskSims/test_frequencies/3F/DL11b/aps/cltable' + str(i+1).zfill(4) + '.dat'  )
output_av = '/home/pablo/codes/data/FlaskSims/test_frequencies/3F/DL11b/av.dat'
output_cov = '/home/pablo/codes/data/FlaskSims/test_frequencies/3F/DL11b/cov.dat'
#####

import numpy as np
import os

def read_Cl ( path ):
    data = np.genfromtxt( path )
    row, col = data.shape 
    data_vec = np.reshape ( data[:,1:] , (col-1)*row  )
    return data_vec


average = np.zeros_like( read_Cl ( sim_list[0] )  )
#cov = np.zeros ( ( len(average), len(average) ) )
for i in range( num_sim ):
    average += read_Cl ( sim_list[i] ) / float( num_sim  )

average_matrix = np.zeros_like( np.genfromtxt( sim_list[0]  )   )
average_matrix[:,0] = np.genfromtxt( sim_list[0]  )[:,0]
average_matrix[:,1:] = np.reshape( average , ( average_matrix[:,1:]   ).shape )

np.savetxt( output_av, average_matrix  )
    
f = open(  output_cov , 'w' )  

dim = len(average)
print( dim )
for a in range( dim ):
    cov_a = np.zeros_like( average  ) 
    line = ''
    for i in range( num_sim ):    
        data_vec = np.array( read_Cl ( sim_list[i]  ) )
        cov_a += ( data_vec[a] - average[a] )*( data_vec - average ) / float( num_sim -1 )
    for b in range( dim ):
        f.write (  f"{ cov_a[b] :.16E}" + ' '  ) 
    f.write ( '\n' ) ; f.flush()
    print( str(a) , end = ' ' ,  flush=True )
f.close()

       


   
    
