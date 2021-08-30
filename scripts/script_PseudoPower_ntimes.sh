#By Pablo Motta
#Script for run the PseudoPower.py

for i in $(seq -w 1 500)
do
#python PseudoPower.py -O /home/pablo/codes/data/FlaskSims/test_bin965/DLtests/DL1/aps/cltable0$i.dat -X /home/pablo/codes/data/FlaskSims/test_bin965/DLtests/DL1/out/out_ -I /home/pablo/codes/data/FlaskSims/test_bin965/organised/sim0$i/map0.fits 
python binClmatrix.py -I /home/pablo/codes/data/FlaskSims/test_bin965/DLtests/DL1/aps/cltable0$i.dat -O /home/pablo/codes/data/FlaskSims/test_bin965/DLtests/DL12/aps/cltable0"$i".dat -DL 12
done


