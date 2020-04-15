from os import system as os
import sys


#Specify compilation settings
settings = ' -ffree-form -ffree-line-length-0 -fdefault-real-8 -O3 -Wunused-label '

#Compile all modules
os("gfortran -J mod/ -c"+settings+"parameters.f -o mod/1.o") 
os("gfortran -J mod/ -c"+settings+"constants.f -o mod/2.o") 
os("gfortran -J mod/ -c"+settings+"tensors.f -o mod/4.o") 
os("gfortran -J mod/ -c"+settings+"derivatives.f -o mod/5.o") 
os("gfortran -J mod/ -c"+settings+"analysis.f -o mod/6.o") 
os("gfortran -J mod/ -c"+settings+"initial_conditions.f -o mod/7.o") 
os("gfortran -J mod/ -c"+settings+"IO.f -o mod/8.o") 
os("gfortran -J mod/ -c"+settings+"rk.f -o mod/9.o") 
os("gfortran -J mod/ -c"+settings+"main.f -o mod/10.o") 



#Link all modules
print ('Starting compilation')
os("gfortran mod/*.o -o GO")
print ('Compiled.') 


#Run the code
os("./GO")


