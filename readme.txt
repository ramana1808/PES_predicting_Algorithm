What this code does : 
----------------------
The code contains the implementation of an algorithm to predict electronic potential energy surface from data points sampled by classical MD on an Azzouz-Bourgis based model potential. Further, quantization of vibrational potential energy surfaces is done and dynamics using Augmented fewest switches surface hopping (AFSSH) is done on these quantizated vibrational surfaces. 

Compotents of the code : 
------------------------
The code is implemented by 3 files. 
1. mod_afssh.f90 - this countains the main code implementing everything 
2. AFSSH.inp - this contains input parameters that are read by the mod_afssh.f90 code
	2.1. Modifying the number of training points sampled, dimensionality of potential, etc. can be done here
3. AFSSH.f90 - this code calls uses mod_afssh.f90 to perform dynamics

Output expected of the code : 
-----------------------------
In the mod_affsh.f90 file, uncommenting these lines in main lead to : 
1. rmse_model - gives the train and test rmse of algorithm 
2. draw_pes - draws the ground state quantized vibrational surface 
3. draw_elec_pot - draws and compares original and predicted electronic pes
4. All commented - runs dynamics using AFSSH
	4.1. In subroutine compute_potential using "call full_pot" will run dynamics on original surface
	4.2. In subroutine compute_potential using "call predicted_pot" will run dynamics on predicted surface   




