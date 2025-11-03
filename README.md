# stellar_rotation
Short code that downloads the SAP lightcurve via lightkurve and computes the rotation rates for age estimates.

Optional functionality to estimate the age of the star based on gyrochronology (via gyrointerp).
This is okay since most of the targets that we search for atmospheric escape around are younger targets anyway (with some exceptions)

To do: 
- Input transit parameters (or auto fetch them) 

Usage: 
1) Input the TARGET in *find_rot_period.py* 
2) Run the following blocks to download the lightcurve 
3) Save the light curve as a csv with the TIC identfier and add it to the *list.txt*  
4) Then run *rotation.py*. This will find the rotation period for all the targets in the list text file
5) Get the rotation period (P_rot) and error (P_roterror). Bot the ACF and the LSP periods are reported but i use the ACF one 
6) Run gyrointerp to generate the age estimate and uncertainites. 