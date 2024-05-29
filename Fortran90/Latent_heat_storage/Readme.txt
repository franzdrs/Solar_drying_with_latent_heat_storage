This folder contains the fortran90 code file to simulate a latent heat storage unit composed of rectangular slabs filled with PCM. 
It also contains some sample input and output data files.

The program takes an input data file. Im the framework of the present code package, this input file is the output file from the solar air heater program. 
Besides some header lines containing the setup of the solar air heater simulation, each row contains the following data:
Date (the date)
Hour_min (the time of day)
Hour_dec (the time of day in decimal format, e.g. 10:15 is 10.25)
t_min (the simulation time in minutes)
Ins (the solar radiation at the time)
Tamb (the ambient temperature at the time)
Wind (the wind speed at the time)
Patm (the atmospheric pressure at the time)
X (the air humidity ratio at the time)
Tout_avg (the outlet air temperature from the solar air heater program)
mdot_a (the mass flow rate at the time)

In reality the only variables for the calculations of the LHS unit are mdot_a and Tout_avg. THe others are read to be passed to the output files. 
