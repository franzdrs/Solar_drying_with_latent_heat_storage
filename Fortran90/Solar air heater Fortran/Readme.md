This folder contains the Fortran 90 code file for the simulation of a back pass, single cover solar air heater (SAH) as well as sample input and output data files. The model is one-dimensional and transient in nature. The main model components are four partial differential equations for cover, absorber plate, air flow and bottom plate. The user must provide the SAH dimensions as well as material properties of the collector components and air. The length of the SAH is divided in a number of sections defined by the user and the equations are solved in the nodes using the finite difference method. The equations are solved iteratively until the convergence criteria (also defined by the user) are met. The user inputs the dimensions and material properties of the SAH parts as well as other parameters. All this is well detailed in the code comments.


The program takes an input data file. In the framework of the present code package, this input file contains weather data and other potentially time-dependent inputs. Each data row contains the following data (the names may vary in the file):

- date (the date)
- hour_min (the time of day)
- time_dec (the time of day in decimal format, e.g. 10:15 is 10.25)
- time_row (the data row number in the file)
- Ins (the solar radiation at the time)
- Tamb (the ambient temperature at the time (dry bulb temperature))
- RH (the relative humidity at the time, in %)
- Wind (the wind speed at the time)
- Patm (the atmospheric pressure at the time)
- Tf_in (the inlet air temperature, which can be the same as ambient or not)
- Xf_in (the inlet air humidity ratio)
- mdot_a (the mass flow rate into the SAH at the time)

The code contains a subroutine called "read_file" to read the data file from a user-defined start row to an end row. The weather data file can therefore represent the weather over an arbitrary time frame and the user can decide the start and end time of the simulation. The user can use a different type of weather data file and change the "read_file" subroutine accordingly.
	
Required information for the calculations of the SAH performance out of this file are time_dec, Ins, Tamb, RH, wind, Patm, Tf_in and mdot_a. RH is used in the read_file function to calculate dew point temperature. If ambient air is what the SAH receives at the inlet, then the inlet air temperature would be that of the ambient, but to make the program more flexible the user can set the inlet conditions different from the ambient. Xf_in is not really needed for the calculations of the SAH but it is passed to the output file in case the downstream
program needs it (for example the fixed bed dryer needs it). The dew point temperature is calculated from the ambient conditions and used to calculate the sky temperature. The time of day in decimal is also used to calculate the sky temperature.

The program writes two output files, one similar to the input file but containing also the outlet air temperature. This file can be used as input for either the latent heat storage or the fixed bed dryer program. The second file contains more detail, with the temperature of the four main SAH parts (cover, absorber, air and bottom plate) at each node in the grid.
