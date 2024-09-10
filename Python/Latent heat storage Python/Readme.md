This folder contains the python code file to simulate a latent heat storage unit composed of parallel rectangular slabs filled with PCM. The slabs are separated by air gaps through which air flows. The model is in 2D, so the length and thickness of the slab is modelled while the conditions along the width are considered constant. There are three main differential equations: for the flowing air, the slab wall and the PCM. The user must give the air, PCM and slab wall properties as well as other parameters. The equations are solved using finite differences. The enthalpy formulation is used for the phase change. The user can input different but constant properties for the liquid and solid phase of the PCM (density, specific heat, and thermal conductivity). The PCM is assumed to melt in a temperature range defined by the user. The user defines the centerpoint of the range and the width of the range. The model can also be used for isothermal phase change by simply using a very narrow temperature range (i.e. 0.1 K) which virtually gives the same results as true isothermal phase change.

The folder also contains sample input and output data files. The user inputs the dimensions and material properties of the LHS parts as well as other parameters. All this is well detailed in the code comments.

The program takes an input data file. In the framework of the present code package, this input file is the output file from the solar air heater program. 
Besides some header lines containing the setup of the solar air heater simulation from which the file comes, each row contains the following data (the names may vary in the file):

- Date (the date)
- Hour_min (the time of day)
- Hour_dec (the time of day in decimal format, e.g. 10:15 is 10.25)
- t_min (the simulation time in minutes)
- Ins (the solar radiation at the time)
- Tamb (the ambient temperature at the time)
- Wind (the wind speed at the time)
- Patm (the atmospheric pressure at the time)
- X (the air humidity ratio at the time)
- Tf_out (the outlet air temperature from the solar air heater program at the time)
- mdot_a (the mass flow rate at the time)


The code contains a function called "extract data" to read the all data rows in the data file.

In reality the only variables needed for the calculations of the LHS unit are mdot_a and Tf_out. The others are read and may be passed to the output files.
