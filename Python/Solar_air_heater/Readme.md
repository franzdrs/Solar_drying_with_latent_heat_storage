This folder contains the python code file for the solar air heater simulation as well as some sample input and output data files.

The program takes an input data file. In the framework of the present code package, this input file is a weather data file with a given format. Besides some header lines containing the setup of the solar air heater simulation, each row contains the following data (the names may vary in the file):

- date (the date)
- hour_min (the time of day)
- time_dec (the time of day in decimal format, e.g. 10:15 is 10.25)
- time_row (the row number)
- Ins (the solar radiation at the time)
- Tamb (the ambient temperature at the time)
- RH (the relative humidity at the time, in %)
- Wind (the wind speed at the time)
- Patm (the atmospheric pressure at the time)
- mdot_a (the mass flow rate into the SAH at the time)

The code contains a subroutine to read the data file from a user-defined start row to an end row. The weather data file can therefore represent the weather over an arbitrary time frame and the 
user can decide the start and end time of the simulation. The user can use a different type of weather data file and change the read_file subroutine accordingly.
Required information for the calculations of the SAH performance are time_dec, Ins, Tamb, wind, Patm, and mdot_a. RH is the usual variable appearing in weather data for the air humidity, but the code
doesn't need it. The read_file subroutine however will use it to calculate the humidity ratio and dew point temperature. The humidity ratio is passed to the output file in case a downstream
program needs it (for example the fixed bed dryer), while the dew point temperature is used to calculate the sky temperature. The time of day in decimal is also used to calculate the sky temperature.



