This folder contains the python code file to simulate a fixed-bed dryer for agricultural products as well as some sample input and output data files.

As it is, the code is set up for the drying of wheat according to wheat's physical and thermal properties. 
The sorption isotherm equation and the thin-layer drying equation necessary for the drying modelling are specific for a crop and make take different forms.
There exist numerous isotherm and thin-layer equations for many crops, and the equation parameters, which are usually a function of variables such as temperature, relative
humidity and/or air velocity, are found from experiments and publishen in the literature. Care must be taken by the user to ensure that the calculations obtained from those
equations found in the literature for a given crop are transformed to the units used in the programs here. For example, in the calculations relative humidity is not expressed as % but
as decimal, and the product moisture is not expressed in % wet basis but in decimal dry basis (kg/kg). In the literature of both sorption isotherms and thin-layer drying other units might have been used, such as 
%wb or %db for moisture, h or min for time, or % for relative humidity.

In the code for the solar air heater and latent heat storage the code was organized in such a way that the user inputs are all concentrated on the top of the code down to a certain line, after which the rest are the internal calculations. In the case of this code for the fixed bed dryer, if the user wants to simulate the drying of crop different than wheat, he/she must modify the code at other few locations to calculate the values based on the sorption isotherm equation and the thin-layer drying equation. 
