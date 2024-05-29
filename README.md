# Solar_drying_with_latent_heat_storage
- This repository contains the code to simulate the performance of a solar air heater (SAH), a latent heat storage unit (LHS) and a fixed bed dryer (FBD).
- The codes implement the numerical solution of mathematical models describing the energy transfer in the SAH and LHS, and the energy and moisture transfer in the FBD. 
- All three programs work by solving the corresponding differential equations using the finite difference method.
- The models offer flexibility in terms of the characteristics of the components. More details are given in the respective folders.
- The user inputs have to be given in the correct SI units (s, kg, kg/m3, J/kg-K, W/m2-K, J/kg, etc). The comments in the code make the expected units obvious.
- The programs work by reading input data from a text file, processing that information together with the user inputs, and produces output file(s) containing the calculation results.
- Although the three programs can be used to simulate the performance of a solar dryer which incorporates a LHS unit, the programs can be used independently to simulate
  the individual components alone or a different arrangement. A conventional solar dryer, for example, could be simulated by applying the output file of the SAH program directly to the
  FBD program.  
- The programs are given "as is" and the user can modify them to suit their needs. This is particularly the case of the FDB. This programs requires some input data for a given agricultural
  product which are not a fixed value but an equation which is specific to a given product. More info on this is given in the folder with the FBD code and in the program comments.
- Reading the comments in the programs from top to bottom should give the user a good idea of how the programs work and where the user is expected to provide input.

- The following published papers give more detail into the mathematical models and their solution procedure:
      10.1016/j.ecmx.2022.100327 (A comparison of steady-state and transient modelling approaches for the performance prediction of solar air heaters)
      10.1016/j.ecmx.2024.100600 (Performance comparison of a fixed-bed solar grain dryer with and without latent heat storage)
