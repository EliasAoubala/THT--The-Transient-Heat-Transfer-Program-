# THT (The Transient Heat Transfer Program)
THT is a One Dimentional Transient Heat Transfer Program which computes and tracks the transient temperature crossection across a rocket engine under the entered combustion conditions near field and wall material properties. This Functioned has been designed to work with composite wall structures which contain more than wall layer.

***Please note this function was written in MATLAB R2020a and may not be compatible with future/earlier versions of MATLAB.***

## Calling the Matlab Function
The MATLAB function can be called in the following format:
```
OneDimentionalHeatTransfer(Tint, Tamb, CombP, matArray, BC)
```
This will return a graph containing the transient temperatures at the set times selected and specfied in the function.

- **Tint**:  This is the Initialisaiton variable which specifies the initial temperature of the Composite wall at time t=0 (K).
- **Tamb**:  This is the Ambient Temperature of the Enviroment form which the heat transfer process takes place  (K).
- **CombP**: This is an Array which includes all required combustion parameters for the bartz analysis.
  - Format: ***[Pr, Tc, gamma, mue, Cp, Pc, c_star, At, dt, M_study, A_study, d_study]***
- **matArray**: This is the material Array containining all the material property parameters and boundaries. Should the material be composite, a new row is created for each layer.
  - Format: ***[name, n, L, K, rho, Cp, T_limit]***
- **BC**: This specifies the Boundary Conditions of the problem to be analysed
  - Format: ***[Fo, t_end, ht, t1, t2, t3, etc...]***

The main Function can be found at ['Matlab Function/OneDimentionalHeatTransfer.m'](OneDimentionalHeatTransfer.m).

An Example Program is also included titled "Example_1.m" found at ['Example Program/Example_1.m'](Untitled2.m).

