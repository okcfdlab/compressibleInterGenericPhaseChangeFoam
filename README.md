compressibleInterGenericPhaseChangeFoam is a solver for OpenFOAM version 2.3 for simulating compressible cryogenic fluids with generic phase change consisting of inertially-driven and thermally-driven phase change.

Solver features:

- Based on the OpenFOAM® solver: compressibleInterFoam 
- Two-fluid mixture flow modeling with separate EOS for each phase 
- Volume of fluid (VOF) phase-fraction based interface capturing approach
- Thermally-driven phase change models of Lee and Scharge-Tanasawa
- Inertially-driven phase change models of Kunz and Merkle 
- Phase change modeling in terms of user defined vaporization and condensation coefficients (Rv and Rc)
- Energy equation solved in terms of static enthalpy 
- Saturation temperature (Tsat) dependence on p (based on Clausius–Clapeyron relation)
- Phase change model coupled with solved temperature and enthalpy fields 
- Ability to read generic thermophysical model libraries (e.g. NG, N2, water) for each phase


Installation/compilation instructions:

- In compilation of the solver, there might be some warnings about unused variables in the code; please disregard these warnings as they are related to some other development tasks. The solver should be compiled with no issue, on typical C++ compiler such as Intel g++, gcc, regardless of the warnings.  
- In case of compilation error regarding the library called twoPhaseMixtureThermo (this should not happen normally, if you use OpenFoam v2.3.x), use the following 

 ```wmake libso twoPhaseMixtureThermo```
 
to compile the library first and then use

```wmake```

to compile the solver. 
