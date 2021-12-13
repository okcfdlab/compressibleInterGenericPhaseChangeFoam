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


Some useful references:

https://link.springer.com/article/10.1007/s00231-017-2031-6
https://www.sciencedirect.com/science/article/pii/S0017931016320051#b0430
https://onlinelibrary.wiley.com/doi/epdf/10.1002/htj.21268
https://www.sciencedirect.com/science/article/pii/S0065271708703344 
https://www.sciencedirect.com/science/article/pii/S0149197016300269?via%3Dihub
https://www.sciencedirect.com/science/article/pii/S2666202719300011#bib0018
https://www.sciencedirect.com/science/article/pii/S001793101300207X
https://www.sciencedirect.com/science/article/pii/S0021999100964817?via%3Dihub
https://www.sciencedirect.com/science/article/pii/S0735193312001807
https://github.com/OpenFOAM/OpenFOAM-2.3.x/tree/master/applications/solvers/multiphase/interPhaseChangeFoam
https://www.sciencedirect.com/science/article/pii/S0017931015002008#b0300
https://www.researchgate.net/profile/Jibran_Haider/publication/259898900_Numerical_Modelling_of_Evaporation_and_Condensation_Phenomena/links/5738d95308ae9f741b2bda90/Numerical-Modelling-of-Evaporation-and-Condensation-Phenomena.pdf
https://www.sciencedirect.com/science/article/pii/S0045793015000808
https://onlinelibrary.wiley.com/doi/full/10.1002/fld.3692
https://digitalcommons.mtu.edu/cgi/viewcontent.cgi?referer=https://www.google.com&httpsredir=1&article=1572&context=etdr
https://aip.scitation.org/doi/10.1063/1.5142739
