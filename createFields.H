    Info<< "Reading thermodynamicProperties\n" << endl;

    IOdictionary thermodynamicProperties
    (
        IOobject
        (
            "thermodynamicProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar Cv_g(thermodynamicProperties.lookup("Cv_g"));

    dimensionedScalar Cv_l(thermodynamicProperties.lookup("Cv_l"));

    dimensionedScalar Cp_g(thermodynamicProperties.lookup("Cp_g"));

    dimensionedScalar Cp_l(thermodynamicProperties.lookup("Cp_l"));

    dimensionedScalar Lv0(thermodynamicProperties.lookup("Lv0"));

    dimensionedScalar rholSat(thermodynamicProperties.lookup("rholSat"));

    dimensionedScalar rhovSat(thermodynamicProperties.lookup("rhovSat"));

    dimensionedScalar p_inf(thermodynamicProperties.lookup("p_inf"));

    dimensionedScalar pSat(thermodynamicProperties.lookup("pSat"));

    dimensionedScalar pBound(thermodynamicProperties.lookup("pBound"));

    dimensionedScalar C_satL(thermodynamicProperties.lookup("C_satL"));

    dimensionedScalar C_satV(thermodynamicProperties.lookup("C_satV"));

    dimensionedScalar C_G(thermodynamicProperties.lookup("C_G"));

    dimensionedScalar mu_G(thermodynamicProperties.lookup("mu_G"));

    scalar b(readScalar(thermodynamicProperties.lookup("b")));  

    scalar PresTol(readScalar(thermodynamicProperties.lookup("PresTol")));  

    label cav_model(readScalar(thermodynamicProperties.lookup("cav_model")));  

    scalar nbubbles(readScalar(thermodynamicProperties.lookup("nbubbles")));  

    dimensionedScalar Cv(thermodynamicProperties.lookup("Cv"));

    dimensionedScalar Cc(thermodynamicProperties.lookup("Cc"));

    dimensionedScalar UInf(thermodynamicProperties.lookup("UInf"));

    dimensionedScalar tInf(thermodynamicProperties.lookup("tInf"));

    dimensionedScalar a(thermodynamicProperties.lookup("a"));

    dimensionedScalar p_ref(thermodynamicProperties.lookup("p_ref"));

    dimensionedScalar p0(thermodynamicProperties.lookup("p0"));

    dimensionedScalar p1(thermodynamicProperties.lookup("p0"));

    dimensionedScalar rhoMin(thermodynamicProperties.lookup("rhoMin"));

    dimensionedScalar p_cav_coef(thermodynamicProperties.lookup("p_cav_coef"));

    dimensionedScalar Rv(thermodynamicProperties.lookup("Rv"));

    dimensionedScalar Rc(thermodynamicProperties.lookup("Rc"));

    dimensionedScalar T0(thermodynamicProperties.lookup("T0"));

    dimensionedScalar TSat(thermodynamicProperties.lookup("TSat"));

    dimensionedScalar T_ref(thermodynamicProperties.lookup("T_ref"));

    dimensionedScalar Hfg(thermodynamicProperties.lookup("Hfg"));

    dimensionedScalar CvTan(thermodynamicProperties.lookup("CvTan"));

    dimensionedScalar R(thermodynamicProperties.lookup("R"));

    dimensionedScalar RvTan(thermodynamicProperties.lookup("RvTan"));

    dimensionedScalar RcTan(thermodynamicProperties.lookup("RcTan"));

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField mdot
    (
       IOobject
       (
          "mdot",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ
       ),
      mesh,
      dimensionedScalar("mdot", dimMass/(dimLength*dimLength*dimLength*dimTime) , 0.),
      "zeroGradient"
    );

    volScalarField mDotcP
                  (
                       IOobject
                          (
                              "mDotcP",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ
                          ),
                      mesh,
                      dimensionedScalar ("mDotcP",dimMass/(dimLength*dimLength*dimLength*dimTime) , 0.)
                  );

    volScalarField mDotvP
                  (
                      IOobject
                          (
                              "mDotvP",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ
                          ),
                      mesh,
                      dimensionedScalar ("mDotvP",dimMass/(dimLength*dimLength*dimLength*dimTime) , 0.)
                  );


    volScalarField vDotcP
                  (
                      IOobject
                          (
                              "mDotcP",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ
                          ),
                      mesh,
                      dimensionedScalar ("vDotcP",dimMass/(dimTime*dimMass) , 0.)
                  );

    volScalarField vDotvP
                  (
                      IOobject
                          (
                              "mDotvP",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ
                          ),
                      mesh,
                      dimensionedScalar ("vDotvP",dimMass/(dimTime*dimMass) , 0.)
                  );

    volScalarField mDotcAlphal
                  (
                      IOobject
                          (
                              "mDotcAlphal",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ
                          ),
                      mesh,
                      dimensionedScalar ("mDotcAlphal",dimMass/(dimLength*dimLength*dimLength*dimTime) , 0.)
                  );

     volScalarField mDotvAlphal
                  (
                      IOobject
                          (
                              "mDotvAlphal",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ
                          ),
                      mesh,
                      dimensionedScalar ("mDotvAlphal",dimMass/(dimLength*dimLength*dimLength*dimTime) , 0.)
                  );


     volScalarField vDotcAlphal
                  (
                      IOobject
                          (
                              "vDotcAlphal",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ
                          ),
                      mesh,
                      dimensionedScalar ("vDotcAlphal",dimMass/(dimTime*dimMass) , 0.)
                  );

     volScalarField vDotvAlphal
                  (
                      IOobject
                          (
                              "vDotvAlphal",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ
                          ),
                      mesh,
                      dimensionedScalar ("vDotvAlphal",dimMass/(dimTime*dimMass) , 0.)
                  );

     volScalarField mDotT
           (
                IOobject
                   (
                       "mDotT",
                       runTime.timeName(),
                       mesh,
                       IOobject::NO_READ
                   ),
                mesh,
                dimensionedScalar ("mDotT",dimDensity/dimTime,0)
           );


     volScalarField mDotcT
           (
                IOobject
                          (
                              "mDotcT",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ
                          ),
                mesh,
                dimensionedScalar ("mDotcT",dimMass/(dimLength*dimLength*dimLength*dimTime) , 0.)
           );

     volScalarField mDotvT
           (
                IOobject
                          (
                              "mDotvT",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ
                          ),
                mesh,
                dimensionedScalar ("mDotvT",dimMass/(dimLength*dimLength*dimLength*dimTime) , 0.)
           );

     volScalarField vDotcT
           (
                IOobject
                          (
                              "mDotcT",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ
                          ),
                mesh,
                dimensionedScalar ("vDotcT",dimMass/(dimTime*dimMass) , 0.)
           );

     volScalarField vDotvT
           (
                IOobject
                          (
                              "mDotvT",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ
                          ),
                mesh,
                dimensionedScalar ("vDotvT",dimMass/(dimTime*dimMass) , 0.)
           );

    #include "createPhi.H"

    Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
    twoPhaseMixtureThermo twoPhaseProperties(mesh);

    volScalarField& alpha1(twoPhaseProperties.alpha1());
    volScalarField& alpha2(twoPhaseProperties.alpha2());

    Info<< "Reading thermophysical properties\n" << endl;

    volScalarField& p = twoPhaseProperties.p();
    volScalarField& T = twoPhaseProperties.T();
    volScalarField& rho1 = twoPhaseProperties.thermo1().rho();
    const volScalarField& psi1 = twoPhaseProperties.thermo1().psi();
    volScalarField& rho2 = twoPhaseProperties.thermo2().rho();
    const volScalarField& psi2 = twoPhaseProperties.thermo2().psi();

    volScalarField alphal = alpha1;
    volScalarField alphav = alpha2;

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1 + alpha2*rho2
    );

    volScalarField pSat1
    (
        IOobject
        (
            "pSat1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ
        ),

       0.0001*p + pSat  
    );

    volScalarField TSat1
    (
        IOobject
        (
            "TSat1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ
        ),

       0.0001*T + TSat  
    );

    dimensionedScalar pMin(twoPhaseProperties.lookup("pMin"));

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.

    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phi
    );

    volScalarField dgdt
    (
        pos(alpha2)*fvc::div(phi)/max(alpha2, scalar(0.0001))
    );

    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);

    // Construct compressible turbulence model
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New(rho, U, rhoPhi, twoPhaseProperties)
    );

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));
