
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

    scalar FmdotC(readScalar(thermodynamicProperties.lookup("FmdotC")));  

    scalar FmdotV(readScalar(thermodynamicProperties.lookup("FmdotV")));  

    label cav_model(readScalar(thermodynamicProperties.lookup("cav_model")));  
	
    scalar nbubbles(readScalar(thermodynamicProperties.lookup("nbubbles")));  

