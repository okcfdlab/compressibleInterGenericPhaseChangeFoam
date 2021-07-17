{

    volScalarField limitedAlpha1 = (min(max(alpha1, scalar(0)), scalar(1)));
    volScalarField rhoCp = rho1*twoPhaseProperties.thermo1().Cp()*alpha1 + rho2*twoPhaseProperties.thermo2().Cp()*(1.0-alpha1);
    volScalarField TCoeff = Hfg/rhoCp;

    volScalarField HdotRho = (vDotcT-vDotvT)*rho;

    volScalarField AbyV(mag(fvc::grad(limitedAlpha1)));  //Tanasawa
    volScalarField Cm1(2.0*CvTan*Hfg*rho2/((2.0-CvTan)*pow(2.0*M_PI*R,0.5)));  //Tanasawa


    if(cav_model == 4)

       { // Lee
           // temperature driven mass transfer term

	    forAll (p, celli) {

            if (T[celli] < TSat1[celli])     
            { // condensation

                    //Lee
                    vDotcT[celli] = TCoeff[celli]*(-Rc.value())*rho2[celli]*(T[celli]-TSat1[celli])/TSat1[celli]*(1.0 - limitedAlpha1[celli]); //Lee

            }

            else
            { // vaporization

                    //Lee
                    vDotvT[celli] = TCoeff[celli]*(Rv.value())*rho1[celli]*(T[celli] - TSat1[celli])/TSat1[celli]*limitedAlpha1[celli]; //Lee
            }}

        }//Lee


    if(cav_model == 5)

       { // Tanasawa
           // temperature driven mass transfer term

	    forAll (p, celli) {

		double tmp = TSat1[celli];

            if (T[celli] < TSat1[celli])         
            { // condensation

                    //Tanasawa
                    vDotcT[celli] = TCoeff[celli]*((-RcTan.value())*Cm1[celli]*AbyV[celli]/Foam::sqrt(Foam::pow(tmp,3.0))*(T[celli] - TSat1[celli])*(1.0-limitedAlpha1[celli]))/T_ref.value();  //Tanasawa

            }

            else
            { // vaporization

                    //Tanasawa
                   vDotvT[celli] = TCoeff[celli]*((RvTan.value())*Cm1[celli]*AbyV[celli]/Foam::sqrt(Foam::pow(tmp,3.0))*limitedAlpha1[celli]*(T[celli] - TSat1[celli]))/T_ref.value(); //Tanasawa

            }}

        }//Tanasawa


    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - fvm::laplacian(twoPhaseProperties.alphaEff(turbulence->mut()), T)
      + (
            fvc::div(fvc::absolute(phi, U), p)
          + fvc::ddt(rho, K) + fvc::div(rhoPhi, K)
        )
       *(
           alpha1/twoPhaseProperties.thermo1().Cv()
         + alpha2/twoPhaseProperties.thermo2().Cv()
        )
    );

    TEqn.relax();

    if(cav_model == 4 || 5)
    {
        solve
        (
                TEqn ==
                fvm::Sp(HdotRho,T) - (vDotcT-vDotvT)*TSat1*rho
        );

    }

    else
    {

        TEqn.solve();

    }

        Info << "TAve = "
        << T.weightedAverage(mesh.V()).value()
        << "  Min(T) = " << min(T).value()
        << "  Max(T) = " << max(T).value()
        << endl;

}
