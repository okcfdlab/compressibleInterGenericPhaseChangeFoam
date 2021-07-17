{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );

    surfaceScalarField phig
    (
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the fixedFluxPressure BCs to ensure flux consistency
     setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
            phiHbyA.boundaryField()
          - (mesh.Sf().boundaryField() & U.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    tmp<fvScalarMatrix> p_rghEqnComp1;
    tmp<fvScalarMatrix> p_rghEqnComp2;


    if (pimple.transonic())
    {
        surfaceScalarField phid1("phid1", fvc::interpolate(psi1)*phi);
        surfaceScalarField phid2("phid2", fvc::interpolate(psi2)*phi);

        p_rghEqnComp1 =
            fvc::ddt(rho1) + fvc::div(phi, rho1) - fvc::Sp(fvc::div(phi), rho1)
          + correction
            (
                psi1*fvm::ddt(p_rgh)
              + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
            );
        deleteDemandDrivenData(p_rghEqnComp1().faceFluxCorrectionPtr());
        p_rghEqnComp1().relax();

        p_rghEqnComp2 =
            fvc::ddt(rho2) + fvc::div(phi, rho2) - fvc::Sp(fvc::div(phi), rho2)
          + correction
            (
                psi2*fvm::ddt(p_rgh)
              + fvm::div(phid2, p_rgh) - fvm::Sp(fvc::div(phid2), p_rgh)
            );
        deleteDemandDrivenData(p_rghEqnComp2().faceFluxCorrectionPtr());
        p_rghEqnComp2().relax();
    }
    else
    {
        p_rghEqnComp1 =
            fvc::ddt(rho1) + psi1*correction(fvm::ddt(p_rgh))
          + fvc::div(phi, rho1) - fvc::Sp(fvc::div(phi), rho1);

        p_rghEqnComp2 =
            fvc::ddt(rho2) + psi2*correction(fvm::ddt(p_rgh))
          + fvc::div(phi, rho2) - fvc::Sp(fvc::div(phi), rho2);
    }

    // Cache p_rgh prior to solve for density update
    
    volScalarField p_rgh_0(p_rgh);

//////////////////////////////////////////////////////////////////////////////////////////////

    volScalarField limitedAlpha1(min(max(alpha1, scalar(0)), scalar(1)));
    volScalarField mcCoeff(Cc*rho2/tInf);
    volScalarField mvCoeff(Cv*rho2/(0.5*max(rho1,rhoMin)*sqr(UInf)*tInf));
    dimensionedScalar mcCoeffMerkle(Cc/(0.5*sqr(UInf)*tInf));   //Merkle
    volScalarField mvCoeffMerkle(Cv*rho1/(0.5*sqr(UInf)*tInf*max(rho2, rhoMin)));  //Merkle
    volScalarField AbyV(mag(fvc::grad(limitedAlpha1)));  //Tanasawa
    volScalarField Cm1(2.0*CvTan*Hfg*rho2/((2.0-CvTan)*pow(2.0*M_PI*R,0.5)));  //Tanasawa
    volScalarField pCoeff(1.0/max(rho1, rhoMin) - 1.0/max(rho2, rhoMin));
    volScalarField vDotHat((vDotcP-vDotvP)/(Cv*p_ref));

/////////////////////////////////////////////////////////////////////////////////////////////

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

////////////////////////////////////////////////////////////////////////////////////////////

    if(cav_model==0) { // no phase-change
	    volScalarField vDotvP(0.0*rho1);
	    volScalarField vDotcP(0.0*rho1);
    } else if(cav_model == 2) { // Kunz
           // pressure driven mass transfer term
	    forAll (p, celli) {
        if (p[celli] > pSat1[celli]) 
        { // condensation
          vDotcP[celli] = pCoeff[celli]*mcCoeff[celli]*sqr(limitedAlpha1[celli])*(1.0 - limitedAlpha1[celli])*(p[celli] - pSat1[celli])/(max(p[celli] - pSat1[celli], 0.01*pSat.value()));  //Kunz
        }
        else 
        { // vaporization
          vDotvP[celli] = pCoeff[celli]*(-mvCoeff[celli])*limitedAlpha1[celli]*(p[celli] - pSat1[celli]);  //Kunz
        }
	    } //forAll loop
    } else if(cav_model == 3) { // Merkle pressure driven mass transfer term
 	    forAll (p, celli) {
        if (p[celli] > pSat1[celli])     
        { // condensation
          vDotcP[celli] = pCoeff[celli]*mcCoeffMerkle.value()*(1.0 - limitedAlpha1[celli])*(p[celli] - pSat1[celli]);  //Merkle
        }
        else
        { // vaporization
          vDotvP[celli] = pCoeff[celli]*(-mvCoeffMerkle[celli])*limitedAlpha1[celli]*(p[celli] - pSat1[celli]);  //Merkle
          }
	    }//forAll loop
    } else if(cav_model == 4) { // Lee temperature driven mass transfer term
	  forAll (p, celli) {
            if (T[celli] < TSat1[celli])     
            { // condensation
                    vDotcP[celli] = pCoeff[celli]*(-Rc.value())*rho2[celli]*min(T[celli]-TSat1[celli],T0.value())/TSat1[celli]*(p[celli]-pSat1[celli])/max(p[celli]-pSat1[celli],1E-06*pSat.value())*(1.0 - limitedAlpha1[celli]);  //Lee
            } else
            { // vaporization
                    vDotvP[celli] = pCoeff[celli]*(-Rv.value())*rho1[celli]*max(T[celli]-TSat1[celli],T0.value())/TSat1[celli]*(p[celli]-pSat.value())/max(pSat.value()-p[celli],1E-06*pSat.value())*limitedAlpha1[celli];  //Lee
            }
        }//forAll loop	
    } else if(cav_model == 5)
       { // Tanasawa temperature driven mass transfer term
	forAll (p, celli) {
        double tmp = TSat1[celli];
   	    if (T[celli] < tmp)     
        { // condensation
          vDotcP[celli] = pCoeff[celli]*(-RcTan.value())*Cm1[celli]*min(T[celli] - TSat1[celli],T0.value())*AbyV[celli]/Foam::sqrt(Foam::pow(tmp,3.0))*(p[celli]-pSat.value())/max(p[celli]-pSat.value(),1E-6*pSat.value())*(1.0-limitedAlpha1[celli]);  //Tanasawa
        }
        else
        { // vaporization
          vDotvP[celli] = pCoeff[celli]*(-RvTan.value())*Cm1[celli]*max(T[celli] - TSat1[celli],T0.value())*AbyV[celli]/Foam::sqrt(Foam::pow(tmp,3.0))*(p[celli]-pSat.value())/max(pSat.value()-p[celli],1E-05*pSat.value())*limitedAlpha1[celli];  //Tanasawa
        }
	}
      }//Tanasawa

///////////////////////////////////////////////////////////////////////////////////////

        if(cav_model != 0)
        {
            solve
            (
                (
                    (max(alpha1, scalar(0))/rho1)*p_rghEqnComp1()
                  + (max(alpha2, scalar(0))/rho2)*p_rghEqnComp2()
                )

              + p_rghEqnIncomp == -(((vDotcP-vDotvP)*(pSat-rho*gh)/(Cv*p_ref)) - fvm::Sp(vDotHat,p_rgh)) ,
              mesh.solver(p_rgh.select(pimple.finalInnerIter()))

            );
        }

        else if(cav_model == 0)
        {
            solve
            (
                p_rghEqnIncomp,
                mesh.solver(p_rgh.select(pimple.finalInnerIter()))
            );

        } 


//////////////////////////////////////////////////////////////////////////////////////

        if (pimple.finalNonOrthogonalIter())
        {
            p = max(p_rgh + (alpha1*rho1 + alpha2*rho2)*gh, pMin);
            p_rgh = p - (alpha1*rho1 + alpha2*rho2)*gh;

            dgdt =
            (
                pos(alpha2)*(p_rghEqnComp2 & p_rgh)/rho2
              - pos(alpha1)*(p_rghEqnComp1 & p_rgh)/rho1
            );

            phi = phiHbyA + p_rghEqnIncomp.flux();

            U = HbyA
              + rAU*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
            U.correctBoundaryConditions();
        }

    }

    // Update densities from change in p_rgh
  
    rho1 += psi1*(p_rgh - p_rgh_0);
    rho2 += psi2*(p_rgh - p_rgh_0);

    rho = alpha1*rho1 + alpha2*rho2;

    K = 0.5*magSqr(U);

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p_rgh) " << min(p_rgh).value() << endl;

//////////////////////////////////////////////////////////////////////////////////////////
}
