{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir(phic*interface.nHatf());

    surfaceScalarField rhof("rhof", fvc::interpolate(rho));
    volScalarField alphalCoeff(1.0/max(rho1, rhoMin) - alpha1*(1.0/max(rho1, rhoMin) - 1.0/max(rho2,rhoMin)));
    volScalarField limitedAlpha1(min(max(alpha1, scalar(0)), scalar(1)));

    volScalarField mcCoeff(Cc*rho2/tInf);
    volScalarField mvCoeff(Cv*rho2/(0.5*max(rho1,rhoMin)*sqr(UInf)*tInf));

    dimensionedScalar mcCoeffMerkle(Cc/(0.5*sqr(UInf)*tInf));   
    volScalarField mvCoeffMerkle(Cv*rho1/(0.5*sqr(UInf)*tInf*max(rho2, rhoMin)));  

    volScalarField AbyV(mag(fvc::grad(limitedAlpha1)));  
    volScalarField Cm1(2.0*CvTan*Hfg*rho2/((2.0-CvTan)*pow(2.0*M_PI*R,0.5)));  
    

    if(cav_model==0) { // no phase-change

        volScalarField vDotcAlpha(0.0*rho1);
        volScalarField vDotcAlphal(0.0*rho1);
        volScalarField vDotvAlphal(0.0*rho1);

    } else if(cav_model == 2) { // Kunz

           // pressure driven mass transfer term
        forAll (p, celli) {
            if (p[celli] > pSat1[celli])        
            { // condensation

                vDotcAlphal[celli] = alphalCoeff[celli]*mcCoeff[celli]*sqr(limitedAlpha1[celli])*max(p[celli] - pSat1[celli], p0.value())/max(p[celli] - pSat1[celli], 0.01*pSat.value());  //Kunz

            }
            else
            { // vaporization

                vDotvAlphal[celli] = alphalCoeff[celli]*mvCoeff[celli]*min(p[celli] - pSat1[celli], p0.value());  //Kunz

            }
        }
    // Kunz
    } else if(cav_model == 3) { // Merkle
           // pressure driven mass transfer term
        forAll (p, celli) {
            if (p[celli] > pSat1[celli])   
            { // condensation

                vDotcAlphal[celli] = alphalCoeff[celli]*mcCoeffMerkle.value()*max(p[celli] - pSat1[celli], p0.value());  //Merkle

            }
            else
            { // vaporization

                vDotvAlphal[celli] = alphalCoeff[celli]*mvCoeffMerkle[celli]*min(p[celli] - pSat1[celli], p0.value());  //Merkle

            }
        }
        //Merkle
    } else if(cav_model == 4) { // Lee
           // temperature driven mass transfer term
        forAll (p, celli) {
            if (T[celli] < TSat1[celli])     
            { // condensation
                vDotcAlphal[celli] = alphalCoeff[celli]*(-Rc.value())*rho2[celli]*min(T[celli]-TSat1[celli],T0.value())/TSat1[celli];  //Lee
            }
            else
            { // vaporization
                vDotvAlphal[celli] = alphalCoeff[celli]*(-Rv.value())*rho1[celli]*max(T[celli]-TSat1[celli],T0.value())/TSat1[celli];  //Lee
            }
        }
        //Lee
    } else if(cav_model == 5) { // Tanasawa
           // temperature driven mass transfer term
        forAll (p, celli) {
            if (T[celli] < TSat1[celli])      
            { // condensation
                vDotcAlphal[celli] = alphalCoeff[celli]*(-RcTan.value())*Cm1[celli]*min(T[celli] - TSat1[celli], T0.value())*AbyV[celli]/Foam::sqrt(Foam::pow(TSat1[celli], 3.0));  //Tanasawa
            }
            else
            { // vaporization
                vDotvAlphal[celli] = alphalCoeff[celli]*(-RvTan.value())*Cm1[celli]*max(T[celli] - TSat1[celli], T0.value())*AbyV[celli]/Foam::sqrt(Foam::pow(TSat1[celli], 3.0));  //Tanasawa
            }
        }
    }//Tanasawa
    
    volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);

    tmp<surfaceScalarField> tphiAlpha;

    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
          - fvm::Sp(divU, alpha1)
         ==
           vDotcAlphal
          + mdot/max(rho, 0.01*rho1)

        );

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;

      tphiAlpha = alpha1Eqn.flux();
    }

    volScalarField alpha10("alpha10", min(max(alpha1, scalar(0.001)),scalar(1)));

    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
    {
        volScalarField::DimensionedInternalField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),

             dimensionedScalar("Sp", dgdt.dimensions(), 0.0) + vDotvAlphal - vDotcAlphal

        );

        volScalarField::DimensionedInternalField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution

            divU*alpha1
          + vDotcAlphal

        );

        forAll(dgdt, celli)
        {
            if (dgdt[celli] > 0.0 && alpha1[celli] > 0.0)
            {
                Sp[celli] -= dgdt[celli]*alpha1[celli];
                Su[celli] += dgdt[celli]*alpha1[celli];
            }
            else if (dgdt[celli] < 0.0 && alpha1[celli] < 1.0)
            {
                Sp[celli] += dgdt[celli]*(1.0 - alpha1[celli]);
            }
        }

        surfaceScalarField phiAlpha1
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha1,
            tphiAlpha,   //phi,
            phiAlpha1,
            Sp,
            Su,
            1,
            0
        );

        alpha1 = min(max(alpha1, scalar(0.0001)), scalar(1));

        surfaceScalarField rho1f(fvc::interpolate(rho1));
        surfaceScalarField rho2f(fvc::interpolate(rho2));
        rhoPhi = phiAlpha1*(rho1f - rho2f) + phi*rho2f;

        alpha2 = max(scalar(1) - alpha1, scalar(0.0001));
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Min(alpha2) = " << min(alpha2).value()
        << endl;
}
