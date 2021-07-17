    alphal = alpha1;
    alphav = alpha2;

    volScalarField alphalCoeff(1.0/max(rho1, rhoMin) - alpha1*(1.0/max(rho1, rhoMin) - 1.0/max(rho2,rhoMin)));
    volScalarField limitedAlpha1(min(max(alpha1, scalar(0)), scalar(1)));
    volScalarField mcCoeff(Cc*rho2/tInf);
    volScalarField mvCoeff(Cv*rho2/(0.5*max(rho1,rhoMin)*sqr(UInf)*tInf));

    dimensionedScalar mcCoeffMerkle(Cc/(0.5*sqr(UInf)*tInf));   //Merkle
    volScalarField mvCoeffMerkle(Cv*rho1/(0.5*sqr(UInf)*tInf*max(rho2, rhoMin)));  //Merkle

    volScalarField AbyV(mag(fvc::grad(limitedAlpha1)));  //Tanasawa
    volScalarField Cm1(2.0*CvTan*Hfg*rho2/((2.0-CvTan)*pow(2.0*M_PI*R,0.5)));  //Tanasawa

    scalar Rb, Rmin, min_divisor;
    
    if(cav_model==0)
      
       { // no phase-change
		
                 mdot.dimensionedInternalField()=0.0;    //
       }

    else if(cav_model == 2)
	{
	forAll(p, celli) 

       { // Kunz
           // pressure driven mass transfer term

            if (p.dimensionedInternalField()[celli] > pSat1.dimensionedInternalField()[celli])   
            { // condensation

                 mDotcAlphal.dimensionedInternalField()[celli] = mcCoeff.dimensionedInternalField()[celli]*sqr(limitedAlpha1.dimensionedInternalField()[celli])*max(p.dimensionedInternalField()[celli] - pSat1.dimensionedInternalField()[celli], p0.value())/max(p.dimensionedInternalField()[celli] - pSat1.dimensionedInternalField()[celli], 0.01*pSat.value());  //Kunz
                 mdot.dimensionedInternalField()[celli] = alphalCoeff.dimensionedInternalField()[celli]*mDotcAlphal.dimensionedInternalField()[celli]; //Kunz   

            }
            else
            { // vaporization

                 mDotvAlphal.dimensionedInternalField()[celli] = mvCoeff.dimensionedInternalField()[celli]*min(p.dimensionedInternalField()[celli] - pSat1.dimensionedInternalField()[celli], p0.value());  //Kunz
                 mdot.dimensionedInternalField()[celli] = alphalCoeff.dimensionedInternalField()[celli]*mDotvAlphal.dimensionedInternalField()[celli];  //Kunz  

            }

        }//Kunz
	}

    if(cav_model == 3)
	{
	forAll(p, celli)
       { // Merkle
           // pressure driven mass transfer term

            if (p.dimensionedInternalField()[celli] > pSat1.dimensionedInternalField()[celli])   
            { // condensation

                 mDotcAlphal.dimensionedInternalField()[celli] = mcCoeffMerkle.value()*max(p.dimensionedInternalField()[celli] - pSat1.dimensionedInternalField()[celli], p0.value()); //Merkle
                 mdot.dimensionedInternalField()[celli] = alphalCoeff.dimensionedInternalField()[celli]*mDotcAlphal.dimensionedInternalField()[celli]; //Merkle   

            }
            else
            { // vaporization

                 mDotvAlphal.dimensionedInternalField()[celli] = mvCoeffMerkle.dimensionedInternalField()[celli]*min(p.dimensionedInternalField()[celli] - pSat1.dimensionedInternalField()[celli], p0.value());  //M$
                 mdot.dimensionedInternalField()[celli] = alphalCoeff.dimensionedInternalField()[celli]*mDotvAlphal.dimensionedInternalField()[celli];  //Merkle  

            }

        }//Merkle
	}

    if(cav_model == 4)
	{
	forAll(p, celli)
       { // Lee
           // temperature driven mass transfer term

            if (T.dimensionedInternalField()[celli] < TSat1.dimensionedInternalField()[celli])     
            { // condensation

                      mDotcAlphal.dimensionedInternalField()[celli] = -Rc.value()*rho2.dimensionedInternalField()[celli]*min(T.dimensionedInternalField()[celli]-TSat1.dimensionedInternalField()[celli],T0.value())/TSat1.dimensionedInternalField()[celli];  //Lee
                      vDotcAlphal.dimensionedInternalField()[celli] = alphalCoeff.dimensionedInternalField()[celli]*mDotcAlphal.dimensionedInternalField()[celli]; //Lee

            }
            else
            { // vaporization


                      mDotvAlphal.dimensionedInternalField()[celli] = -Rv.value()*rho1.dimensionedInternalField()[celli]*max(T.dimensionedInternalField()[celli]-TSat1.dimensionedInternalField()[celli],T0.value())/TSat1.dimensionedInternalField()[celli];  //Lee
                      vDotvAlphal.dimensionedInternalField()[celli] = alphalCoeff.dimensionedInternalField()[celli]*mDotvAlphal.dimensionedInternalField()[celli];  //Lee

            }

        }//Lee
	}

    if(cav_model == 5)
	{
	forAll(p, celli)
       { // Tanasawa
           // temperature driven mass transfer term

            if (T.dimensionedInternalField()[celli] < TSat1.dimensionedInternalField()[celli])   
            { // condensation

                      mDotcAlphal.dimensionedInternalField()[celli] = -RcTan.value()*Cm1.dimensionedInternalField()[celli]*min(T.dimensionedInternalField()[celli] - TSat1.dimensionedInternalField()[celli], T0.value())*AbyV.dimensionedInternalField()[celli]/Foam::sqrt(Foam::pow(TSat1.dimensionedInternalField()[celli], 3.0));  //Tanasawa
                      vDotcAlphal.dimensionedInternalField()[celli] = alphalCoeff.dimensionedInternalField()[celli]*mDotcAlphal.dimensionedInternalField()[celli]; //Tanasawa

            }
            else
            { // vaporization


                      mDotvAlphal.dimensionedInternalField()[celli] = -RvTan.value()*Cm1.dimensionedInternalField()[celli]*max(T.dimensionedInternalField()[celli] - TSat1.dimensionedInternalField()[celli], T0.value())*AbyV.dimensionedInternalField()[celli]/Foam::sqrt(Foam::pow(TSat1.dimensionedInternalField()[celli], 3.0));  //Tanasawa
                      vDotvAlphal.dimensionedInternalField()[celli] = alphalCoeff.dimensionedInternalField()[celli]*mDotvAlphal.dimensionedInternalField()[celli];  //Tanasawa

            }

        }//Tanasawa
    }//forAll loop



