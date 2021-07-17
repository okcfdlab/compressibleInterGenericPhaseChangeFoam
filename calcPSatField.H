
{


const dimensionedScalar t30_11("30.11", dimensionSet(0,0,0,1,0,0,0), 30.11);
const dimensionedScalar t273_15("273.15", dimensionSet(0,0,0,1,0,0,0), 273.15);
const dimensionedScalar t1("1", dimensionSet(0,0,0,1,0,0,0), 1);
const dimensionedScalar p610_94("610.94", dimensionSet(1,-1,-2,0,0,0,0), 610.94);

const dimensionedScalar t443_028("443.028", dimensionSet(0,0,0,1,0,0,0), 443.028);
const dimensionedScalar t0_49("0.49", dimensionSet(0,0,0,1,0,0,0), 0.49);
const dimensionedScalar p_pa("1e5", dimensionSet(1,-1,-2,0,0,0,0), 1e5);


pSat1 == p_pa*pow(10, 3.9895 - (t443_028 / (T - t0_49)));


}


