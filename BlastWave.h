#ifndef BLASTWAVE_H
#define BLASTWAVE_H

// Usage instructions:
/*
    If you're using this in a ROOT macro make sure to both '#include "BlastWave.h"' and
    add this to your code to load the shared library libBlastWave.so:

    gSystem->Load("libBlastWave");

    In your code, make the class object first, then pass it and the address to the 
    RLBlastWave function to the TF1 constructor. Make sure you fix the necessary parameters:

    TF1 *BW = new TF1("BW",blastWave,&BlastWave::RLBlastWave,pTStart,pTStop,nParameters);
    BW->SetParameters(parameters);
    BW->FixParameter(9,m0);
    BW->FixParameter(10,iqstat);
    BW->FixParameter(11,norm);

    "nParameters" = 13
    "parameters" must be a vector of Double_t that corresponds to the list below.

    !! Don't use "normalization" and "amplitude" at the same time !!
    If fitting spectra, or using "amplitude", just fix "normalization" = 0.0 to turn it off.
    Only nomalization = 1.0 has any effect currently.

    !! Don't use too many decimal places for the rest mass 'm0' !!
    The integration doesn't work if there are too many decimal places, but three seems to be safe.

    Index   Parameter
    -----   ---------
    par[0]   Ry (fm)                Source transverse major axis radius.
    par[1]   Rx (fm)                Source transverse minor axis radius.
    par[2]   T (GeV)                Source temperature
    par[3]   rho_0                  Zeroth order term in source transverse flow rapidity.
    par[4]   rho_2                  Second order term in source transverse flow rapidity.
    par[5]   tau (fm/c)             Peak freezeout time in proper time.
    par[6]   a_skin                 *fractional* skin depth
    par[7]   del-tau (fm/c)         Duration of freezeout in proper time.
    par[8]   nMax (integer >= 1)    Max for the summation in the emission function.
    par[9]   m0 (GeV)           [FIXED]     Rest mass
    par[10]  iqstat             [FIXED]     1 (boson) or -1 (fermion)
    par[11]  normalization      [FIXED]
    par[12]  amplitude
*/


class BlastWave
{
private:
    /* data */
public:
    BlastWave(/* args */);
    ~BlastWave();
    static Double_t BW_R_Integrand(Double_t *x, Double_t *par);
    static Double_t BW_Phi_S_Integrand(Double_t *x, Double_t *par);
    static Double_t BW_Phi_P_Integrand(Double_t *x, Double_t *par);
    Double_t RLBlastWave(Double_t *x, Double_t *par);
};

#endif