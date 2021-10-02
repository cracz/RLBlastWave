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

    !!Don't use "normalization" and "amplitude" at the same time!!
    If fitting spectra, or using "amplitude", just fix "normalization" = 0.0 to turn it off.
    Only nomalization = 1.0 has any effect currently.

    nMax is the maximum for the summation in the emission function S(x,K).

    Index   Parameter
    -----   ---------
    par[0]   Ry (fm)
    par[1]   Rx (fm)
    par[2]   T (GeV)
    par[3]   rho_0
    par[4]   rho_2
    par[5]   tau (fm/c)
    par[6]   a_skin
    par[7]   del-tau (fm/c)
    par[8]   nMax (integer >= 1)
    par[9]   m0 (GeV rest mass)                 [FIXED]
    par[10]  iqstat (1 (boson) or -1 (fermion)) [FIXED]
    par[11]  normalization                      [FIXED]
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