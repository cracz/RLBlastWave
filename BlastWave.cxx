#include "TF1.h"
#include "TMath.h"
#include "BlastWave.h"

BlastWave::BlastWave(/* args */)
{
}

BlastWave::~BlastWave()
{
}

Double_t BlastWave::BW_R_Integrand(Double_t *x, Double_t *par)
{
    //=== Here x[0] = r ===//

    /*
    Index   Parameter
    -----   ---------
    par[0]   Ry (fm)
    par[1]   Rx (fm)
    par[2]   T (GeV)
    par[3]   rho_0
    par[4]   rho_2
    par[5]   tau (fm/c)
    par[6]   a_skin
    par[7]   del_tau (fm/c)
    par[8]   m0 (GeV rest mass)         [fixed]
    par[9]   pT (GeV/c)
    par[10]  n (summation iterator)     [fixed]
    par[11]  phi_p
    par[12]  phi_s
    */

    // (6) Sixth, create the full integrand
    Double_t r_ellipse = TMath::Sqrt(TMath::Power(x[0]*TMath::Cos(par[12])/par[1], 2) + TMath::Power(x[0]*TMath::Sin(par[12])/par[0], 2));

    Double_t spatial_density = 1.0;
    if (par[6] != 0.0) { spatial_density = 1.0 / (1.0 + TMath::Exp((r_ellipse - 1.0) / par[6])); }

    Double_t eta2 = TMath::Power(par[0]/par[1],2); // (Ry/Rx)^2
    Double_t phi_b = TMath::ATan2(TMath::Tan(par[12]), eta2);

    Double_t rho = r_ellipse * (par[3] + par[4]*TMath::Cos(2*phi_b));

    //Double_t pT = TMath::Sqrt(par[9]*par[9] - par[8]*par[8]);
    Double_t mT = TMath::Sqrt(par[8]*par[8] + par[9]*par[9]);
    
    Double_t alpha = (par[9]/par[2]) * TMath::SinH(rho);
    Double_t beta  = (mT/par[2]) * TMath::CosH(rho);

    Double_t fvalue = TMath::Exp(par[10] * alpha * TMath::Cos(phi_b - par[11])) * spatial_density;
    Double_t fvalue_r = x[0] * fvalue; // since the 'r' integral has 'rdr' in it

    Double_t beta_N = beta * par[10];
    Double_t K1  = TMath::BesselK(1,beta_N);
    Double_t G00 = 2.0 * K1;

    return fvalue_r * G00;
}

Double_t BlastWave::BW_Phi_S_Integrand(Double_t *x, Double_t *par)
{
    //=== Here x[0] = phi_s ===//

    /*
    Index   Parameter
    -----   ---------
    par[0]   Ry (fm)
    par[1]   Rx (fm)
    par[2]   T (GeV)
    par[3]   rho_0
    par[4]   rho_2
    par[5]   tau (fm/c)
    par[6]   a_skin
    par[7]   del_tau (fm/c)
    par[8]   m0 (GeV rest mass)
    par[9]   pT (GeV/c)
    par[10]  n (summation iterator)
    par[11]  phi_p
    */

    //                  Ry  >   Rx    ?   Ry   :   Rx
    Double_t Rmax = (par[0] > par[1]) ? par[0] : par[1];
    Rmax *= (1.0 + 8.0*par[6]);

    // (5) Fifth, take the 'r' integral of the full integrand from 0 to Rmax
    // Note phi_s is added:
    const int nParameters = 13;
    Double_t parameters[nParameters] = {par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8],par[9],par[10],par[11],x[0]};

    TF1 *r_integrand = new TF1("r_integrand",BW_R_Integrand,0.0,Rmax,nParameters);
    r_integrand->SetParameters(parameters);
    r_integrand->FixParameter(8,par[8]);
    r_integrand->FixParameter(10,par[10]);

    Double_t r_integrand_result = r_integrand->Integral(0, Rmax);
    delete r_integrand;
    return r_integrand_result;
}

Double_t BlastWave::BW_Phi_P_Integrand(Double_t *x, Double_t *par)
{
    //=== Here x[0] = phi_p ===//

    /*
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
    par[9]   m0 (GeV rest mass)
    par[10]  iqstat (boson or fermion)
    par[11]  pT (GeV/c)
    */
    
    Double_t nSum = 0;
    // (3) Third, this function has a for loop over the 'n' terms of the 'phi_s' and 'r' integrals
    for (int n = 1; n <= (int)par[8]; n++)
    {
        // (4) Fourth, take the 'phi_s' integral (for a single 'n') from 0 to 2pi
        Double_t signFactor = TMath::Power(par[10], n+1);

        // Note that nMax (par[8]) and iqstat (par[10]) are dropped, n (iterator) and phi_p (x[0]) are added:
        const int nParameters = 12;
        Double_t parameters[nParameters] = {par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[9],par[11],(Double_t)n,x[0]};

        TF1 *phi_s_integrand = new TF1("phi_s_integrand",BW_Phi_S_Integrand,0.0,TMath::TwoPi(),nParameters);
        phi_s_integrand->SetParameters(parameters);
        phi_s_integrand->FixParameter(8,par[9]);        // First argument is the NEW index in the 'parameters' array!
        phi_s_integrand->FixParameter(10,(Double_t)n);

        nSum += signFactor * phi_s_integrand->Integral(0, TMath::TwoPi());
        delete phi_s_integrand;
    }
    return nSum;
}

Double_t BlastWave::RLBlastWave(Double_t *x, Double_t *par)
{
    //=== Here x[0] = pT ===//

    /*
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
    par[9]   m0 (GeV rest mass)
    par[10]  iqstat (boson or fermion)
    par[11]  norm (normalization)
    par[12]  amp (amplitude)
    */

    Double_t mT = TMath::Sqrt(par[9]*par[9] + x[0]*x[0]);
    
    // (2) Second, create the integrand for the 'phi_p' integral and integrate it from 0 to 2pi for one particular pT
    const int nParameters = 12;
    Double_t parameters[nParameters] = {par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8],par[9],par[10],x[0]};
    
    TF1 *phi_p_integrand = new TF1("phi_p_integrand",BW_Phi_P_Integrand,0.0,TMath::TwoPi(),nParameters);

    Double_t normFactor = 1.0;
    if (par[11] == 1.0) 
    {
        // Get the value at pT = 0
        Double_t initialParameters[nParameters] = {par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8],par[9],par[10],0};
        phi_p_integrand->SetParameters(initialParameters);
        phi_p_integrand->FixParameter(9,par[9]);
        phi_p_integrand->FixParameter(10,par[10]);
        normFactor = par[9] * phi_p_integrand->Integral(0, TMath::TwoPi());
    }

    phi_p_integrand->SetParameters(parameters);
    phi_p_integrand->FixParameter(9,par[9]);
    phi_p_integrand->FixParameter(10,par[10]);

    Double_t phi_p_integral_result = phi_p_integrand->Integral(0, TMath::TwoPi());
    delete phi_p_integrand;

    return par[12] * mT * phi_p_integral_result / normFactor;
}
