#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "BlastWave.h"

void classTestMacro()
{
    gSystem->Load("libBlastWave");

    TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1000);
    canvas->SetGrid();
    canvas->SetTicks();
    canvas->SetLogy();
    canvas->SetLeftMargin(0.13);

    Double_t m0_pion = 0.140;   // Don't use too many decimal places or the integration doesn't work.
    Double_t m0_kaon = 0.494;
    Double_t m0_proton = 0.938;

    int boson = 1;
    int fermion = -1;

    /*
    Index   Parameter
    -----   ---------
    0       Ry (fm)
    1       Rx (fm)
    2       T (GeV)
    3       rho_0
    4       rho_2
    5       tau (fm/c)
    6       a_skin
    7       del-tau (fm/c)
    8       nMax (integer >= 1)
    9       m0 (GeV rest mass)
    10      iqstat (boson or fermion)
    */
   
    Double_t Ry = 12.04;//13.0;//12.04;
    Double_t Rx = 12.04;//11.0;//12.04;
    Double_t T  = 0.1;
    Double_t rho_0 = 0.9;
    Double_t rho_2 = 0.0;//0.05;//0.0;
    Double_t tau   = 9.0;
    Double_t a_skin = 0.0;
    Double_t del_tau = 2.0;
    Double_t nMax = 4.0;
    Double_t m0 = m0_proton;
    Double_t iqstat = (Double_t)fermion;
    Double_t norm = 1.0;    // Only 1 or 0
    Double_t amp = 1.0;     // Amplitude -- LEAVE THIS AT 1 FOR NOW

    BlastWave *blastWave = new BlastWave();
    
    Double_t pTStart = 0.0;
    Double_t pTStop  = 2.0;
    const int nParameters = 13;
    Double_t parameters[nParameters] = {Ry, Rx, T, rho_0, rho_2, tau, a_skin, del_tau, nMax, m0, iqstat, norm, amp};
    TF1 *BW = new TF1("BW",blastWave,&BlastWave::RLBlastWave,pTStart,pTStop,nParameters);
    BW->SetParameters(parameters);
    BW->FixParameter(9,m0);
    BW->FixParameter(10,iqstat);
    BW->FixParameter(11,norm);

    BW->SetParName(0,"R_{y} (fm)");
    BW->SetParName(1,"R_{x} (fm)");
    BW->SetParName(2,"T (GeV)");
    BW->SetParName(3,"#rho_{0}");
    BW->SetParName(4,"#rho_{2}");
    BW->SetParName(5,"#tau (fm/c)");
    BW->SetParName(6,"a_{s}");
    BW->SetParName(7,"#delta#tau (fm/c)");
    BW->SetParName(8,"N");
    BW->SetParName(9,"m_{0} (GeV/c^{2})");
    BW->SetParName(10,"iqstat");
    BW->SetParName(11,"Normalization");

    BW->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    BW->GetYaxis()->SetTitle("dN/p_{T}dp_{T} (arb. units)");
    BW->SetMinimum(10e-04);

    Double_t pionParameters[nParameters] = {Ry, Rx, T, rho_0, rho_2, tau, a_skin, del_tau, nMax, m0_pion, (Double_t)boson, norm, amp};
    TF1 *BW_pion = new TF1("BW_pion",blastWave,&BlastWave::RLBlastWave,pTStart,pTStop,nParameters);
    BW_pion->SetParameters(pionParameters);
    BW_pion->FixParameter(9,m0_pion);
    BW_pion->FixParameter(10,(Double_t)boson);
    BW_pion->FixParameter(11,norm);

    BW_pion->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    BW_pion->GetYaxis()->SetTitle("dN/p_{T}dp_{T} (arb. units)");
    BW_pion->SetMinimum(10e-04);

    BW->Draw();
    BW_pion->Draw("SAME");

    //canvas->SaveAs("BW_proton_N4_round.png");
    //canvas->SaveAs("BW_proton_N4_nonRound.png");
    //canvas->SaveAs("BW_pion_N1_round.png");
    //canvas->SaveAs("BW_pion_N1_nonRound.png");
    //canvas->SaveAs("BW_pion_N4_round.png");
    
    //canvas->SaveAs("BW_N4_nonRound.png");
    canvas->SaveAs("BW_N4_round.png");
    delete canvas;
    delete BW;
    delete BW_pion;
    delete blastWave;
}