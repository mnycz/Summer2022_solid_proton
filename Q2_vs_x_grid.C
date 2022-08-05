//exit cuts: W > 2, trigger eff. cut
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include "math.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TPaveText.h"
#include "TText.h"
#include "TSystem.h"
#include "TArc.h"
#include "TString.h"
#include <vector>
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TFile.h"
#include "TMarker.h"

//#include "PVDIS_tri_eff.h"
#include "PVDIS_tri_eff_Rakitha1.h"
using namespace std;

const double runtime = 90.*86400; //LD2 90 days 11GeV
//const double runtime = 120.*86400; //LD2 120 days 6.6GeV
const double pb = 0.85; //beam polarization
bool Is_Rad = false; //whether we use the rate w/o ineternal and pre-vertex rad effects
                     //for eDIS generator, there is no internal and pre-vertex rad effects
                     //only the eAll generator has

const char* kKeyList[]={
  "eAll_commited4fe_rod_6mm_11_LH2_100files_1e6","eDIS_11G", "pDIS_11G", "eDIS_11G_bg", "pDIS_11G_bg",
    "eDIS_nobaffle_11G", "pDIS_nobaffle_11G", "eDIS_nobaffle_11G_bg", "pDIS_nobaffle_11G_bg",
    "eDIS_6.6G", "pDIS_6.6G", "eDIS_6.6G_bg", "pDIS_6.6G_bg",
    "eDIS_nobaffle_6.6G", "pDIS_nobaffle_6.6G", "eDIS_nobaffle_6.6G_bg", "pDIS_nobaffle_6.6G_bg",
    "_eAll_11G", "unknown"
};// "eAll_commited4fe_rod_6mm_11_LH2_100files_1e6" // 11 GeV
//"eAll_commitd6fc41a0_rod_6mm_22_LH_norad_1e6" // OLD 22 GeV
//"solid_PVDIS_LH2_moved_full_eAll_filenum100_22GeV_Z10cm_1e6" // Fixed 22 GeV

const int Nbin=15;
double bin[15][4]={
 0.20,0.30,     0.0,14.0,
 0.30,0.35,     0.0,14.0,
 0.35,0.40,     0.0, 5.8,
 0.35,0.40,     5.8,14.0,
 0.40,0.45,     0.0, 6.4,
 0.40,0.45,     6.4,14.0,
 0.45,0.50,     0.0, 7.0,
 0.45,0.50,     7.0,14.0,
 0.50,0.55,     0.0, 7.6,
 0.50,0.55,     7.6,14.0,
 0.55,0.60,     0.0, 8.2,
 0.55,0.60,     8.2,14.0,
 0.60,0.67,     0.0, 8.8,
 0.60,0.67,     8.8,14.0,
 0.67,0.80,     0.0,14.0
};

void TH2toTxt(TH2* h2, const char* key, bool skipZero=false)
{
    if(h2->GetEntries()<1.0E-4) {
        cout<<"TH2F "<<h2->GetName()<<"  has too few events, quit ...\n";
        return;
    }
    FILE * pFile;
    char buf[255];
    // . "/w/eic-scshelf2104/users/gsevans/thirdWeekSULIs22/files_11GeV/%s_%s.txt"
    pFile = fopen (Form("/w/eic-scshelf2104/users/gsevans/8thWeekSULIs22/11GeV_files/%s_%s.txt",h2->GetName(),key), "w");
    fprintf(pFile,"     x     Q2  rate(Hz)\n");
    TAxis *xAxis = h2->GetXaxis();
    TAxis *yAxis = h2->GetYaxis();
    int nx = xAxis->GetNbins();
    int ny = yAxis->GetNbins();

    for(int i=1;i<=nx;i++) {
        double xx=xAxis->GetBinCenter(i);
        for(int j=1;j<=ny;j++) {
            double yy=yAxis->GetBinCenter(j);
            double vv=h2->GetBinContent(i,j);
            if(skipZero && vv==0.0) continue;
            fprintf(pFile,"%6.3f %6.2f  %6.4f\n", xx, yy, vv);
        }
    }
    fclose (pFile);
}

void analysis_PVDIS_FOM_sim(const char* infile, const char* key, double beam=11.0, double current_uA=3.0)
{
  current_uA=50.0;
  gStyle->SetOptStat(0);
  cout << "GETTING FILE: " << infile << endl;
    TFile* f = new TFile(infile, "READ");
    cout << "GOTTEN FILE" << endl;
    TTree* t = (TTree*)f->Get("T");
    double Ei = 0;
    double Q2 = 0;
    double W = 0;
    double x = 0;
    double y = 0;
    double Abeam = 0;
    double px_gen = 0;
    double py_gen = 0;
    double pz_gen = 0;
    double vx_gen = 0;
    double vy_gen = 0;
    double vz_gen = 0;
    double p_gen, theta_gen, phi_gen;

    double rate = 0;
    double rateRad = 0;
    int    ecPID = 0;
    double ecPhi = 0;
    double ecR = 0;
    double ecP = 0;
    cout << "SETTING BRANCHES" << endl;
    t->SetBranchAddress("Ei",       &Ei    );
    t->SetBranchAddress("Q2",       &Q2     );
    t->SetBranchAddress("W",        &W      );
    t->SetBranchAddress("x",        &x      );
    t->SetBranchAddress("y",        &y      );
    t->SetBranchAddress("Abeam",    &Abeam  );

    t->SetBranchAddress("vx",       &vx_gen );
    t->SetBranchAddress("vy",       &vy_gen );
    t->SetBranchAddress("vz",       &vz_gen );
    t->SetBranchAddress("px",       &px_gen );
    t->SetBranchAddress("py",       &py_gen );
    t->SetBranchAddress("pz",       &pz_gen );
    t->SetBranchAddress("p",        &p_gen  );
    t->SetBranchAddress("theta", &theta_gen );
    t->SetBranchAddress("phi",     &phi_gen );

    t->SetBranchAddress("rate",     &rate   );
    t->SetBranchAddress("rateRad",  &rateRad);

    t->SetBranchAddress("ecPID",    &ecPID  );
    t->SetBranchAddress("ecPhi",    &ecPhi  );
    t->SetBranchAddress("ecR",      &ecR    );
    t->SetBranchAddress("ecP",      &ecP    );

    double Pi = atan(1) * 4;

    double thatrate[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double Abeam_sum[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double Q2_sum[15] = {0};
    double y_sum[15] = {0};

    int nBinQ2=14; // 14 for 11 GeV, 30 for 22 GeV
    double Q2Min=0.0,Q2Max=14; //14 for 11 GeV, 30 for 22 GeV
    double nBinx=10;  // 10 for analytic calculations, 100 for PDF
    cout << "ALMOST TO TH2F STUFF" << endl;
    if(beam>1.0 && beam<7.0) {
     nBinQ2=70;
     Q2Min=0.0;
     Q2Max=7.0;
    }
    //TH2F *hAbeamErr_Q2x_acc=new TH2F(Form("AbeamErr_Q2x_acc_%.0fuA",current_uA),"AbeamErr_Q2x_acc",100,0,1,nBinQ2,Q2Min,Q2Max);
    TH2F *hAbeamErr_Q2x_acc=new TH2F(Form("AbeamErr_Q2x_acc_%.0fuA",current_uA),"AbeamErr_Q2x_acc",nBinx,0,1,nBinQ2,Q2Min,Q2Max);
    TH2F *hAbeam_Q2x_acc=new TH2F(Form("Abeam_Q2x_acc_%.0fuA",current_uA),"Abeam_Q2x_acc",nBinx,0,1,nBinQ2,Q2Min,Q2Max);
    TH2F *hrate_Q2x_acc=new TH2F(Form("rate_Q2x_acc_finebin_%.0fuA",current_uA),Form("rate @ %.0fuA (Hz);x;Q^{2} [GeV^{2}]",current_uA),nBinx,0.0,1,nBinQ2,Q2Min,Q2Max);
    TH2F *hrate_Q2x_NoTrigEff=new TH2F(Form("rate_Q2x_finebin_NoTrigEff_%.0fuA",current_uA),Form("rate @ %.0fuA (Hz), No trigger eff. cut;x;Q^{2} [GeV^{2}]",current_uA),nBinx,0,1,nBinQ2,Q2Min,Q2Max);
    //this histogram just get the rate in Q2 and x bin
    TH2F* tmpRate_acc = new TH2F(Form("rate_Q2x_acc_%.0fuA",current_uA), Form("rate @ %.0f uA (Hz);x;Q^{2} [GeV^{2}]",current_uA), 10, 0, 1, nBinQ2/10,Q2Min,Q2Max);
    TH2F* tmpRate_NoTrigEff = new TH2F(Form("rate_Q2x_NoTrigEff_%.0fuA",current_uA), Form("rate @ %.0f uA (Hz), No trigger eff. cut;x;Q^{2} [GeV^{2}]",current_uA), 10, 0, 1, nBinQ2/10,Q2Min,Q2Max);
    cout << "PAST TH2F STUFF" << endl;
    for (unsigned int entry = 0; entry < t->GetEntries(); entry++) {
        t->GetEntry(entry);
	double theta_degree = theta_gen*(180.0/Pi);
        //if (W < 2) continue; //W < 2 GeV cut
	if (W>2.0 && theta_degree>22.0 && theta_degree<35.0) { 

	  //scale from originally 50 uA to 1 uA
	  double thisrate = rate*current_uA/50.0;
	  if (Is_Rad) thisrate = rateRad*current_uA/50.0;
	  
	  //get the trigger efficiency from EC
	  double eff = GetElectronTriggerEffi(GetRadiusIndex(ecR), GetMomentumIndex(ecP));
	  //for acceptance, we only consider whether the scattered e-/e+ hit the virtual
	  //plane in front of the EC or not. This has already been ensured in the fileReducer
	  //so acc is always 1 for now.
	  double acc = 1.;
	  
	  hAbeam_Q2x_acc->Fill(x,Q2,-Abeam*thisrate*acc*eff);
	  hrate_Q2x_acc->Fill(x,Q2,thisrate*acc*eff);
	  hrate_Q2x_NoTrigEff->Fill(x,Q2,thisrate*acc);
	  
	  tmpRate_acc->Fill(x,Q2,thisrate*acc*eff); //Get the rate in Hz
	  tmpRate_NoTrigEff->Fill(x,Q2,thisrate*acc); //Get the rate in Hz
	  
	  if (acc*eff !=0) {
	    //hAbeamErr_Q2x_acc->Fill(x,Q2,1./sqrt(thisrate*runtime*acc*eff)/pb*100);
	    hAbeamErr_Q2x_acc->Fill(x,Q2,(thisrate*runtime*acc*eff));
	  }
	  
	  for(int k = 0; k < Nbin; k++){
            if (bin[k][0] <= x && x < bin[k][1] && bin[k][2] <= Q2 && Q2 < bin[k][3]){
	      thatrate[k] += thisrate*acc*eff;
	      Abeam_sum[k] += -Abeam*thisrate*acc*eff;
	      Q2_sum[k] += Q2*thisrate*acc*eff;
	      y_sum[k]  += y*thisrate*acc*eff;
            }
	  }
	}
    }
    cout << "ALMOST DIVIDING" << endl;
    hAbeam_Q2x_acc->Divide(hAbeam_Q2x_acc,hrate_Q2x_acc);
    cout << "DIVIDED ONE" << endl;
    hAbeamErr_Q2x_acc->Divide(hAbeamErr_Q2x_acc,hAbeam_Q2x_acc);
    cout << "DONE DIVIDING" << endl;
    /* This section is commented out since it isn't needed when creating grids
    TCanvas *c_AbeamErr_Q2x_acc = new TCanvas("AbeamErr_Q2x_acc","AbeamErr_Q2x_acc",900,600);
    gPad->SetGrid();

    hAbeamErr_Q2x_acc->SetMarkerColor(kGreen);
    hAbeamErr_Q2x_acc->SetTitle("");
    hAbeamErr_Q2x_acc->GetXaxis()->SetTitle("x");
    hAbeamErr_Q2x_acc->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
    hAbeamErr_Q2x_acc->GetYaxis()->CenterTitle();
    hAbeamErr_Q2x_acc->GetXaxis()->CenterTitle();
    hAbeamErr_Q2x_acc->Draw();

    //relative statistical uncertainty of the parity violating asymmetry
    double AbeamErr[Nbin];
    cout << "y_ave"<<"\t"<<"Q2_ave" << "\t" << "Abeam_ave" << "\t" << "rate" << "\t" <<  "AbeamErr" << endl;
    for(int k = 0; k < Nbin; k++){
    //     double Apv = 0.84e-4*(bin[k][2]+bin[k][3])/2.;
        double Q2_ave = Q2_sum[k]/thatrate[k];
        double y_ave = y_sum[k]/thatrate[k];
        double Abeam_ave=Abeam_sum[k]/thatrate[k];
        AbeamErr[k] = 1./sqrt(thatrate[k]*runtime)/Abeam_ave/pb*100;
        cout << y_ave<<"\t"<<Q2_ave << "\t" << Abeam_ave << "\t" << int(thatrate[k]) << "\t" <<  AbeamErr[k] << " "<<thatrate[k]<<endl;
    }


    double x_cor[15]={0.250,0.325,0.375,0.375,0.425,0.425,0.475,0.475,0.525,0.525,0.575,0.575,0.635,0.635,0.735};
    double Q2_cor[15]={4.2,5.0,5.5,6.3,6.0,7.0,6.5,7.8,7.1,8.5,7.6,9.1,8.2,9.8,9.8};
    
    for(int k = 0; k < Nbin; k++){
        TMarker marker;
        marker.SetMarkerStyle(20);
        marker.SetMarkerColor(kRed);
        marker.DrawMarker(x_cor[k],Q2_cor[k]);
        TText *label = new TText(x_cor[k],Q2_cor[k],Form("%.02f",AbeamErr[k]));
        //TText *label = new TText(x_cor[k],Q2_cor[k],Form("%.02f",thatrate[k]/1000.));
        label->SetTextColor(kBlack);
        label->SetTextSize(0.03);
        label->Draw();
    }
    c_AbeamErr_Q2x_acc->SaveAs(Form("AbeamErr_Q2x_acc_%.0fuA_%s.png",current_uA,key));

    //this is the rate for different Q2 and x bin
    //set the draw-option "text" format
    gStyle->SetPaintTextFormat(".0f");
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    tmpRate_acc->Draw("colztext");
    gPad->SetLogz(1);
    c1->SaveAs(Form("rate_Q2x_acc_%.0fuA_%s.png",current_uA,key));

    TH2F* tmpRate_acc_copy = (TH2F*) tmpRate_acc->Clone();
    tmpRate_acc_copy->SetTitle("");
    tmpRate_acc_copy->Draw("colztext");
    c1->SaveAs(Form("rate_Q2x_acc_%.0fuA_%s.eps",current_uA,key));
    c1->SaveAs(Form("rate_Q2x_acc_%.0fuA_%s.C",current_uA,key));

    tmpRate_NoTrigEff->Draw("colztext");
    gPad->SetLogz(1);
    c1->SaveAs(Form("rate_Q2x_NoTrigEff_%.0fuA_%s.png",current_uA,key));

    TH2F* tmpRate_NoTrigEff_copy = (TH2F*) tmpRate_NoTrigEff->Clone();
    tmpRate_NoTrigEff_copy->SetTitle("");
    tmpRate_NoTrigEff_copy->Draw("colztext");
    c1->SaveAs(Form("rate_Q2x_NoTrigEff_%.0fuA_%s.eps",current_uA,key));
    c1->SaveAs(Form("rate_Q2x_NoTrigEff_%.0fuA_%s.C",current_uA,key));
    */
    cout << "ALMOST TO FILE PLACE" << endl;
    TH2toTxt(hrate_Q2x_acc,key,true);
    cout << "PAST FILE PLACE" << endl;
    ///TH2toTxt(hrate_Q2x_NoTrigEff,key,true);
    ///TH2toTxt(hAbeamErr_Q2x_acc,key,true);
}

void analysis_PVDIS_FOM_sim(int job, double beam, double current)
{
    char infile[255];
    //sprintf(infile,"../nt_PVDIS_%s.root",kKeyList[job]);
    sprintf(infile,"/lustre19/expphy/volatile/halla/triton/mnycz/SOLID/nt_PVDIS_%s.root",kKeyList[job]); 

    analysis_PVDIS_FOM_sim(infile,kKeyList[job], beam, current);
    
}

void ana()
{
  //analysis_PVDIS_FOM_sim(4,11,3.0);
  // analysis_PVDIS_FOM_sim(5,11,3.0);

  //analysis_PVDIS_FOM_sim(12,6.6,3.0);
  //analysis_PVDIS_FOM_sim(13,6.6,3.0);
}

void ana_baffle()
{
  analysis_PVDIS_FOM_sim(0,11,3.0); 
  // analysis_PVDIS_FOM_sim(1,11,3.0);

  //analysis_PVDIS_FOM_sim(8,6.6,3.0);
  // analysis_PVDIS_FOM_sim(9,6.6,3.0);
}

void Q2_vs_x_grid()//Analysis_Test_George()
{
  //ana();
    ana_baffle();
}

void testEi()
{
    analysis_PVDIS_FOM_sim("nt_PVDIS_eAll_Ei_11G.root","eAll_Ei_11G",11,1.0);
    analysis_PVDIS_FOM_sim("nt_PVDIS_eAll_11G.root","eAll_11G",11,1.0);
}
