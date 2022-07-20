#include "TMinuit.h"
#include "TRandom.h"
#include "LHAPDF/LHAPDF.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <utility> 
#include <map> //Error Matrix --- use map to store the entires --- then fill the matrix
#include "TComplex.h"
#include "TMatrixT.h"
#include <TCanvas.h> // . Added this and below T things for plotting. Some may not be necessary
#include <TGraph.h> 
#include <TLegend.h>
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TApplication.h"
using namespace std;
using namespace LHAPDF;

// USER INPUTS:
const Double_t RUNTIME_DAYS = 90.0;  // Number of days the experiment runs
const Double_t BEAM_POLARIZ  =0.85;  // Beam polarization
const Double_t SIN2_TH = 0.235;  // Sin^2(theta_w), (theta_w = weak mixing angle)
const string PDF_NAME = "CT18NLO";
Double_t BEAM_E = 11;//11;
const vector<string> PDF_names = {"CT18NLO", "CT18NLO"};//"PDF4LHC21_40"};  // Vector of  PDFs to analyze

// Create a new type that can return two values  // uod
struct Asym_Values{
  double Fit_Val;
  double Fixed_Val;
};
double r_prime=0.0;

// Constants
const Double_t GF = 1.16637e-5;  // GeV^2, Fermi constant
const Double_t ALPHA = 7.2973525e-3;  // Fine structure constant
const Double_t MP = 0.93828;  // GeV/c^2 Mass of a proton

const Int_t MAXBINS = 500;  // . changed to 2000 to work with the fine grid

//NOTE: The arrays are declared with fixed size, larger than what is needed for any file, to avoid overly complicated dynamic memory allocation
Double_t Gx[MAXBINS];
Double_t GQ2[MAXBINS];
Double_t Grate[MAXBINS]; // . Added this
Double_t Gu[MAXBINS];
Double_t Gubar[MAXBINS];
Double_t Gd[MAXBINS];
Double_t Gdbar[MAXBINS];
Double_t Gc[MAXBINS];
Double_t Gcbar[MAXBINS];
Double_t Gs[MAXBINS];
Double_t Gsbar[MAXBINS];
Double_t Gy[MAXBINS];
Double_t GY[MAXBINS];  // . added this

Double_t GA_clean[MAXBINS];  // Predicted A with no noise
Double_t GA_noisy[MAXBINS];  // Predicted A with noise added
// .  Some of these might not need to be global:
Double_t GA_cor_shift[MAXBINS];  // Pred. A w/ uncor. noise shifted by 1sigma cor
Double_t GA_clean_sim[MAXBINS];  // Predicted A ignoring sea, s, and c quarks
Double_t GA_fit[MAXBINS];  // Container for the fitted A value

// . Add the G for Global
// Global so can be used in error matrix
Double_t dA_stat[MAXBINS];  // Statistical error in A
Double_t dA_sys_uncor[MAXBINS];  // Uncorrelated systematic error in A
Double_t GdA_tot_uncor[MAXBINS];  // Total uncorrelated error (stat+uncor sys) in A
Double_t dA_sys_cor[MAXBINS];  // Correlated systematic error in A


Double_t Ggoodness[MAXBINS];  // . Use or delete

// . As time turn these to not global variables
Double_t GA_const[MAXBINS];  // . added this


Int_t start_i = 0;  // . added Values for when to start the next x
Int_t stop_i = 0;

// . Added below 3 lines
std::vector<std::vector<double>> du_fitted(sizeof(PDF_names), vector<double>(0));
std::vector<std::vector<double>> du_f_err(sizeof(PDF_names), vector<double>(0));
std::vector<std::vector<double>> du_x(sizeof(PDF_names), vector<double>(0));

std::vector<double> sinq2_val[6];  // . uod (use or delete)
std::vector<double> chi_square[6];  // . uod

double temporary_val =0.0;
//double PDF_O_D[90][90][90];
std::map<pair<int,int>, double> PDF_O_D;  // uod
//int counter = 0; // How many entries  // . uod
//PDF_O_D.clear();
ofstream Outfile;
int error_type;
int lum_sets=0;

// Prep for exporting data
ofstream all_bins_file;
ofstream x_cols_file;
ofstream fits_file;

void calcAsym(string fileName, string pdfName, double beamE) {
  cout << "PDFNAME=" << pdfName << endl;
  
  all_bins_file.open (PDF_NAME + "_" + to_string(beamE) + 
		      "GeV_PDFStuff_Analytic_Values_for_All_Bins.csv");
  all_bins_file << "x,Q2 (GeV^2),rate (Hz),A_clean,A_noisy,A_simplified,";
  all_bins_file << "dA_stat,dA_sys_uncor,dA_sys_cor,";
  all_bins_file << "d/u noisy,d/u clean,d/u PDF,d/u simplified,";
  all_bins_file << "d/u shifted by 1 sig sys cor,";
  all_bins_file << "d(d/u) from A_stat,d(d/u) from A_sys_uncor,";
  all_bins_file << "d(d/u) from A_sys_cor,d(d/u) from sin^2(theta_w),";
  all_bins_file << "d(d/u) from model\n";
  /*  all_bins_file << "a,b,c,d,g,f,ffffds\n";
  all_bins_file << "howdie\n";
  all_bins_file << "and, another one , bits, the ,dust";
  all_bins_file << "yeep, yep,\n";
  all_bins_file << "test, 2, 3,";
  all_bins_file << "and4";
  all_bins_file.close();*/

  x_cols_file.open (PDF_NAME + "_" + to_string(beamE) +
		    "GeV_PDFStuff_Analytic_Values_for_each_x.csv");
  x_cols_file << "x,Q2 (GeV^2),rate (Hz),";
  x_cols_file << "d/u noisy,d/u clean,d/u PDF,d/u simplified,";
  x_cols_file << "d/u shifted by 1 sig sys cor,";
  x_cols_file << "d(d/u) from A_stat,d(d/u) from A_sys_uncor,";
  x_cols_file << "d(d/u) from A_sys_cor,d(d/u) from sin^2(theta_w),";
  x_cols_file << "d(d/u) from model\n";


  Int_t i = 0;
  const PDF* pdf = mkPDF(pdfName, 0);
  fstream inputFile;
  cout << "fileName=" << fileName << endl;
  inputFile.open(fileName,fstream::in);
  
  // Set up random variables for use in generating pseudo-data
  TRandom *random = new TRandom();
  r_prime = random->Gaus();  // Random variable to be used with correlated systematic uncertainties

  // Get coupling constants
  double C_1d = 1/2. - 2 * SIN2_TH / 3;
  double C_1u = -1/2. + 4 * SIN2_TH / 3;
  double C_2d = 1/2. - 2 * SIN2_TH;
  double C_2u = -1/2. + 2 * SIN2_TH;
  //cout << "SIN2_TH=" <<  SIN2_TH << ", C_1d=" << C_1d;
  //cout << ", C_1u=" << C_1u << ", C_2d=" << C_2d << ", C_2u=";
  //cout << C_2u << endl;

  // Sums to be averaged when calculating analytic uncertainty for a value of x
  int num = 0;  // Number of bins for a given value of x
  double xSum = 0;  // x
  double Q2Sum = 0;  // Q2 (GeV^2)
  double rateSum = 0;  // rate (Hz)
  double duSum_noisy = 0;  // d/u calculated using noisy pseodo-data
  double duSum_clean = 0;  // d/u calculated using psuedo-data (no noise) 
  double duSum_exact = 0;  // d/u calculated using quark PDF values

  double weights = 0;
  double dduSum_A_sys_uncor = 0;
  double dduSum_sin2th = 0;

  double duSum_sys_cor_shift = 0;

  double duSum_sim = 0;

  // Read data from input file and use it (along with pdfs) to create pseudo-data
  string file_line;
  getline(inputFile, file_line);  // Remove header before looping
  // Loop through the input file one line at a time
  while (getline(inputFile, file_line)) {
    //cout << "TOP OF WHILE LOOP" << endl;
    //cout << "_______________________________________" << endl;

    // Read the current line one piece of data at a time.
    // Assumes file is formatted as follows: x,Q2,rate
    vector<string> tokens; // storage for tokens 
    string token; // current token    
    istringstream iss (file_line);
    while (getline(iss, token, ',')) {
      tokens.push_back(token);
    }
    Gx[i] = stof(tokens[0]);  // x (Bjorken scaling variable)
    GQ2[i] = stof(tokens[1]);  // GeV^2 
    Grate[i] = stof(tokens[2]);  // Hz (rate)

    // If just finished a value of x and are moving onto a new x,
    if (Gx[i] != Gx[i-1] and i != 0) {
      // Print out the average values for this x
      double xCom = xSum / weights;
      double Q2Com = Q2Sum / weights;
      double duCom_noisy = duSum_noisy / weights;
      double duCom_clean = duSum_clean / weights;
      double duCom_exact = duSum_exact / weights;
      cout << "x: " << xCom << ", Q2: " << Q2Com << ", rate:" << rateSum << endl;
      cout << "Combined noisy d/u: " << duCom_noisy << endl;
      cout << "Combined clean d/u: " << duCom_clean << endl;
      cout << "Combined exact d/u: " << duCom_exact << endl;
 
      double dduCom_A_stat = 1 / sqrt(weights);
      double dduCom_A_sys_uncor = dduSum_A_sys_uncor / weights;
      double dduCom_sin2th = dduSum_sin2th / weights;
      cout << "Combined d/u statistical uncertainty: ";
      cout << dduCom_A_stat << endl;
      cout << "Combined d/u A_sys_uncor uncertainty: "; 
      cout << dduCom_A_sys_uncor << endl;
      cout << "Combined d/u sin2(theta_w) uncertainty: ";
      cout << dduCom_sin2th << endl;

      double duCom_sys_cor_shift = duSum_sys_cor_shift / weights;
      double dduCom_A_sys_cor = duSum_noisy / weights - duCom_sys_cor_shift;
      cout << "Noisy d/u - Correlated Sys. Shifted d/u: ";
      cout << dduCom_A_sys_cor << endl;

      double duCom_sim = duSum_sim / weights;
      double dduCom_model = duCom_sim - duSum_clean / weights;
      cout << "Simplified_A d/u - clean_A d/u: ";
      cout << dduCom_model << endl;


      x_cols_file << xCom << "," << Q2Com << "," << rateSum;

      x_cols_file << "," << duCom_noisy << "," << duCom_clean;
      x_cols_file << "," << duCom_exact;

      x_cols_file << "," << duCom_sys_cor_shift << "," << duCom_sim;

      x_cols_file << "," << dduCom_A_stat << "," << dduCom_A_sys_uncor;
      x_cols_file << "," << abs(dduCom_A_sys_cor) << "," << dduCom_sin2th;
      x_cols_file << "," << abs(dduCom_model) << "\n";


      cout << "_____________________________________" << endl;
      // Reset sums to zero for the next x value
      xSum = 0;
      Q2Sum = 0;
      rateSum = 0;
      num = 0;
      duSum_noisy = 0;
      duSum_clean = 0;
      duSum_exact = 0;
      weights = 0;
      dduSum_A_sys_uncor = 0;
      dduSum_sin2th = 0;

      duSum_sys_cor_shift = 0;

      duSum_sim = 0;
    }
    
    // Obtain quark PDF values
    Gd[i] = (pdf->xfxQ2(1, Gx[i], GQ2[i]));// / Gx[i];  // d
    Gdbar[i] = (pdf->xfxQ2(-1, Gx[i], GQ2[i]));// / Gx[i];  // dbar
    Gu[i] = (pdf->xfxQ2(2, Gx[i], GQ2[i]));// / Gx[i];  // u
    Gubar[i] = (pdf->xfxQ2(-2, Gx[i], GQ2[i]));// / Gx[i];  // ubar
    Gs[i] = (pdf->xfxQ2(3, Gx[i], GQ2[i]));// / Gx[i];  // s
    Gsbar[i] = (pdf->xfxQ2(-3, Gx[i], GQ2[i]));// / Gx[i];  // sbar
    Gc[i] = (pdf->xfxQ2(4, Gx[i], GQ2[i]));// / Gx[i];  // c
    Gcbar[i] = (pdf->xfxQ2(-4, Gx[i], GQ2[i]));// / Gx[i];  // cbar
      
    // Obtain valance and plus quark PDF values
    double d_p = Gd[i] + Gdbar[i];
    double u_p = Gu[i] + Gubar[i];
    double s_p = Gs[i] + Gsbar[i];
    double c_p = Gc[i] + Gcbar[i];
    double d_V = Gd[i] - Gdbar[i];
    double u_V = Gu[i] - Gubar[i];
    double s_V = Gs[i] - Gsbar[i];
    double c_V = Gc[i] - Gcbar[i];
    
    // Obtain y and Y
    Gy[i] = GQ2[i] / (2 * MP * Gx[i] * beamE);  // Inelasticity
    GY[i] = (1 - (1-Gy[i]) * (1-Gy[i])) / (1 + (1-Gy[i]) * (1-Gy[i]));
      
    // Create theoretical A_{RL,p}^{e^{-},PVDIS}
    GA_const[i] = BEAM_POLARIZ * 3 * sqrt(2) * GF * GQ2[i] / 
      (4 * M_PI * ALPHA);  // Proportionality constant for A
    // Top left part of A's equation
    double A_tleft = 2 * (u_p + c_p) * C_1u - (d_p + s_p) * C_1d;
    // Top right part of A's equation
    double A_tright = GY[i] * (2 * (u_V + c_V) * C_2u - (d_V + s_V) * C_2d);
    // Divisor of A's equation
    double A_div = 4 * (u_p + c_p) + (d_p + s_p); 
    GA_clean[i] = GA_const[i] * (A_tleft + A_tright) / A_div;  // A_theory
    //    all_bins_file << "," << GA_clean[i];

    // Calculate uncertainties

    // Statistical uncertainty in A
    dA_stat[i] = 1 / sqrt(Grate[i] * RUNTIME_DAYS * 86400) / BEAM_POLARIZ;
    // Uncorrelated systematic uncertainty: 
    //   radiative corrections = 0.2%; event reconstruction = 0.2%
    dA_sys_uncor[i] = abs(GA_clean[i]) * sqrt(0.002*0.002 + 0.002*0.002);
    // Correlated systematic uncertainty: polarimetry = 0.4%, Q2 = 0.2%
    dA_sys_cor[i] = abs(GA_clean[i]) * sqrt(0.004*0.004 + 0.002*0.002); 
    // Total uncorrelated uncertainty
    GdA_tot_uncor[i] = sqrt(dA_stat[i]*dA_stat[i] + 
			  dA_sys_uncor[i]*dA_sys_uncor[i]);
    // Create the psuedo-data that has noise based on uncertainty
    GA_noisy[i] = GA_clean[i] + random->Gaus()*GdA_tot_uncor[i] 
      + r_prime*dA_sys_cor[i];
    //    all_bins_file << "," << GA_clean[i] << "," << GA_noisy[i]


    // d/u error stuff
    //cout << "x=" << Gx[i] << ", 
    cout << "Q2=" << GQ2[i];
    // Error in d/u due to statistical uncertainty in A 
    double ddu_A_stat = abs(dA_stat[i] * 
      ((-2)*GA_const[i] * (2*C_1d + C_1u + (2*C_2d + C_2u)*GY[i]) /
       pow((GA_clean[i] + C_1d*GA_const[i] + C_2d*GY[i]*GA_const[i]), 2)));
    cout << ", ddu_A_stat=" << ddu_A_stat;
    double weight = 1 / (ddu_A_stat*ddu_A_stat);
    weights += weight;

    // Error in d/u due to uncorrelated systematic uncertainty in A
    double ddu_A_sys_uncor = abs(dA_sys_uncor[i] * 
      ((-2)*GA_const[i] * (2*C_1d + C_1u + (2*C_2d + C_2u)*GY[i]) /
       pow((GA_clean[i] + C_1d*GA_const[i] + C_2d*GY[i]*GA_const[i]), 2)));
    cout << ", ddu_A_sys_uncor=" << ddu_A_sys_uncor;
    dduSum_A_sys_uncor += ddu_A_sys_uncor * weight;

    // Error in d/u due to uncertainty in SIN2_TH
    double ddu_sin2th = abs(0.0006 * 24*GA_const[i] * 
		    (GA_const[i] - 6*GA_clean[i]*GY[i] + GA_const[i]*GY[i]) /
                    pow(6*GA_clean[i] + 3*GA_const[i]*(1+GY[i]) -
                        4*GA_const[i]*SIN2_TH*(1+3*GY[i]), 2));
    cout << ", ddu_sin2th=" << ddu_sin2th << endl;
    dduSum_sin2th += ddu_sin2th * weight;


    // d/u stuff

    // Calculate d/u using pseudo-data that has noise
    Double_t du_numerator = -4*GA_noisy[i] + 2*C_1u*GA_const[i] + 
      2*C_2u*GY[i]*GA_const[i];
    Double_t du_denominator = GA_noisy[i] + C_1d*GA_const[i] + 
      C_2d*GY[i]*GA_const[i];
    double du_noisy = du_numerator/du_denominator;
    duSum_noisy += weight * du_noisy;

    // Calculate d/u using psuedo-data that has no noise
    du_numerator = -4*GA_clean[i] + 2*C_1u*GA_const[i] 
      + 2*C_2u*GY[i]*GA_const[i];
    du_denominator = GA_clean[i] + C_1d*GA_const[i] + C_2d*GY[i]*GA_const[i];
    double du_clean = du_numerator/du_denominator;
    duSum_clean += weight * du_clean;

    // Calculate the precise d/u value
    double du_exact = Gd[i]/Gu[i];
    duSum_exact += weight * du_exact;



    // Correlated uncertainty stuff

    // Shift the noisy data by 1 sigma of systematic correlated uncertainty
    GA_cor_shift[i] = GA_noisy[i] + dA_sys_cor[i];
    // Calculate d/u using the shifted pseudo-data
    du_numerator = -4*GA_cor_shift[i] + 2*C_1u*GA_const[i] 
      + 2*C_2u*GY[i]*GA_const[i];
    du_denominator = GA_cor_shift[i] + C_1d*GA_const[i] + C_2d*GY[i]*GA_const[i];
    double du_sys_cor_shift = du_numerator/du_denominator;
    duSum_sys_cor_shift += weight * du_sys_cor_shift;
    double ddu_sys_cor_shift = du_noisy - du_sys_cor_shift; // . Is this useful?



    // Checking uncertainty due to ignoring sea, s, and c quarks

    // Calculate A assuming ubar=dbar=sbar=cbar=s=c=0
    double A_const = BEAM_POLARIZ * 3 * sqrt(2) * GF * GQ2[i] / 
      (4 * M_PI * ALPHA);
    Double_t A_numerator = A_const *
      (2.0*C_1u + 2.0*GY[i]*C_2u - Gd[i]/Gu[i]*(C_1d + GY[i]*C_2d));
    Double_t A_denominator = Gd[i]/Gu[i] + 4;
    GA_clean_sim[i] = (A_numerator/A_denominator);
    // Calculate d/u using psuedo-data that has uncro noise and is shifted by cor
    du_numerator = -4*GA_clean_sim[i] + 2*C_1u*A_const + 2*C_2u*GY[i]*A_const;
    du_denominator = GA_clean_sim[i] + C_1d*A_const + C_2d*GY[i]*A_const;
    //cout << "d/u theory: " << du_numerator/du_denominator << endl;
    double du_sim = du_numerator/du_denominator;
    duSum_sim += weight * du_sim;
    double ddu_model = du_sim - du_clean;


    all_bins_file << Gx[i] << "," << GQ2[i] << "," << Grate[i];

    all_bins_file << "," << GA_clean[i] << "," << GA_noisy[i];
    all_bins_file << "," << GA_clean_sim[i];

    all_bins_file << "," << dA_stat[i] << "," << dA_sys_uncor[i]; 
    all_bins_file << "," << dA_sys_cor[i];

    all_bins_file << "," << du_noisy << "," << du_clean << "," << du_exact;
    all_bins_file << "," << du_sys_cor_shift << "," << du_sim;

    all_bins_file << "," << ddu_A_stat << "," << ddu_A_sys_uncor;
    all_bins_file << "," << ddu_sys_cor_shift << "," << ddu_sin2th;
    all_bins_file << "," << ddu_model << "\n";

    xSum += Gx[i] * weight;
    Q2Sum += GQ2[i] * weight;
    rateSum += Grate[i];
    num += 1;

    i++;
    //counter++;
  }
  double xCom = xSum / weights;
  double Q2Com = Q2Sum / weights;
  double duCom_noisy = duSum_noisy / weights;
  double duCom_clean = duSum_clean / weights;
  double duCom_exact = duSum_exact / weights;
  cout << "x: " << xCom << ", Q2: " << Q2Com << ", rate:" << rateSum << endl;
  cout << "Combined noisy d/u: " << duCom_noisy << endl;
  cout << "Combined clean d/u: " << duCom_clean << endl;
  cout << "Combined exact d/u: " << duCom_exact << endl;

  double dduCom_A_stat = 1 / sqrt(weights);
  double dduCom_A_sys_uncor = dduSum_A_sys_uncor / weights;
  double dduCom_sin2th = dduSum_sin2th / weights;
  cout << "Combined d/u statistical uncertainty: ";
  cout << dduCom_A_stat << endl;
  cout << "Combined d/u A_sys_uncor uncertainty: ";
  cout << dduCom_A_sys_uncor << endl;
  cout << "Combined d/u sin2(theta_w) uncertainty: ";
  cout << dduCom_sin2th << endl;

  double duCom_sys_cor_shift = duSum_sys_cor_shift / weights;
  double dduCom_A_sys_cor = duSum_noisy / weights - duCom_sys_cor_shift;
  cout << "Noisy d/u - Correlated Sys. Shifted d/u: ";
  cout << dduCom_A_sys_cor << endl;

  double duCom_sim = duSum_sim / weights;
  double dduCom_model = duCom_sim - duSum_clean / weights;
  cout << "Simplified_A d/u - clean_A d/u: ";
  cout << dduCom_model << endl;

  x_cols_file << xCom << "," << Q2Com << "," << rateSum;
  x_cols_file << "," << duCom_noisy << "," << duCom_clean;
  x_cols_file << "," << duCom_exact;
  x_cols_file << "," << duCom_sys_cor_shift << "," << duCom_sim;
  x_cols_file << "," << dduCom_A_stat << "," << dduCom_A_sys_uncor;
  x_cols_file << "," << abs(dduCom_A_sys_cor) << "," << dduCom_sin2th;
  x_cols_file << "," << abs(dduCom_model) << "\n";


  cout << "AFTER WHILE LOOP "<< endl;
  cout << "_______________________________________" << endl;

  // Clean up
  delete pdf;
  inputFile.close();
}
//______________________________________________________________________________
/* Asym_Values func(Double_t x, Double_t q2, Double_t y,
     Double_t u, Double_t ubar,Double_t d, Double_t dbar,Double_t s, Double_t sbar,Double_t c, Double_t cbar,
     Double_t *par)*/
// . Made return type double instead of blank
double func(Double_t *par, Int_t i) {

  // . Potentially delete new structure type Asym-Values or find a way to use it
  //Asym_Values Asym; // can access to parts of Asym
  
  // Coupling constants
  double C_1d = 1/2. - 2 * SIN2_TH / 3;
  double C_1u = -1/2. + 4 * SIN2_TH / 3;
  double C_2d = 1/2. - 2 * SIN2_TH;
  double C_2u = -1/2. + 2 * SIN2_TH;

  // Calculate A assuming ubar=dbar=sbar=cbar=s=c=0
  double A_const = BEAM_POLARIZ * 3 * sqrt(2) * GF * GQ2[i] / (4 * M_PI * ALPHA);
  Double_t A_numerator = A_const * 
    (2.0*C_1u + 2.0*GY[i]*C_2u - par[0]*(C_1d + GY[i]*C_2d));
  Double_t A_denominator = par[0] + 4;

  return (A_numerator/A_denominator); //return type needs to be of the form of double
  //return Asym;
}
//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {

  // Create and fill error matrix
  //Asym_Values Calculation;  // . Commented out. Use or delete this line
  int nbins = stop_i - start_i;  
  TMatrixT<double> Err_matrix(nbins, nbins);
  TMatrixT<double> Stat_Err_matrix(nbins, nbins);
  for (int i=0; i<nbins; i++) {
    for (int j=0; j<nbins; j++) {
      if (i==j) {
	Stat_Err_matrix(i,j)= (pow(dA_stat[i+start_i],2)
			       + pow(dA_sys_uncor[i+start_i],2) + pow(dA_sys_cor[i+start_i],2));
      }
      else{	 	 	
	Stat_Err_matrix(i,j) = 0.0; 
	// . should this equal pow(dA_sys_cor[i+start_i],2)?
	// This assumes that there are no off-diagonal elements
      }
    }
  }

  TMatrixT<double> Diff(1, nbins);
  TMatrixT<double> Diff_T(nbins, 1);
  Err_matrix = Stat_Err_matrix;

  //calculate chisquare
  //Double_t chisq = 0;
  Double_t delta;
  //Double_t Pseduo=0.0;
  //TRandom *random_2 = new TRandom();
  for(int i = start_i; i < stop_i; ++i) { 
    GA_fit[i] = func(par, i);  
    
    delta = (GA_noisy[i] - GA_fit[i]);///GdA_t_uncor[i];
    //delta = (Pseduo-GA_fit[i])/GdA_t_uncor[i]; 
    //chisq += (abs(Ggoodness[i] - 1.0000) < 0.1)*delta*delta;
    //Diff(0,i) = (abs(Ggoodness[i] - 1.0000) < 0.1)*pow((GA_noisy[i]-GA_fit[i]),1);
    //Diff(0,i) = (abs(Ggoodness[i] - 1.0000) < 0.1)*pow((Pseduo-GA_fit[i]),1); 
    // Diff(0,i) =pow((Pseduo-GA_fit[i]),1);
    Diff(0,i-start_i) = delta;//*delta;
  }

  TMatrixT<double> Invt(nbins, nbins);
  
  Err_matrix.SetTol(1.e-26);
  //double det = Err_matrix.Determinant(); 
  Double_t det_2 = 0.0;
  Invt = Err_matrix.Invert(&det_2);
  //Invt = Err_matrix.InvertFast(&det_2);
  Diff_T.Transpose(Diff);
  TMatrixT<double> Temp(nbins, nbins);
  TMatrixT<double> chi(1, nbins);
  //Temp.Mult(Diff,Err_matrix);
  Temp.Mult(Diff, Invt);  
  chi.Mult(Temp, Diff_T);
  f = chi.Sum();
}
//______________________________________________________________________________

int main(Int_t argc, char *argv[]) {

  if(argc < 1) {
    cout << "Please call with appropriate arguments:" << endl;
    cout << "./codeName beamType luminosity withUnfolding(true/false) inputFile";
    cout << endl;
    return -1;
  }
  //TCanvas* c1 = new TCanvas("c1", "c1", 800, 600); //.//
  //int Np = 0;//.//
  string fileName = argv[1];
  
  //.//string pdfName = "CT18NLO";
  //for (vector<string>::iterator t=PDF_names.begin(); t!=PDF_names.end(); ++t){
  //  string pdfName = *t;
  //for (auto &pdfName: PDF_names) {
  //for (int Np = 0; Np < PDF_names.size(); ++Np) {
  //string pdfName = PDF_names[Np];
  double beamE = BEAM_E; // GeV, beamline energy
  //string pdfName = "NNPDF31_nlo_as_0118";
  //string pdfName = "MMHT2014nlo68cl";
  cout << endl << "Simulation File Parameters:" << endl;
  cout << "File name: " << fileName << endl;
  cout << "Without unfolding" << endl;   // . What does this mean?
  cout << "Polarization: " << BEAM_POLARIZ << "%" << endl << endl;
  
  //calcAsym(fileName, pdfName, flagUnfolding, luminosity,beamTypeStr,Ordered_values);
  calcAsym(fileName, PDF_NAME, beamE);   // . was pdfName
  int i = 0;// . added the below loop and if statement and stop_i
  cout << "Gx[i]" << Gx[i]; 
  while (Gx[i] > 0) {
    //cout << endl << "Gx[i]=" << Gx[i] << ", i=" << i << ", Gx[i+1]=" << Gx[i+1];

    if (Gx[i] != Gx[i+1]) {
      stop_i = i + 1;
      //cout << ", stop_i=" << stop_i << ", start_i=" << start_i;
       
      // Instantiate Minuit for 1 parameter
      TMinuit *gMinuit = new TMinuit(1);//4); 
      gMinuit->SetFCN(fcn);  // Set the address of the function to be minimized
      
      //Perform fit with betas fixed
      //--------------------------------------------------------------------
      //cout << endl << "Betas fixed" << endl;
      Double_t arglist[10];
      Int_t ierflg = 0;
      
      // . I'm considering commenting out the next two lines because https://llr.in2p3.fr/activites/physique/glast/workbook/pages/advanced_GRUG/GRUGminuit.pdf  doesn't have them and I can't find a good explanation of what they do
      arglist[0] = 1;  // 1 because we want a chisquare fit
      gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
      
      static Double_t vstart[1] = {0.01};  // Initial parameter guess
      static Double_t step[1] = {0.001};
      // Set information for the parameter (like bounds)
      gMinuit->mnparm(0, "d/u", vstart[0], step[0], 0, 1, ierflg); // or 0, 0 to have no bounds
      
      // . I Commented out the next four lines because they seem like they don't need to be fixed, but I only kind of understand what is going on. I uncommented the 500 lines since it's in https://llr.in2p3.fr/activites/physique/glast/workbook/pages/advanced_GRUG/GRUGminuit.pdf
      arglist[0] = 500; // Max calls
      arglist[1] = .001;  // Tolerance  (was 1)
      //gMinuit->FixParameter(1);
      
      //gMinuit->SetPrintLevel(1);
      //double tmp_sin=0.0;
      //double tmp_error=0.0;
      //gMinuit -> GetParameter(0,tmp_sin,tmp_error);
      //cout<<"tmp_sin"<<","<<tmp_sin<<","<<tmp_error<<endl;
      // arglist is a list of arguments that MIGRAD takes. 2 means that it takes two arguments from arglist (I think based on https://root-forum.cern.ch/t/understanding-minuit/34836/5 )
      // Minimize the function
      gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
      //gMinuit->mnexcm("SCAN", arglist, 0, ierflg);   
      //cout<<"tmp_sin"<<","<<tmp_sin<<","<<tmp_error<<endl;
      //--------------------------------------------------------------------
      
      // Print results
      /*
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	//gMinuit->mnprin(3,amin);*/
      
      /*
	for(int i = 0; i < NBINS; ++i)
	cout << GA_clean[i] << "   " << GA_fit[i] << endl;
      */
      double du_out, err;
      gMinuit -> GetParameter(0, du_out, err);
      Outfile << du_out << "," << err << endl;
      du_fitted[0].push_back(du_out);
      du_f_err[0].push_back(err);
      du_x[0].push_back(Gx[i]);
      
      delete gMinuit;
    }
    i++;  // . Added this
    cout << "start_i=" << start_i << ", stop_i=" << stop_i;
    start_i = stop_i;  // . Added this
  }

  fits_file.open (PDF_NAME + "_" + to_string(beamE) + 
		  "GeV_PDFStuff_Fitted_Values_for_each_x.csv");
  fits_file << "x,d/u,d/u err";

  int length = du_x[0].size();
  cout << endl <<endl << du_x[0].size() << endl;
  for(int i=0; i < length; i++) {
    std::cout << " (du_x[0])[i]=" << (du_x[0])[i];
    std::cout << " (du_fitted[0])[i]=" << (du_fitted[0])[i];
    std::cout << " (du_f_err[0])[i]=" << (du_f_err[0])[i] << endl;
    fits_file << "\n" << (du_x[0])[i] << "," << (du_fitted[0])[i];
    fits_file << "," << (du_f_err[0])[i];
  }
  cout <<endl << endl;
  
  double du_x_a[du_x[0].size()];
  std::copy((du_x[0]).begin(), (du_x[0]).end(), du_x_a);
  double du_f_a[(du_fitted[0]).size()];
  std::copy((du_fitted[0]).begin(), (du_fitted[0]).end(), du_f_a);
  double du_f_err_a[(du_f_err[0]).size()];
  std::copy((du_f_err[0]).begin(), (du_f_err[0]).end(), du_f_err_a);
  double du_x_err_a[length];
  memset(du_x_err_a, 0, length*sizeof(int));

  /*
  cout << endl <<endl << du_x[0].size() << endl;
  for(int i=0; i < du_x[0].size(); i++) {
    std::cout << " (du_x_a)[i]=" << (du_x_a)[i];
    std::cout << " (du_f_a)[i]=" << (du_f_a)[i];
    std::cout << " (du_f_err[0])[i]=" << (du_f_err_a)[i];
    std::cout << " (du_x_err[0])[i]=" << (du_x_err_a)[i] << endl;
    }*/

  auto gr = new TGraphErrors((du_x[0]).size(), du_x_a, du_f_a, du_x_err_a, du_f_err_a);

  auto lineplot = new TGraph((du_x[0]).size(), du_x_a, du_f_a);
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);

  //gr->SetFillColor(2);
  gr->SetFillColorAlpha(2, 0.35);
  gr->GetXaxis()->SetTitle("x");
  gr->GetYaxis()->SetTitle("d/u");
  gr->SetTitle("");
  gr->Draw("a3");
  lineplot->SetMarkerStyle(7);
  lineplot->SetMarkerColorAlpha(2, 0.35); 
  lineplot->SetLineColor(2);
    lineplot->SetTitle("CT18NLO");//"Fitted d/u values");
  lineplot->Draw("same PL");
  //gPad->BuildLegend();
  auto legend = new TLegend(.1, .1, .3, .3, "legend");
  legend->AddEntry(lineplot,"CT18NLO","PL");
  legend->Draw("same");
  c1->Update();
  //string plot_name = "du_plot_" + PDF_NAME + "_" + std::to_string(BEAM_E) + "GeV.png";
  c1->SaveAs("du_plot.png");//plot_name);//"du_plot_" + PDF_NAME + "_" + std::to_string(BEAM_E) + "GeV.png");
  //  Np++
  //.//}
  /*
  vector<string> data={"Hello World!","Goodbye World!"};
  for (vector<string>::iterator t=data.begin(); t!=data.end(); ++t) 
    {
      cout<<*t<<endl;
    }
  */
  return 0;
}
