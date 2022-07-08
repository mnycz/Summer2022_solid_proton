#include "TMinuit.h"
#include "TRandom.h"
#include "LHAPDF/LHAPDF.h"
#include<cmath>
#include<fstream>
#include <sstream>
#include <iostream>
#include <utility> 
#include <map> //Error Matrix --- use map to store the entires --- then fill the matrix
#include "TComplex.h"
#include "TMatrixT.h"
#include <TCanvas.h> // . Added this and below T things for plotting
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
const Double_t SIN2_TH = 0.235;  // Sin^2(theta_w), (theta_w is the weak mixing angle)
const vector<string> PDF_names = {"CT18NLO", "CT18NLO"};//"PDF4LHC21_40"};  // Vector of  PDFs to analyze

// Create a new type that can return two values
struct Asym_Values{
  double Fit_Val;
  double Fixed_Val;
};
double r_prime=0.0;

const Double_t GF = 1.16637e-5;  // GeV^2, Fermi constant
const Double_t ALPHA = 7.2973525e-3;  // Fine structure constant
const Double_t MP = 0.93828;  // GeV/c^2 Mass of a proton
const Double_t MN = 0.93957;
const Double_t MZ = 91.188;
const Int_t MAXBINS = 2000;  // . added the 2

const Double_t Qu = 2.0/3.0;
const Double_t Qc = 2.0/3.0;
const Double_t Qd = -1.0/3.0;
const Double_t Qs = -1.0/3.0;
const Double_t Qe = -1.0;
const Double_t I3u = 1.0/2.0;
const Double_t I3c = 1.0/2.0;
const Double_t I3d = -1.0/2.0;
const Double_t I3s = -1.0/2.0;
const Double_t I3e = -1.0/2.0;

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
Double_t GA_const[MAXBINS];  // . added this
Double_t Ga_rand[MAXBINS];
Double_t Gda[MAXBINS];
Double_t GdApvu_rel[MAXBINS];  // . use or delete
Double_t GdApv_rel[MAXBINS];  // . use or delete (check for other unused vars)
Double_t Ggoodness[MAXBINS];

Int_t start_i = 0;  // . added Values for when to start the next x
Int_t stop_i = 0;

// . Added below 5 lines
Double_t du_A_err[MAXBINS];  // Error in d/u due to A
Double_t du_sin2th_err[MAXBINS];  // Error in d.u due to sin^2(theta)
std::vector<std::vector<double>> du_fitted(sizeof(PDF_names), vector<double>(0));
std::vector<std::vector<double>> du_f_err(sizeof(PDF_names), vector<double>(0));
std::vector<std::vector<double>> du_x(sizeof(PDF_names), vector<double>(0));

Double_t GApv[MAXBINS];
Double_t GcalcA[MAXBINS];

Bool_t IS_ED = false;
Int_t NBINS = 0;

std::vector<double> sinq2_val[6];
std::vector<double> chi_square[6];

std::vector<double> Error_1; // make global ---- then use to create Error Matrix
std::vector<double> Error_Sys;
std::vector<double> Error_corr; //
std::vector<double> PDF_Diagonal;
std::vector<double> PDF_Off_Diagonal;
double temporary_val =0.0;
//double PDF_O_D[90][90][90];
std::map<pair<int,int>, double> PDF_O_D;
int counter = 0; // How many entries
//PDF_O_D.clear();
ofstream Outfile;
ofstream Outfile_chi;
ofstream Outfile_Sinsq;
ofstream Outfile_PDF;
int error_type;
int lum_sets=0;
//Outfile.open("/w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/Error_Comparison/Previous_Stat");
//Outfile<<"sin2"<<","<<"Uncertainty"<<endl;
//void calcAsym(string fileName, string pdfName, Bool_t flag_unf, Double_t lum, string beamType,Int_t Order)
void calcAsym(string fileName, string pdfName, double beamE) {
  cout << "PDFNAME=" << pdfName << endl;
  //ofstream Outfile; 
  // ****** I want to store the results is a csv file. Create Three directorys -- CT18, NNPDF, MMHT ---> and two subdirectories (Hydrogen and D2) ******//

  //change the location of the outfile to your directory. This will sort them based on name (if you want to store the results of the fit into a file)
  //     if(pdfName.compare("CT18NLO")==0){
  //Outfile.open(Form("/w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/Error_Comparison/CT18/%s/Updated_Stat_Sys_Diag.csv",Nuclei.Data()),std::ios_base::app); 
  //}
  //else if(pdfName.compare("NNPDF31_nlo_as_0118")==0){
  //Outfile.open(Form("/w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/Error_Comparison/NNPDF/%s/Updated_Stat_Sys_Diag.csv",Nuclei.Data()),std::ios_base::app);
  //}
  //else{
  //Outfile.open(Form("/w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/Error_Comparison/MMHT/%s/Updated_Stat_Sys_Diag.csv",Nuclei.Data()),std::ios_base::app);
  //}   
   
  Int_t i = 0;
  Double_t dA_stat, dA_sys_uncor, dA_sys_cor;

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
  double xsum = 0;  // x
  double Q2sum = 0;  // Q2
  double dusumR = 0;  // d/u calculated using noisy pseodo-data
  double dusumT = 0;  // d/u calculated using psuedo-data (no noise) 
  double dusumE = 0;  // d/u calculated using quark rates directly from PDFs
  double dusumR_err = 0;  // (Not yet used)
  double dusumT_err = 0;  // (Not yet used)
  double dusumE_err = 0;  // error in d/u calculated analytically. Not yet dependent on d/u
  double weights = 0;

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
    if (Gx[i] != Gx[i-1]) {
      // Print out the average values for this x
      cout << "x: " << xsum / num << ", Q2: " << Q2sum / num;
      cout << ", duSum_err: " << 1 / sqrt(weights) << endl;
      //cout << ", duSum_err: " << sqrt(dusumE_err) / num << endl;
      cout << "dusumR: " << dusumR / weights << ", dusumT: " << dusumT / weights;
      cout << ", dusumE: " << dusumE / weights << endl;
      //cout << "dusumR: " << dusumR / num << ", dusumT: " << dusumT / num;
      //cout << ", dusumE: " << dusumE / num << endl;
      // Reset sums to zero for the next x value
      xsum = 0;
      Q2sum = 0;
      num = 0;
      //dusumE_err = 0;
      dusumR = 0;
      dusumT = 0;
      dusumE = 0;
      weights = 0;
    }
    
    // Obtain quark rates
    Gd[i] = (pdf->xfxQ2(1, Gx[i], GQ2[i]));// / Gx[i];  // d
    Gdbar[i] = (pdf->xfxQ2(-1, Gx[i], GQ2[i]));// / Gx[i];  // dbar
    Gu[i] = (pdf->xfxQ2(2, Gx[i], GQ2[i]));// / Gx[i];  // u
    Gubar[i] = (pdf->xfxQ2(-2, Gx[i], GQ2[i]));// / Gx[i];  // ubar
    Gs[i] = (pdf->xfxQ2(3, Gx[i], GQ2[i]));// / Gx[i];  // s
    Gsbar[i] = (pdf->xfxQ2(-3, Gx[i], GQ2[i]));// / Gx[i];  // sbar
    Gc[i] = (pdf->xfxQ2(4, Gx[i], GQ2[i]));// / Gx[i];  // c
    Gcbar[i] = (pdf->xfxQ2(-4, Gx[i], GQ2[i]));// / Gx[i];  // cbar
      
    // Obtain valance and plus quark rates
    double d_p = Gd[i] + Gdbar[i];
    double u_p = Gu[i] + Gubar[i];
    double s_p = Gs[i] + Gsbar[i];
    double c_p = Gc[i] + Gcbar[i];
    double d_V = Gd[i] - Gdbar[i];
    double u_V = Gu[i] - Gubar[i];
    double s_V = Gs[i] - Gsbar[i];
    double c_V = Gc[i] - Gcbar[i];
    //cout << "d_p=" <<  d_p << ", u_p=" << u_p; 
    //cout << ", s_p=" << s_p << ", c_p=" << c_p << endl;
    
    // Obtain y and Y
    Gy[i] = GQ2[i] / (2 * MP * Gx[i] * beamE);  // Inelasticity
    GY[i] = (1 - (1-Gy[i]) * (1-Gy[i])) / (1 + (1-Gy[i]) * (1-Gy[i]));
    //cout << "Gy[i]=" <<  Gy[i] << ", GY[i]=" << GY[i] << endl;
      
    // Create theoretical A_{RL,p}^{e^{-},PVDIS}
    GA_const[i] = BEAM_POLARIZ * 3 * sqrt(2) * GF * GQ2[i] / 
      (4 * M_PI * ALPHA);  // Proportionality constant for A
    // Top left part of A's equation
    double A_tleft = 2 * (u_p + c_p) * C_1u - (d_p + s_p) * C_1d;
    // Top right part of A's equation
    double A_tright = GY[i] * (2 * (u_V + c_V) * C_2u - (d_V + s_V) * C_2d);
    // Divisor of A's equation
    double A_div = 4 * (u_p + c_p) + (d_p + s_p); 
    //cout << "GA_const=[i]" << GA_const[i] << ", A_tleft=" << A_tleft;
    //cout << ", A_tright=" << A_tright << ", A_div=" << A_div << endl;
    GApv[i] = GA_const[i] * (A_tleft + A_tright) / A_div;  // A_theory
    //cout << "GApv=" << GApv[i] << ", GQ2[i]=" << GQ2[i];
    //cout << ", Gx[i]=" << Gx[i] << endl << endl;

    // Calculate uncertainties

    // Statistical uncertainty in A
    dA_stat = 1 / sqrt(Grate[i] * RUNTIME_DAYS * 86400) / BEAM_POLARIZ;
    Error_1.push_back(dA_stat);

    // Uncorrelated systematic uncertainty: 
    //   radiative corrections = 0.2%; event reconstruction = 0.2%
    dA_sys_uncor = abs(GApv[i]) * sqrt(0.002*0.002 + 0.002*0.002);
    Error_Sys.push_back(dA_sys_uncor);

    // Correlated systematic uncertainty: polarimetry = 0.4%, Q2 = 0.2%
    dA_sys_cor = abs(GApv[i]) * sqrt(0.004*0.004 + 0.002*0.002); 
    Error_corr.push_back(dA_sys_cor);

    // Total uncorrelated uncertainty
    Gda[i] = sqrt(dA_stat*dA_stat + dA_sys_uncor*dA_sys_uncor);
    //cout << "dA_stat=" << dA_stat << ", dA_sys_uncor=" << dA_sys_uncor << endl;
    //cout << "Gda[i]=" << Gda[i] << ", dA_sys_cor=" << dA_sys_cor << endl;
      
    // Create the psuedo-data that has noise based on uncertainty
    Ga_rand[i] = GApv[i] + random->Gaus()*Gda[i] + r_prime*dA_sys_cor;

    // d/u error stuff

    // Error in d/u due to statistical uncertainty in A  // . change to du_A_stat
    /*du_A_err[i] =
      abs(Gda[i] * (-2)*GA_const[i] * (2*C_1d + C_1u + (2*C_2d + C_2u)*GY[i]) /
          pow((GApv[i] + C_1d*GA_const[i] + C_2d*GY[i]*GA_const[i]), 2));
    */
    double du_A_stat = abs(dA_stat) * 
      ((-2)*GA_const[i] * (2*C_1d + C_1u + (2*C_2d + C_2u)*GY[i]) /
       pow((GApv[i] + C_1d*GA_const[i] + C_2d*GY[i]*GA_const[i]), 2));
    double du_A_sys_uncor = abs(dA_sys_uncor) * 
      ((-2)*GA_const[i] * (2*C_1d + C_1u + (2*C_2d + C_2u)*GY[i]) /
       pow((GApv[i] + C_1d*GA_const[i] + C_2d*GY[i]*GA_const[i]), 2));
    du_A_err[i] = sqrt(du_A_stat*du_A_stat + du_A_sys_uncor*du_A_sys_uncor);

    //cout << "du_A_err[i]=" << du_A_err[i];

    // Error in d/u due to uncertainty in SIN2_TH
    du_sin2th_err[i] = abs(0.0006 * 24*GA_const[i] *
                           (GA_const[i] - 6*GApv[i]*GY[i] + GA_const[i]*GY[i]) /
                           pow(6*GApv[i] + 3*GA_const[i]*(1+GY[i]) -
                               4*GA_const[i]*SIN2_TH*(1+3*GY[i]), 2));
    //cout << ", du_sin2th_err[i]=" << du_sin2th_err[i] << endl;
    //double du_err = 1/(du_A_err[i]*du_A_err[i] + du_sin2th_err[i]*du_sin2th_err[i]);
    //weights += du_err;
    //dusumE_err += du_A_err[i]*du_A_err[i] + du_sin2th_err[i]*du_sin2th_err[i];
    double weight = 1 / 
      (du_A_err[i]*du_A_err[i] + du_sin2th_err[i]*du_sin2th_err[i]);
    weights += weight;

    // d/u stuff

    // Calculate d/u using pseudodata that has noise
    Double_t du_numerator = -4*Ga_rand[i] + 2*C_1u*GA_const[i] + 
      2*C_2u*GY[i]*GA_const[i];
    Double_t du_denominator = Ga_rand[i] + C_1d*GA_const[i] + C_2d*GY[i]*GA_const[i];
    //cout << "d/u rand: " << du_numerator/du_denominator << endl;
    //dusumR += du_numerator/du_denominator;  //
    dusumR += weight*du_numerator/du_denominator;

    // Calculate d/u using psuedodata that has no noise
    du_numerator = -4*GApv[i] + 2*C_1u*GA_const[i] + 2*C_2u*GY[i]*GA_const[i];
    du_denominator = GApv[i] + C_1d*GA_const[i] + C_2d*GY[i]*GA_const[i];
    //cout << "d/u theory: " << du_numerator/du_denominator << endl;
    dusumT += weight*du_numerator/du_denominator;

    // Calculate the precise d/u value
    //cout << "d/u: " << Gd[i]/Gu[i] << endl;
    dusumE += weight*Gd[i]/Gu[i];

    xsum += Gx[i];
    Q2sum += GQ2[i];
    num += 1;
    //dusumE_err += du_A_err[i]*du_A_err[i] + du_sin2th_err[i]*du_sin2th_err[i]; 

    i++;
    counter++;
  }
  // Print out the averaged values for the last x value
  cout << "x: " << xsum / num << ", Q2: " << Q2sum / num;
  cout << ", duSum_err: " << 1 / sqrt(weights) << endl;
  cout << "dusumR: " << dusumR / weights << ", dusumT: " << dusumT / weights;
  cout << ", dusumE: " << dusumE / weights << endl;
  //  cout << "x: " << xsum / num << ", Q2: " << Q2sum / num;
  //cout << ", duSum_err: " << sqrt(dusumE_err) / num << endl;
  //cout << "dusumR: " << dusumR / num << ", dusumT: " << dusumT / num;
  //cout << ", dusumE: " << dusumE / num << endl;  

  cout << "AFTER WHILE LOOP "<< endl;
  cout << "_______________________________________" << endl;

  NBINS = counter;
  //cout << "NBINS: " << NBINS;

  // Clean up
  delete pdf;
  inputFile.close();
}
//______________________________________________________________________________
/* Asym_Values func(Double_t x, Double_t q2, Double_t y,
     Double_t u, Double_t ubar,Double_t d, Double_t dbar,Double_t s, Double_t sbar,Double_t c, Double_t cbar,
     Double_t *par)*/
// . Made return type double instead of blank
double func(Double_t x, Double_t q2, Double_t y, Double_t u, Double_t ubar, Double_t d, Double_t dbar, Double_t s, Double_t sbar, Double_t c, Double_t cbar, Double_t *par, Int_t i) {

  // . Potentially delete new structure type Asym-Values or find a way to use it
  //Asym_Values Asym; // can access to parts of Asym

  //Double_t uplus = u + ubar;
  //Double_t uv = u - ubar;
  //Double_t dplus = d + dbar;
  //Double_t dv = d - dbar;
  //Double_t cplus = c + cbar;
  //Double_t splus = s + sbar;

  //double Y = (1 - (1-Gy[i]) *  (1-Gy[i])) / (1 +  (1-Gy[i]) *  (1-Gy[i]));
  
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
  counter = stop_i - start_i;  // . ADDED THIS LINE
  //cout << "counter=" << counter << ", stop_i=" << stop_i;
  //cout << ", start_i=" << stop_i << endl;
  TMatrixT<double> Err_matrix(counter, counter);
  TMatrixT<double> Stat_Err_matrix(counter, counter);
  for (int i=0; i<counter; i++) {
    for (int j=0; j<counter; j++) {
      if (i==j) {
	//Stat_Err_matrix(i,j)=pow(Error_1[i],2);
	Stat_Err_matrix(i,j)= (pow(Error_1[i+start_i],2)
			       + pow(Error_Sys[i+start_i],2) + pow(Error_corr[i+start_i],2));
	//Stat_Err_matrix(i,j)= (pow(Error_1[i],2)  + pow(Error_Sys[i],2));
	//cout << i << " " << j;
      }
      else{	 	 	
	Stat_Err_matrix(i,j) = 0.0; 
	// . should this equal pow(Error_corr[i+start_i],2)?
	// This assumes that there are no off-diagonal elements
      }
    }
  }

  //Stat_Err_matrix.Print();
  //cout<<"Check Stat Matrix"<<","<<Stat_Err_matrix(0,0)<<","
  //cout<<Stat_Err_matrix(1,1) <<endl;

  TMatrixT<double> Diff(1, counter);
  TMatrixT<double> Diff_T(counter, 1);
  Err_matrix = Stat_Err_matrix;

  //calculate chisquare
  //Double_t chisq = 0;
  Double_t delta;
  //Double_t Pseduo=0.0;
  //TRandom *random_2 = new TRandom();
  for(int i = start_i; i < stop_i; ++i) { // . Replaced 0 w/ start_i and NBINS w/ stop_i
    GcalcA[i] = func(Gx[i], GQ2[i], Gy[i], Gu[i], Gubar[i], Gd[i], Gdbar[i], Gs[i], Gsbar[i], Gc[i], Gcbar[i], par, i);  // . added the i
    //Calculation=func(Gx[i],GQ2[i],Gy[i],Gu[i],Gubar[i],Gd[i],Gdbar[i],Gs[i],Gsbar[i],Gc[i],Gcbar[i],par);  
    //GcalcA[i]=Calculation.Fit_Val;
    //Pseduo= Calculation.Fixed_Val + random_2->Gaus()*Gda[i] + r_prime*(sqrt(pow(0.01*Calculation.Fit_Val*1,2)));
   
    //cout << "" << GcalcA[i] << "Gx[i]=" << Gx[i] << ", GQ2[i]" << GQ2[i] << " par" << par[0] << endl;//  << ", " << Gx[i] << ", " <<  GQ2[i] << ", " <<  Gy[i] << ", " <<  Gu[i] << ", " <<  Gubar[i] << ", " <<  Gd[i] << ", " <<  Gdbar[i] << ", " <<  Gs[i] << ", " <<  Gsbar[i] << ", " <<  Gc[i] << ", " <<  Gcbar[i] << ", " << i << endl;
    //cout << Ga_rand[i] << ", " << GcalcA[i] << ", " << Gda[i] << "A terms" <<endl;
    delta = (Ga_rand[i] - GcalcA[i]);///Gda[i];
    //delta = (Pseduo-GcalcA[i])/Gda[i]; 
    //chisq += (abs(Ggoodness[i] - 1.0000) < 0.1)*delta*delta;
    //Diff(0,i) = (abs(Ggoodness[i] - 1.0000) < 0.1)*pow((Ga_rand[i]-GcalcA[i]),1);
    //Diff(0,i) = (abs(Ggoodness[i] - 1.0000) < 0.1)*pow((Pseduo-GcalcA[i]),1); 
    // Diff(0,i) =pow((Pseduo-GcalcA[i]),1);
    Diff(0,i-start_i) = delta;//*delta;
    //cout << "delta " << delta;
    //cout << " Diff" << Diff(0,i) << endl;
  }
  //Diff.Print();
  TMatrixT<double> Invt(counter, counter);
  
  Err_matrix.SetTol(1.e-26);
  //double det = Err_matrix.Determinant(); 
  // cout << "det: " << det << endl;
  Double_t det_2 = 0.0;
  Invt = Err_matrix.Invert(&det_2);
  //Invt = Err_matrix.InvertFast(&det_2);
  //cout<<"Valid"<<","<<Invt.IsValid()<<","<<det<<","<<det_2<<endl;
  //cout << "Invt:";
  //Invt.Print();
  Diff_T.Transpose(Diff);
  TMatrixT<double> Temp(counter, counter);
  TMatrixT<double> chi(1, counter);
  //Temp.Mult(Diff,Err_matrix);
  Temp.Mult(Diff, Invt);  
  chi.Mult(Temp, Diff_T);
  //TMatrix chi = Temp.Mult
  //chisq += (Diff
  //cout << "Temp: ;";
  //Temp.Print();
  //cout << "chi:";
  //chi.Print();
  // .chi_square[lum_sets].push_back(chi.Sum());
  //cout<<"chisq"<<","<<chisq<<endl;
  //}
  //f = chisq;
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
  int Np = 0;//.//
  string fileName = argv[1];
  
  string pdfName = "CT18NLO";
  //for (vector<string>::iterator t=PDF_names.begin(); t!=PDF_names.end(); ++t){
  //  string pdfName = *t;
  //for (auto &pdfName: PDF_names) {
  //for (int Np = 0; Np < PDF_names.size(); ++Np) {
  //string pdfName = PDF_names[Np];
  double beamE = 11; // GeV, beamline energy
  //string pdfName = "NNPDF31_nlo_as_0118";
  //string pdfName = "MMHT2014nlo68cl";
  cout << endl << "Simulation File Parameters:" << endl;
  cout << "File name: " << fileName << endl;
  cout << "Without unfolding" << endl;   // . What does this mean?
  cout << "Polarization: " << BEAM_POLARIZ << "%" << endl << endl;
  
  //calcAsym(fileName, pdfName, flagUnfolding, luminosity,beamTypeStr,Ordered_values);
  calcAsym(fileName, pdfName, beamE);   
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
	cout << GApv[i] << "   " << GcalcA[i] << endl;
      */
      double du_out, err;
      gMinuit -> GetParameter(0, du_out, err);
      Outfile << du_out << "," << err << endl;
      du_fitted[0].push_back(du_out);
      du_f_err[0].push_back(err);
      du_x[0].push_back(Gx[i]);
      /*double  Output_Sinsq[sinq2_val[0].size()];
	double Output_chi_square[chi_square[0].size()];
      */
      //std::vector<double> Output_Sinsq;
      
      delete gMinuit;
      //cout<<"SIZE"<<sinq2_val[0].size()<<endl;
      
      /*// . COMMENTED OUT THE BELOW B/C DON" 
      for (auto it = sinq2_val[lum_sets].begin();
	   it != sinq2_val[lum_sets].end(); it++) {
	//cout << *it <<endl;
	//Output_Sinsq.push_back(*it);
	Outfile_Sinsq<<lum_sets<<"," <<*it<<endl;
      }  
      
      /*for (int k=0;k<Output_Sinsq.size();k++){
	cout<<"Output_Sinsq"<<","<<Output_Sinsq[k]<<endl;
	Outfile_Sinsq<<Output_Sinsq[k]<<endl;
	}*/
      /*
      for (auto ip = chi_square[lum_sets].begin();
	   ip != chi_square[lum_sets].end(); ip++) {
	//cout<< *ip << endl;
	Outfile_chi<<lum_sets<<","<<std::setprecision(9)<<*ip<<endl;
	} */
    }
    i++;  // . Added this
    cout << "start_i=" << start_i << ", stop_i=" << stop_i;
    start_i = stop_i;  // . Added this
  }
  
  cout << endl <<endl << du_x[0].size() << endl;
  for(int i=0; i < du_x[0].size(); i++) {
    std::cout << " (du_x[0])[i]=" << (du_x[0])[i];
    std::cout << " (du_fitted[0])[i]=" << (du_fitted[0])[i];
    std::cout << " (du_f_err[0])[i]=" << (du_f_err[0])[i] << endl;
  }
  cout <<endl << endl;
  
  int length = du_x[0].size();
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
  c1->SaveAs("file_name.png");
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
