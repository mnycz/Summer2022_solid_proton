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
using namespace std;
using namespace LHAPDF;

// Create a new type that can return two values
struct Asym_Values{
  double Fit_Val;
  double Fixed_Val;
};
Double_t x = 0.0;
Double_t Q2 = 0.0;
Double_t rate = 0.0;  // . Delete this? And maybe Q2 and x?
Double_t SIN2_TH = 0.235;  

double r_prime=0.0;

const Double_t GF = 1.16637e-5;  // GeV^2, Fermi constant
const Double_t ALPHA = 7.29735e-3;  // Fine structure constant
const Double_t MP = 0.93828;  // GeV/c^2 Mass of a proton
const Double_t MN = 0.93957;
const Double_t MZ = 91.188;
const Double_t RUNTIME_DAYS = 90.0;
const Double_t BEAM_POLARIZ  =0.85;  // . edited from .8 to .85
const Int_t MAXBINS = 500;

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
Double_t Ga_rand[MAXBINS];
Double_t Gda[MAXBINS];
Double_t GdApvu_rel[MAXBINS];
Double_t GdApv_rel[MAXBINS];
Double_t Ggoodness[MAXBINS];

Double_t GApv[MAXBINS];
Double_t GcalcA[MAXBINS];

Bool_t IS_ED = false;
Int_t NBINS = 0;

//std::vector<double> sinq2_val;
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
void calcAsym(string fileName, string pdfName, double beamE)
//void calcAsym(string fileName, string pdfName, Bool_t flag_unf, Double_t lum, string beamType,Int_t Order)
 {


   cout<<"PDFNAME="<<pdfName<<","<<endl;
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
   Int_t n;
   Double_t dApv, dApvu, dA_stat;//, dA_sys, uwQ2, uwx, uwy, uwY, wY;
   Double_t dA_sys_uncor, dA_sys_cor; // . Added this but haven't yet deleted dA_sys   

   Double_t PDF_error=0.0;

   const PDF* pdf = mkPDF(pdfName, 0);
   // . change from fstream to ifstream and commented out fstream::in from the .open
   fstream inputFile;
   cout<<"fileName="<<","<<fileName<<endl;
   inputFile.open(fileName,fstream::in);
   
   TRandom *random = new TRandom();

   r_prime = random->Gaus();
   //cout<<"TESTING TO MAKE THIS THING WORK "<<inputFile<< " and "<< inputFile.is_open()<<endl<<"______________________________	\_________"<<endl;
   //string testing_file = "";
   //cout<<getline(inputFile, testing_file)<<" yup"<<endl;
   //string mystring ="";
   //inputFile >> mystring;
   //cout << mystring << " first line"<<endl;
   //double xx = 0.0;
   //double Q2Q2 = 0.0;
   //double raterate = 0.0;
   //inputFile >> xx >> Q2Q2 >> raterate;
   //cout << xx <<"," << Q2Q2 << ", "<<  raterate<< " values"<<endl;
   string file_line;
   getline(inputFile, file_line);  // remove header
   while (getline(inputFile, file_line)){  
      // The arguments in the while argument above (in number and type) should match the format in the grid file
      //istringstream iss (file_line);
      vector<string> tokens; // here we store the tokens 
      //cout<<"TESTING "<< file_line<<" !!!!!!!!!! "<<endl<<"__________________ \____________\_________"<<endl;
      string token; // current token    
      istringstream iss (file_line);
      while (getline(iss, token, ','))
       {
	 //cout << "token " << token << endl;
	 tokens.push_back(token);
	 //cout << "token " << token << endl;
	 //cout << " and tokens " << tokens.back() << endl;
       }
      Gx[i] = stof(tokens[0]);
      GQ2[i] = stof(tokens[1]);
      Grate[i] = stof(tokens[2]);
      cout<<"TESTING "<<x << ", " << Q2<<", "<< rate<<endl<<"__________________ \
____________\
_________"<<endl;
     


     /*
       Code to calculate the statistical uncertainty
       Can reuse / adapt parts of your Python code
      */

     

     /*
       Below you can store systematic uncertainties. Error_Sys,Error_1, and Error_corr are examples. They type vector (c++). Unline arrays you do not need to define their length beforehand. You technically know the size (size of the grid file) so you could translate these into arrays if you prefer ---> double Error_Sys[N]; where N is the size of grid file. Maybe useful to initialize these arrays in a for loop before using them
      */

      //Error_Sys.push_back(0);
      ////Error_Sys.push_back(0.01*GApv[i]);
      //Error_1.push_back((dApvu / BEAM_POLARIZ)*sqrt(100.0)) / lum));
      //Error_corr.push_back(0.0); // 0 for these values

    
      Gd[i] = (pdf->xfxQ2(1, Gx[i], GQ2[i]));// / Gx[i];  // d
      Gdbar[i] = (pdf->xfxQ2(-1, Gx[i], GQ2[i]));// / Gx[i];  // dbar
      Gu[i] = (pdf->xfxQ2(2, Gx[i], GQ2[i]));// / Gx[i];  // u
      Gubar[i] = (pdf->xfxQ2(-2, Gx[i], GQ2[i]));// / Gx[i];  // ubar
      Gs[i] = (pdf->xfxQ2(3, Gx[i], GQ2[i]));// / Gx[i];  // s
      Gsbar[i] = (pdf->xfxQ2(-3, Gx[i], GQ2[i]));// / Gx[i];  // sbar
      Gc[i] = (pdf->xfxQ2(4, Gx[i], GQ2[i]));// / Gx[i];  // c
      Gcbar[i] = (pdf->xfxQ2(-4, Gx[i], GQ2[i]));// / Gx[i];  // cbar
      //cout<<"TESTING SOME MORE"<<endl<<"______________________________ \_________"<<endl;
       
      double d_p = Gd[i] + Gdbar[i];
      double u_p = Gu[i] + Gubar[i];
      double s_p = Gs[i] + Gsbar[i];
      double c_p = Gc[i] + Gcbar[i];
      double d_V = Gd[i] - Gdbar[i];
      double u_V = Gu[i] - Gubar[i];
      double s_V = Gs[i] - Gsbar[i];
      double c_V = Gc[i] - Gcbar[i];
      cout << d_p << ", " << u_p << ", " << s_p << ", " <<c_p << endl;
      double C_1d = 1/2. - 2 * SIN2_TH / 3;
      double C_1u = -1/2. + 4 * SIN2_TH / 3;
      double C_2d = 1/2. - 2 * SIN2_TH;
      double C_2u = -1/2. + 2 * SIN2_TH;
      cout << SIN2_TH << ", " << C_1d << ", " << C_1u << ", " << C_2d << ", " << C_2u << endl;

      Gy[i] = GQ2[i] / (2 * MP * Gx[i] * beamE);
      double Y = (1 - (1-Gy[i]) * (1-Gy[i])) / (1 + (1-Gy[i]) * (1-Gy[i]));
      cout << Gy[i] << ", " << Y << endl;
      
      double A_const = BEAM_POLARIZ * 3 * sqrt(2) * GF * GQ2[i] / (4 * M_PI * ALPHA);  // Proportionality constant for A
      double A_tleft = 2 * (u_p + c_p) * C_1u - (d_p + s_p) * C_1d;  // Top left part of A's equation
      double A_tright = Y * (2 * (u_V + c_V) * C_2u - (d_V + s_V) * C_2d);  // Top right part of A's equation
      double A_div = 4 * (u_p + c_p) + (d_p + s_p);  // Divisor of A's equation
      cout << A_const << ", " << A_tleft << ", " << A_tright << ", " << A_div << endl;
      GApv[i] = A_const * (A_tleft + A_tright) / A_div;  // A_theory
      cout<<"GApv="<<""<<GApv[i]<<endl<<endl;

      dA_stat = 1 / sqrt(Grate[i] * RUNTIME_DAYS * 86400) / BEAM_POLARIZ;  // Statistical uncertainty in A
      Error_1.push_back(dA_stat);
      //dA_sys = 0.01*GApv[i] + ((0.01*GApv[i]*BEAM_POLARIZ) / 1);  // Systematic error
      //Gda[i]=sqrt(dA_stat*dA_stat + pow(0.01*GApv[i],2));
      dA_sys_uncor = abs(GApv[i]) * sqrt(0.002*0.002 + 0.002*0.002); // Uncorrelated systematic uncertainty: radiative corrections = 0.2%; event reconstruction = 0.2%
      Error_Sys.push_back(dA_sys_uncor);
      dA_sys_cor = abs(GApv[i]) * sqrt(0.004*0.004 + 0.002*0.002); // Correlated systematic uncertainty: polarimetry = 0.4%, Q2 = 0.2%
      Error_corr.push_back(dA_sys_cor);
      Gda[i] = sqrt(dA_stat*dA_stat + dA_sys_uncor*dA_sys_uncor);
      cout <<dA_stat << ", " << dA_sys_uncor << ", " << Gda[i] << ", " << dA_sys_cor << endl;
      Ga_rand[i] = GApv[i] + random->Gaus()*Gda[i] + r_prime*dA_sys_cor;
      //Ga_rand[i]=GApv[i] + random->Gaus()*Gda[i] //+ r_prime*(sqrt(pow(0.01*GApv[i]*1,2)));


      Double_t du_numerator = -4*Ga_rand[i] + 2*C_1u*A_const + 2*C_2u*Y*A_const;
      Double_t du_denominator = Ga_rand[i] + C_1d*A_const + C_2d*Y*A_const;
      cout << "d/u rand " << du_numerator/du_denominator << endl;

      du_numerator = -4*GApv[i] + 2*C_1u*A_const + 2*C_2u*Y*A_const;
      du_denominator = GApv[i] + C_1d*A_const + C_2d*Y*A_const;
      cout << "d/u theory " << du_numerator/du_denominator << endl;

      cout << "d/u " << Gd[i]/Gu[i] << endl;

      i++;
      counter++;
   }
   cout<<"TESTING EVEN MORE"<<endl<<"______________________________\
_________"<<endl;

   NBINS = counter;// . Changed from n + 1;
   cout << "NBINS: " << NBINS << ", " << n;
   delete pdf;
   inputFile.close();
 }
 //______________________________________________________________________________
/* Asym_Values func(Double_t x, Double_t q2, Double_t y,
     Double_t u, Double_t ubar,Double_t d, Double_t dbar,Double_t s, Double_t sbar,Double_t c, Double_t cbar,
     Double_t *par)*/
// . Made return type double instead of blank
double func(Double_t x, Double_t q2, Double_t y, Double_t u, Double_t ubar,Double_t d, Double_t dbar,Double_t s, Double_t sbar,Double_t c, Double_t cbar,Double_t *par, Int_t i)//{

 {


  Asym_Values Asym; // can access to parts of Asym
   
 
   /*
     All of these need to be change to reflect the actual pdfs need as well the what is needed for your asymmetry calculation
    */

  Double_t uplus = u+ubar;
  Double_t uv = u-ubar;
  Double_t dplus = d+dbar;
  Double_t dv = d-dbar;
  Double_t cplus = c+cbar;
  Double_t splus = s+sbar;
  
  //Double_t Y1 = 1.0+(1.0-y)*(1.0-y)-2.0*MP*MP*x*x*y*y/q2;
  //Double_t Y3 = 1.0-(1.0-y)*(1.0-y);
  //Double_t nu = (GF*q2/(sqrt(8.0)*acos(-1.0)*ALPHA)) / (1.0+q2/(MZ*MZ));

  //  Double_t gvu=0.0;
  //  Double_t gvc=0.0;
  //  Double_t gvd=0.0;
  //  Double_t gve=0.0;
  //  Double_t gvs=0.0;

  //Double_t gau = I3u; Double_t gac = I3c; Double_t gad = I3d; Double_t gas = I3s; Double_t gae = I3e;
  //Double_t gvu = I3u-2.0*Qu*par[0]; Double_t gvc = I3c-2.0*Qc*par[0]; Double_t gvd = I3d-2.0*Qd*par[0]; Double_t gvs = I3s-2.0*Qs*par[0]; Double_t gve = I3e-2.0*Qe*par[0];
    

  //Proton SFs
  //Double_t F1g_p = (Qu*Qu*uplus+Qd*Qd*dplus+Qc*Qc*cplus+Qs*Qs*splus) / 2.0;
  //Double_t F1gZ_p = Qu*gvu*uplus+Qd*gvd*dplus+Qc*gvc*cplus+Qs*gvs*splus;
  //Double_t F1Z_p = ((gvu*gvu+gau*gau)*uplus+(gvd*gvd+gad*gad)*dplus+(gvc*gvc+gac*gac)*cplus+(gvs*gvs+gas*gas)*splus) / 2.0;
  //Double_t F2g_p = 2.0*x*F1g_p;
  //Double_t F2gZ_p = 2.0*x*F1gZ_p;
  //Double_t F2Z_p = 2.0*x*F1Z_p;
  //Double_t F3gZ_p = 2.0*(Qu*gau*uv+Qd*gad*dv);
  //Double_t F3Z_p = 2.0*(gvu*gau*uv+gvd*gad*dv);
  
  //Neutron SFs
  //Double_t F1g_n = (Qd*Qd*uplus+Qu*Qu*dplus+Qc*Qc*cplus+Qs*Qs*splus) / 2.0;
  //Double_t F1gZ_n = Qd*gvd*uplus+Qu*gvu*dplus+Qc*gvc*cplus+Qs*gvs*splus;
  //Double_t F1Z_n = ((gvd*gvd+gad*gad)*uplus+(gvu*gvu+gau*gau)*dplus+(gvc*gvc+gac*gac)*cplus+(gvs*gvs+gas*gas)*splus) / 2.0;
  //Double_t F2g_n = 2.0*x*F1g_n;
  //Double_t F2gZ_n = 2.0*x*F1gZ_n;
  //Double_t F2Z_n = 2.0*x*F1Z_n;
  //Double_t F3gZ_n = 2.0*(Qd*gad*uv+Qu*gau*dv);
  //Double_t F3Z_n = 2.0*(gvd*gad*uv+gvu*gau*dv);

  double A_const = BEAM_POLARIZ * 3 * sqrt(2) * GF * GQ2[i] / (4 * M_PI * ALPHA);  // Proportionality constant for A
  double Y = (1 - (1-Gy[i]) *  (1-Gy[i])) / (1 +  (1-Gy[i]) *  (1-Gy[i]));
  
  
  Double_t A_numerator;
  Double_t A_denominator;
  double C_1d = 1/2. - 2 * SIN2_TH / 3;
  double C_1u = -1/2. + 4 * SIN2_TH / 3;
  double C_2d = 1/2. - 2 * SIN2_TH;
  double C_2u = -1/2. + 2 * SIN2_TH;



  /*
    Can use A_numerator and A_denominator to seperate the Asym
   */
  A_numerator = A_const * (2.0*C_1u + 2.0*Y*C_2u - par[0]*(C_1d + Y*C_2d)); //nu*(2.0*gae*Y1*F1gZ_p+gve*Y3*F3gZ_p);
  A_denominator = par[0] + 4;//2.0*Y1*F1g_p-nu*(2.0*gve*Y1*F1gZ_p+gae*Y3*F3gZ_p);
  Double_t du_numerator = -4*Ga_rand[i] + 2*C_1u*A_const + 2*C_2u*Y*A_const;
  Double_t du_denominator = Ga_rand[i] + C_1d*A_const + C_2d*Y*A_const;
  /* Fix above equations to the new expression
   */
  //cout << A_numerator << ", and " << A_denominator << endl;
  cout << "d/u " << du_numerator/du_denominator << endl;

    return (A_numerator/A_denominator); //return type needs to be of the form of double
  //return Asym;
 }
 //______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
 {
   
   Asym_Values Calculation;
   TMatrixT<double> Err_matrix(counter,counter);
   TMatrixT<double> Stat_Err_matrix(counter,counter);
   //cout<<"counter size"<<counter<<endl;
   for (int i=0;i<counter;i++){
     for (int j=0;j<counter;j++){
       if (i==j){
	 //Stat_Err_matrix(i,j)=pow(Error_1[i],2);
	 Stat_Err_matrix(i,j)= (pow(Error_1[i],2)  + pow(Error_Sys[i],2) + pow(Error_corr[i],2));
	 //Stat_Err_matrix(i,j)= (pow(Error_1[i],2)  + pow(Error_Sys[i],2));      
       }
       else{	 	 	
	 Stat_Err_matrix(i,j) = 0.0;  // . might be dA_sys_cor^2
	   /* This assumes that there are no off-diagonal elements
	    */
       }
     }
   }

   //Stat_Err_matrix.Print();
   //cout<<"Check Stat Matrix"<<","<<Stat_Err_matrix(0,0)<<","<<Stat_Err_matrix(1,1) <<endl;

   TMatrixT<double> Diff(1,counter);
   TMatrixT<double> Diff_T(counter,1);
   Err_matrix = Stat_Err_matrix;

 //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    Double_t Pseduo=0.0;
    TRandom *random_2 = new TRandom();
    for(int i = 0; i < NBINS; ++i){
      GcalcA[i] = func(Gx[i],GQ2[i],Gy[i],Gu[i],Gubar[i],Gd[i],Gdbar[i],Gs[i],Gsbar[i],Gc[i],Gcbar[i],par, i);  // . added the i
      //Calculation=func(Gx[i],GQ2[i],Gy[i],Gu[i],Gubar[i],Gd[i],Gdbar[i],Gs[i],Gsbar[i],Gc[i],Gcbar[i],par);  
      //GcalcA[i]=Calculation.Fit_Val;
      //Pseduo= Calculation.Fixed_Val + random_2->Gaus()*Gda[i] + r_prime*(sqrt(pow(0.01*Calculation.Fit_Val*1,2)));
      
      //cout << Ga_rand[i] << ", " << GcalcA[i] << ", " << Gda[i] << "A terms" <<endl;
      delta  = (Ga_rand[i]-GcalcA[i]);///Gda[i];
      //delta = (Pseduo-GcalcA[i])/Gda[i]; 
      //chisq += (abs(Ggoodness[i] - 1.0000) < 0.1)*delta*delta;
      //Diff(0,i) = (abs(Ggoodness[i] - 1.0000) < 0.1)*pow((Ga_rand[i]-GcalcA[i]),1);
      //Diff(0,i) = (abs(Ggoodness[i] - 1.0000) < 0.1)*pow((Pseduo-GcalcA[i]),1); 
      // Diff(0,i) =pow((Pseduo-GcalcA[i]),1);
      Diff(0,i) = delta;//*delta;
      //cout << "delta " << delta;
      //cout << " Diff" << Diff(0,i) << endl;
    }
    //Diff.Print();
    TMatrixT<double> Invt(counter,counter);
    
    Err_matrix.SetTol(1.e-26);
    double det = Err_matrix.Determinant(); 
    cout << det << "<- det" << endl;
    Double_t det_2=0.0;
    Invt = Err_matrix.Invert(&det_2);
    //Invt = Err_matrix.InvertFast(&det_2);
    //cout<<"Valid"<<","<<Invt.IsValid()<<","<<det<<","<<det_2<<endl;
    //cout << "Invt:";
    //Invt.Print();
    Diff_T.Transpose(Diff);
    TMatrixT<double> Temp(counter,counter);
    TMatrixT<double> chi(1,counter);
    //Temp.Mult(Diff,Err_matrix);
    Temp.Mult(Diff,Invt);  
    chi.Mult(Temp,Diff_T);
    ///////////cout << "temp" << Temp << endl << "chi" << chi << endl;
   //TMatrix chi = Temp.Mult
    //chisq += (Diff
    //cout << "temp;";
    //Temp.Print();
    //cout << "chi:";
    //chi.Print();
    // .chi_square[lum_sets].push_back(chi.Sum());
    //cout<<"chisq"<<","<<chisq<<endl;
    //}
    //f = chisq;
    f=chi.Sum();
 }
 //______________________________________________________________________________
 
int main(Int_t argc, char *argv[])
{ // Added the 1 \/
    if(argc < 1){
      cout << "Please call with appropriate arguments:" << endl << "./codeName beamType luminosity withUnfolding(true/false) inputFile" << endl;
      return -1;
    }
    
    string fileName = argv[1];
    
    string pdfName = "CT18NLO";
    double beamE = 11; // GeV, beamline energy
    //string pdfName = "NNPDF31_nlo_as_0118";
    //string pdfName = "MMHT2014nlo68cl";

    cout << endl << "Simulation File Parameters:" << endl;
    cout << "File name: " << fileName << endl;
    //cout << "Beam luminosity: " << luminosity << endl;
    cout << "Without unfolding" << endl;
    cout << "Polarization: 80%" << endl << endl;
    
    //calcAsym(fileName, pdfName, flagUnfolding, luminosity,beamTypeStr,Ordered_values);
    calcAsym(fileName, pdfName, beamE);   
    //return 0; // . ADDED THIS TO STOP FROM GOING FURHTER AND BREAKING THINGS
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Instantiate Minuit for 1 parameter
    TMinuit *gMinuit = new TMinuit(1);//4); 
    gMinuit->SetFCN(fcn);  // Set the address of the minimization function
    
    //Perform fit with betas fixed
    //--------------------------------------------------------------------
    cout << endl << "Betas fixed" << endl;
    Double_t arglist[10];
    Int_t ierflg = 0;
 
    // . I commented out the next two lines because https://llr.in2p3.fr/activites/physique/glast/workbook/pages/advanced_GRUG/GRUGminuit.pdf  doesn't have them and I can't find a good explanation of what they do
    arglist[0] = 1;  // 1 because chisquare fit
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
 
    static Double_t vstart[1] = {0.0};//0.1,0,0,0}; here and right below was 4, not 1, and here had 0.1, not 0.0
    static Double_t step[1] = {0.001};//,0.001,0.001,0.001};
    gMinuit->mnparm(0, "d/u", vstart[0], step[0], 0, 1, ierflg); // or 0,0
    //gMinuit->mnparm(0, "sin2th", vstart[0], step[0], 0,0,ierflg);
    ////gMinuit->mnparm(1, "Polarization", vstart[1], step[1], 0,0,ierflg);
    //gMinuit->mnparm(1, "beta_ht", vstart[2], step[2], 0,0,ierflg);
    ////gMinuit->mnparm(2, "Polarization", vstart[3], step[3], 0,1,ierflg);
    //gMinuit->mnparm(2, "beta_csv", vstart[3], step[3], 0,0,ierflg);
    ////gMinuit->mnparm(3, "Polarization", vstart[3], step[3], 0,1,ierflg);
    
    // . I Commented out the next four lines because they seem like they don't need to be fixed, but I only kind of understand what is going on. I uncommented the 500 lines since it's in https://llr.in2p3.fr/activites/physique/glast/workbook/pages/advanced_GRUG/GRUGminuit.pdf
    arglist[0] = 500; // Max calls
    arglist[1] = .001;  // Tolerance  (was 1)
    //gMinuit->FixParameter(1);
    //gMinuit->FixParameter(2);
    
    //gMinuit->FixParameter(3);
    //gMinuit->SetPrintLevel(1);
    double tmp_sin=0.0;
    double tmp_error=0.0;
    //gMinuit -> GetParameter(0,tmp_sin,tmp_error);
    //cout<<"tmp_sin"<<","<<tmp_sin<<","<<tmp_error<<endl;
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
    double sin2_out,err;
    gMinuit -> GetParameter(0,sin2_out,err);
    Outfile<<sin2_out<<","<<err<<endl;

    /*double  Output_Sinsq[sinq2_val[0].size()];
    double Output_chi_square[chi_square[0].size()];
    */
    //std::vector<double> Output_Sinsq;

   delete gMinuit;
   //cout<<"SIZE"<<sinq2_val[0].size()<<endl;
   
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


   for (auto ip = chi_square[lum_sets].begin();
        ip != chi_square[lum_sets].end(); ip++) {
     //cout<< *ip << endl;
     Outfile_chi<<lum_sets<<","<<std::setprecision(9)<<*ip<<endl;
   }

   

   return 0;
 }
