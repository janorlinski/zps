#include<vector>
#include<fstream>


int GroupProject () {




std::vector <double> BIN ;
std::vector <double> LICZ ; 	
	
	
//wczytywanie z pliku do wektorów	

std::ifstream myStr ; 

double binx, countsx;
myStr.open("20umRun47no0125Ge12.txt");

if(!myStr) { // file couldn't be opened
      std::cerr << "Error: file could not be opened" <<std::endl;
      exit(1);
   }



double countsy=0;

while (myStr>>countsx)
{

	
	countsy=countsx ;
	
	LICZ.push_back(countsy);


  }
  myStr.close (); 
  
   std::cout<<LICZ.size()<<std::endl;
   
   double dim = LICZ.size() ;
   
//wpisywanie z wektorów do histogramu .
  
TH1I* Spec_uncal= new TH1I ("Spec_uncal","myhist",dim,0.,dim-1.);

  for (int i=0; i<LICZ.size(); i++) {
  
    Spec_uncal->Fill(i,LICZ[i]);
    
    
    }
    
    //rysowanie nieskalibrowanego widma
    Spec_uncal->Draw("hist");
 	

 Int_t i, NFound, binnum;
  Int_t NBins = Spec_uncal->GetNbinsX ();
  Double_t* source = new Double_t [NBins];
  Double_t* destin = new Double_t [NBins];

  for (i = 200; i < NBins-200; i++)  source[i] = Spec_uncal->GetBinContent (i+1);

  TSpectrum* Analyser1dim = new TSpectrum (2); // note: default maxpositions = 100

  // Searching for Candidate Peaks' positions

  NFound = Analyser1dim->SearchHighRes (source, destin, NBins, 
                                        10., 50.0, kFALSE, 100, kFALSE, 0);




/*
  Int_t TSpectrum::SearchHighRes (
     Double_t* source,           : pointer to the vector of source spectrum
     Double_t* destVector,       : pointer to the vector of resulting deconvolved spectrum
     Int_t     ssize,            : length of source spectrum
     Double_t  sigma,            : sigma of searched peaks (see manual)
     Double_t  threshold,        : threshold value in % for selected peaks. Peaks with
                                   amplitude < threshold*highest_peak/100 are ignored, see manual
     bool      backgroundRemove, : logical variable, set if the removal of background 
                                   before deconvolution is desired
     Int_t     deconIterations,  : number of iterations in deconvolution operation
     bool      markov,           : true = first the source spectrum is replaced by new spectrum
                                   calculated using Markov chains method
     Int_t     averWindow        : (Markov: ) averaging window of searched peaks, (see manual)
  )

  The goal of this function is to identify automatically the peaks in spectrum with the
  presence of the continuous background and statistical fluctuations - noise.
  
  1. Background is removed (if desired), 
  2. Markov smoothed spectrum is calculated (if desired), 
  3. Response function is generated according to given sigma + Deconvolution is carried out.

  See root.cern.ch/doc/v608/classTSpectrum.html#a5ba181a45b117934b48c4ef5f78d0b2b
*/


  cout << "\n* SearchHighRes has found: " << NFound << " peaks.\n\n";


Double_t* Posit  = Analyser1dim->GetPositionX ();

Double_t* Amplit = new Double_t [NFound];
  

  Bool_t* FixPosit  = new Bool_t [NFound];          // should these parameters be fixed? [T/F]
  Bool_t* FixAmplit = new Bool_t [NFound];

  for (i = 0; i < NFound; i++) 
  {
    FixAmplit[i] = kFALSE;
    FixPosit [i] = kFALSE;

    binnum = 1 + Int_t (Posit[i] + 0.5);            // the "nearest" bin
    Amplit   [i] = Spec_uncal->GetBinContent (binnum);   // estimating the peak amplitude
  }

  cout << "* Press Enter to fit the spectra using the AWMI approach.";
  cin.ignore();

TSpectrumFit* Fitter1dim = new TSpectrumFit (NFound);
 Fitter1dim->SetFitParameters (200, NBins-1 , 1000, 0.1, Fitter1dim->kFitOptimChiCounts,
                                Fitter1dim->kFitAlphaHalving, Fitter1dim->kFitPower2,
                                Fitter1dim->kFitTaylorOrderFirst );
                                

  Fitter1dim->SetPeakParameters (2., kFALSE, Posit, FixPosit, Amplit, FixAmplit);                                

/* 
  void TSpectrumFit::SetPeakParameters	(
      Double_t        sigma,        : initial value of sigma parameter
      Bool_t          fixSigma,     : true if sigma should be fixed
      const Double_t* positionInit, : array of initial values of peaks positions
      const Bool_t*   fixPosition,  : array for peaks. True if a given position should be fixed
      const Double_t* ampInit,      : array of initial values of peaks amplitudes
      const Bool_t*   fixAmp        : array for amplitudes. True if a given amp should be fixed
  )
  See: 
*/



 Fitter1dim->FitAwmi (source);
 
 
  Double_t* Positions        = Fitter1dim-> GetPositions ();       // Caution: Xscale is Bin No.
  Double_t* PositionsErrors  = Fitter1dim-> GetPositionsErrors (); // Caution: Xscale is Bin No.
  Double_t* Amplitudes       = Fitter1dim-> GetAmplitudes ();
  Double_t* AmplitudesErrors = Fitter1dim-> GetAmplitudesErrors ();
  Double_t* Areas            = Fitter1dim-> GetAreas ();           // Caution: Xscale is Bin No.
  Double_t* AreasErrors      = Fitter1dim-> GetAreasErrors ();     // Caution: Xscale is Bin No.
  
   Double_t x1 = Spec_uncal->GetBinCenter(1), dx = Spec_uncal->GetBinWidth (1);

  Double_t sigma, sigmaErr;
  Fitter1dim->GetSigma (sigma, sigmaErr);                          // Caution: Xscale is Bin No.

  // Current TSpectrumFit needs a sqrt(2) correction factor for sigma

  sigma /= TMath::Sqrt2(); sigmaErr /= TMath::Sqrt2();

  // Convert "bin numbers" into "x-axis values"

  sigma *= dx; sigmaErr *= dx;

  cout << "\n * Overall sigma found by fit: " << sigma << " (+-" << sigmaErr << ")" << endl;
  cout << " * Fit chi^2 = " << Fitter1dim->GetChi () << "\n\n";
 
 
   for (i = 0; i < NFound; i++)
  {
    // Convert "bin numbers" into "x-axis values"
    Positions      [i] = x1 + Positions[i] * dx;
    PositionsErrors[i] *= dx;
    Areas          [i] *= dx;
    AreasErrors    [i] *= dx;

    cout << "Found: "
         << Positions [i] << " (+-" << PositionsErrors [i] << ") "
         << Amplitudes[i] << " (+-" << AmplitudesErrors[i] << ") "
         << Areas     [i] << " (+-" << AreasErrors     [i] << ")" << endl;
  }
 
 
 
 
 /* 
  
Double_t Chan_vector[3] = {138.666,261.389,360.817} ; //this table holds the initial peak positions in a uncallibrated spectrum
Double_t En_vector[3] = {121.7824,244.6989,344.2811} ;
Double_t Chan_vector_errors[3] = {1.01,1.01,1.01} ;
Double_t En_vector_errors[3] = {0.0004,0.001,0.009} ;




	
//initial fit to the provided peak energies


	TF1 *linear = new TF1("f3", "[0]+[1]*x",0,4095);
	linear->SetParameters(1,1);
auto gr = new TGraphErrors(3,Chan_vector,En_vector,Chan_vector_errors,En_vector_errors);
        gr->Fit(linear);
        
TCanvas *c2 = new TCanvas ("c2");	
 
   gr->SetTitle("TGraphErrors initial");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("AP");

*/



  return 0;
}
