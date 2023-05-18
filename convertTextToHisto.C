// callibrator.C

#include <iostream>
#include <fstream>
#include <vector>

#include "TH1I.h"
#include "TFile.h"
#include "TSpectrum.h"
#include "TSpectrumFit.h"

TH1I getHistoFromTxtFile (TString fileName) {
	
	//setup
	
	Bool_t verbatim = true; //choose whether you want to see all debug prompts
	std::vector <Int_t> value;	
	
	//wczytywanie z pliku do wektorów	

	std::ifstream myStr;
	Int_t value_temp;
	myStr.open(fileName);

	if(!myStr) { // file couldn't be opened
		std::cerr << "Error: file could not be opened" <<std::endl;
		exit(1);
	}
	
	string dummy;
	getline(myStr, dummy);
	//cout << dummy << endl;
	
	while (myStr >> value_temp) value.push_back(value_temp);
	myStr.close (); 
   
	//wpisywanie z wektorów do histogramu 
	
	Int_t nChannels = value.size();
	if (verbatim) std::cout << "Loaded file " << fileName << "; nChannels = " << nChannels << ". \n"; 
	
	TH1I Spec_uncal("Spec_uncal", "myhist", nChannels, 0.5, nChannels+0.5);
	for (int i=0; i<nChannels; i++) Spec_uncal.SetBinContent(i+1, value[i]);
    
	//editing the uncalibrated spectrum 
	
	Spec_uncal.SetTitle("Uncalibrated spectrum '" + fileName + "'");
	
	return Spec_uncal;
	
}

vector<TH1I> getHistosFromTxtFileList (TString fileList) {
	
	//setup
	
	Bool_t verbatim = false; //choose whether you want to see all debug prompts
	TString fileNameRoot = "output.root";
	if (verbatim) std::cout << "hello" << std::endl;
	std::vector < TH1I > histos;
	//wczytywanie z pliku do wektorów	

	std::ifstream myStr;
	TString file;
	myStr.open(fileList);

	if(!myStr) { // file couldn't be opened
		std::cerr << "Error: file could not be opened" <<std::endl;
		exit(1);
	}
	
	TFile* fout = new TFile(fileNameRoot, "recreate");
	
	while (myStr >> file) {
		
		TString fileName = TString(file);
		
		if (verbatim) cout << "Loading txt file '" << fileName << "'. \n"; 
		TH1I histo = getHistoFromTxtFile(fileName);
		histo.SetName(fileName);
		histo.Write();
		histos.push_back(histo);
	}
	
	myStr.close (); 
	fout->Close();
	return histos;
}

void fitTwoGaussPeaksToHisto (TH1I Histogram) {
	
	//std::cout << "hello this is fitTwoGauss..." << std::endl;
	
	int i, NFound, NBins, binnum;
	NBins = Histogram.GetNbinsX();
	int loCut = 200;
	int hiCut = 30000;
	//NBinsTrue = NBins - (loCut + hiCut);
	
	Double_t* source = new Double_t [NBins];
	Double_t* destin = new Double_t [NBins];
	
	for (i = 0; i < NBins; i++) {
		source[i] = Histogram.GetBinContent (i+1);
		if (i<loCut || i>hiCut) source[i] = 0;
		//cout << i << ": " << source[i] << endl;
		
	}
	
	TSpectrum* Analyser1dim = new TSpectrum (2); 
	NFound = Analyser1dim->SearchHighRes (source, destin, NBins, 10.0, 50.0, kFALSE, 100, kFALSE, 0);

	cout << "SearchHighRes has found: " << NFound << " peaks \n";

	Bool_t* FixPosit  = new Bool_t [NFound];          // should these parameters be fixed? [T/F]
	Bool_t* FixAmplit = new Bool_t [NFound];

	Double_t* Posit  = Analyser1dim->GetPositionX ();
	Double_t* Amplit = new Double_t [NFound];

	for (i = 0; i < NFound; i++) {
		
		FixAmplit[i] = kFALSE;
		FixPosit [i] = kFALSE;
		binnum = int(Posit[i]);            // the "nearest" bin
		Amplit[i] = Histogram.GetBinContent(binnum);   // estimating the peak amplitude
		
	}

	TSpectrumFit* Fitter1dim = new TSpectrumFit (NFound);
	Fitter1dim->SetFitParameters (loCut, hiCut, 1000, 0.1, Fitter1dim->kFitOptimChiCounts,
                                Fitter1dim->kFitAlphaHalving, Fitter1dim->kFitPower2,
                                Fitter1dim->kFitTaylorOrderFirst );
                                
	Fitter1dim->SetPeakParameters (10., kFALSE, Posit, FixPosit, Amplit, FixAmplit);                                
	Fitter1dim->FitAwmi (source);
		
	Double_t* Positions        = Fitter1dim-> GetPositions ();       // Caution: Xscale is Bin No.
	//Double_t* PositionsErrors  = Fitter1dim-> GetPositionsErrors (); // Caution: Xscale is Bin No.
	//Double_t* Amplitudes       = Fitter1dim-> GetAmplitudes ();
	//Double_t* AmplitudesErrors = Fitter1dim-> GetAmplitudesErrors ();
	//Double_t* Areas            = Fitter1dim-> GetAreas ();           // Caution: Xscale is Bin No.
	//Double_t* AreasErrors      = Fitter1dim-> GetAreasErrors ();     // Caution: Xscale is Bin No.
  
	//Double_t x1 = Histogram.GetBinCenter(1);
	//Double_t dx = Histogram.GetBinWidth (1);
	
	//cout << "\n * Overall sigma found by fit: " << sigma << " (+-" << sigmaErr << ")" << endl;
	//cout << " * Fit chi^2 = " << Fitter1dim->GetChi () << "\n\n";
 
 
   for (i = 0; i < NFound; i++)
  {
    // Convert "bin numbers" into "x-axis values"
    //Positions      [i] = x1 + Positions[i] * dx;
    //PositionsErrors[i] *= dx;
    //Areas          [i] *= dx;
    //AreasErrors    [i] *= dx;

    //cout << "Found: "
         //<< Positions [i] << " (+-" << PositionsErrors [i] << ") "
         //<< Amplitudes[i] << " (+-" << AmplitudesErrors[i] << ") "
         //<< Areas     [i] << " (+-" << AreasErrors     [i] << ")" << endl;
         
		std::cout << "Found position: " << Positions[i] << std::endl;
  }
  
  
  
	
}

int convertTextToHisto () {

	TString listName = "lista_jedenplik.txt";
	std::cout << "Selected listName '" << listName << "'\n";
	
	std::vector<TH1I> histos = getHistosFromTxtFileList(listName);
	Int_t nHistos = histos.size();
	//std::cout << histos[0].GetNbinsX() << std::endl;

	for (Int_t i=0; i<nHistos; i++) {
		fitTwoGaussPeaksToHisto(histos[i]);
	}
	
	return 0;
	
}
