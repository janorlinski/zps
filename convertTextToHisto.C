// callibrator.C

#include <iostream>
#include <fstream>
#include <vector>

#include "TH1I.h"
#include "TFile.h"
#include "TSpectrum.h"
#include "TSpectrumFit.h"
#include "TF1.h"
#include "TStopwatch.h"

Double_t gausswithlinearbkg (Double_t *xarg, Double_t *par) {
	
	Double_t x = xarg[0];
	Double_t result = 0.0;
 
	result = par[0]/(par[2]*TMath::Sqrt(2*TMath::Pi())) * TMath::Exp((-TMath::Power((x-par[1]),2))/(2*TMath::Power(par[2],2))) + par[3]*x + par[4];
	return result;
	
}


vector <double> fitTwoGaussPeaksToHisto (TH1I Histogram, TFile* file) {

	// 279.01 keV
	// 547.5  keV
	
	std::cout << "fitting two Gauss peaks to file " << Histogram.GetName() << std::endl;
	
	int i, NFound, NBins, binnum;
	NBins = Histogram.GetNbinsX();
	int loCut = 2000;
	int hiCut = 5500;
	double fitrange = 50.0;
	std::vector < double > controlValues(4);
	controlValues[0] = 0.0;
	controlValues[1] = 0.0;
	controlValues[2] = 0.0;
	controlValues[3] = 0.0;
	
	Double_t* source = new Double_t [NBins];
	Double_t* destin = new Double_t [NBins];
	
	for (i = 0; i < NBins; i++) {
		source[i] = Histogram.GetBinContent (i+1);
		if (i<loCut || i>hiCut) source[i] = 0;
	}
	
	TSpectrum* Analyser1dim = new TSpectrum (2); 
	NFound = Analyser1dim->SearchHighRes (source, destin, NBins, 10.0, 50.0, kFALSE, 100, kFALSE, 0);

	cout << "SearchHighRes has found: " << NFound << " peaks \n";

	file->cd();

	if (NFound < 2 || NFound > 2) {
		std::cout << "N peaks not eq to 2, abort fitting \n";
		Histogram.Write();
		return controlValues;
	}

	Double_t* Posit  = Analyser1dim->GetPositionX ();

	TF1* dopasowanie1 = new TF1("dopasowanie1", gausswithlinearbkg, Posit[0]-fitrange, Posit[0]+fitrange, 5);
	TF1* dopasowanie2 = new TF1("dopasowanie2", gausswithlinearbkg, Posit[1]-fitrange, Posit[1]+fitrange, 5);
           
	dopasowanie1->SetParameters(1000, Posit[0], 10, 0.0000000001, 0.00001);
	dopasowanie1->SetParNames("Amplituda", "Srednia", "Sigma", "A", "B");

	dopasowanie2->SetParameters(1000, Posit[1], 10, 0.0000000001, 0.00001);
	dopasowanie2->SetParNames("Amplituda", "Srednia", "Sigma", "A", "B");
	  
	Histogram.Fit("dopasowanie1","MIQR+","",Posit[0]-fitrange, Posit[0]+fitrange); //MIQR
	Histogram.Fit("dopasowanie2","MIQR+","",Posit[1]-fitrange,Posit[1]+fitrange); //MIQR
	 
	Double_t sigma1 = dopasowanie1->GetParameter(2);
	Double_t sigma2 = dopasowanie2->GetParameter(2);
	 
	Double_t chi21 = dopasowanie1->GetChisquare();
	Double_t chi22 = dopasowanie2->GetChisquare();
	 
	Double_t NDF1 = dopasowanie1->GetNDF();
	Double_t NDF2 = dopasowanie2->GetNDF();
	 
	Double_t chi2NDF1 = chi21 / NDF1;
	Double_t chi2NDF2 = chi22 / NDF2;
	
	controlValues[0] = sigma1;
	controlValues[1] = sigma2;
	controlValues[2] = chi2NDF1;
	controlValues[3] = chi2NDF1;

	std::cout << "Found position: " << dopasowanie1->GetParameter(1) << " +/- " << dopasowanie1->GetParError(1) << "\n";
	std::cout << "Found position: " << dopasowanie2->GetParameter(1) << " +/- " << dopasowanie2->GetParError(1) << "\n";
	
	Histogram.Write();
	
	// obliczanie kalibracji
	
	Double_t y1 = 279.01;
	Double_t y2 = 547.5;
	
	Double_t x1 = dopasowanie1->GetParameter(1);
	Double_t x2 = dopasowanie2->GetParameter(1);
	
	Double_t calibFunctionSlope = (y1 - y2)/(x1 - x2);
	Double_t calibFunctionYIntercept  = (-1.0) * calibFunctionSlope * x1 + y1;
	
	TString name = Histogram.GetName();
    name.ReplaceAll(".txt", ".cal");
    
	ofstream calibOut(name);
	calibOut << calibFunctionYIntercept << "\n";
	calibOut << calibFunctionSlope << "\n";
    calibOut.close();

	return controlValues;

}


TH1I getHistoFromTxtFile (TString fileName) {
	
	//setup
	
	Bool_t verbatim = false; //choose whether you want to see all debug prompts
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
	TString fileNameRoot = fileList;
	fileNameRoot.ReplaceAll(".list", ".root");
	if (verbatim) std::cout << "hello" << std::endl;
	
	std::vector < TH1I > histos;
	std::vector < TString > names;
	std::vector < double > sigma1vec;
	std::vector < double > sigma2vec;
	std::vector < double > chi2NDF1vec;
	std::vector < double > chi2NDF2vec;
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
		vector<Double_t> controlValues = fitTwoGaussPeaksToHisto(histo, fout);
		
		histos.push_back(histo);
		names.push_back(fileName);
		sigma1vec.push_back(controlValues[0]);
		sigma2vec.push_back(controlValues[1]);
		chi2NDF1vec.push_back(controlValues[2]);
		chi2NDF2vec.push_back(controlValues[3]);
	}
	
	Double_t NHistos = histos.size();
	
	TH1F* sigma1 = new TH1F("sigma1", "sigma1", NHistos, 0.5, NHistos+0.5);
	TH1F* sigma2 = new TH1F("sigma2", "sigma2", NHistos, 0.5, NHistos+0.5);
	TH1F* chi2NDF1 = new TH1F("chi2NDF1", "chi2NDF1", NHistos, 0.5, NHistos+0.5);
	TH1F* chi2NDF2 = new TH1F("chi2NDF2", "chi2NDF2", NHistos, 0.5, NHistos+0.5);
	
	for (int i=0; i<NHistos; i++) {
		
		sigma1->SetBinContent(i+1, sigma1vec[i]);
		sigma2->SetBinContent(i+1, sigma2vec[i]);
		chi2NDF1->SetBinContent(i+1, chi2NDF1vec[i]);
		chi2NDF2->SetBinContent(i+1, chi2NDF2vec[i]);
		
		sigma1->GetXaxis()->SetBinLabel(i+1, names[i]);
		sigma2->GetXaxis()->SetBinLabel(i+1, names[i]);
		chi2NDF1->GetXaxis()->SetBinLabel(i+1, names[i]);
		chi2NDF2->GetXaxis()->SetBinLabel(i+1, names[i]);
		
	}
	
	myStr.close (); 
	sigma1->Write();
	sigma2->Write();
	chi2NDF1->Write();
	chi2NDF2->Write();
	fout->Close();
	return histos;
	
}

int convertTextToHisto () {

	TStopwatch t;
	t.Start();

	TString listName = "lista_trzypliki.list";
	std::cout << "Selected listName '" << listName << "'\n";
	
	std::vector<TH1I> histos = getHistosFromTxtFileList(listName);
	Int_t nHistos = histos.size();
	
	Double_t time = t.RealTime();
	std::cout << "Done analyzing " << nHistos << " spectra in " << time << " seconds. Av = " << time/nHistos << " seconds/spec. \n";
	return 0;

}
