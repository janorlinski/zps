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

Double_t gausswithlinearbkg (Double_t *xarg, Double_t *par)
{
 Double_t x = xarg[0] , result = 0.;
 result = par[0]/(par[2]*TMath::Sqrt(2*TMath::Pi())) * TMath::Exp((-TMath::Power((x-par[1]),2))/(2*TMath::Power(par[2],2))) + par[3]*x + par[4];
 return result;
}


void fitTwoGaussPeaksToHisto (TH1I Histogram, TFile* file) {
	
	std::cout << "fitting two Gauss peaks to file " << Histogram.GetName() << std::endl;
	
	int i, NFound, NBins, binnum;
	NBins = Histogram.GetNbinsX();
	int loCut = 2000;
	int hiCut = 5500;
	double fitrange = 50.0;
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

	if (NFound < 2 || NFound > 2) {
		std::cout << "N peaks not eq to 2, abort fitting \n";
		return;
	}

	Double_t* Posit  = Analyser1dim->GetPositionX ();
	//Double_t* Amplit = new Double_t [NFound];

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

	std::cout << "Found position: " << dopasowanie1->GetParameter(1) << " +/- " << dopasowanie1->GetParError(1) << "\n";
	std::cout << "Found position: " << dopasowanie2->GetParameter(1) << " +/- " << dopasowanie2->GetParError(1) << "\n";
	
	dopasowanie1->Draw();
	dopasowanie2->Draw();
	
	file->cd();
	
	Histogram.Write();
	
	//dopasowanie1->Write();
	//dopasowanie2->Write();

}


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
		fitTwoGaussPeaksToHisto(histo, fout);
		//histo.Write();
		histos.push_back(histo);
	}
	
	myStr.close (); 
	fout->Close();
	return histos;
}

int convertTextToHisto () {

	TStopwatch t;
	t.Start();

	TString listName = "lista_jedenplik.txt";
	std::cout << "Selected listName '" << listName << "'\n";
	
	std::vector<TH1I> histos = getHistosFromTxtFileList(listName);
	Int_t nHistos = histos.size();
	//std::cout << histos[0].GetNbinsX() << std::endl;

	//for (Int_t i=0; i<nHistos; i++) {
		//fitTwoGaussPeaksToHisto(histos[i]);
	//}
	
	Double_t time = t.RealTime();
	std::cout << "Done analyzing " << nHistos << " spectra in " << time << " seconds. Av = " << time/nHistos << " seconds/spec. \n";
	return 0;

}
