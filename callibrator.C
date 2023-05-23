// callibrator.C

#include <iostream>
#include <fstream>
#include <vector>

#include "TH1I.h"
#include "TH1F.h"
#include "TFile.h"
#include "TSpectrum.h"

Double_t gausswithlinearbkg (Double_t *xarg, Double_t *par)
{
 Double_t x = xarg[0] , result = 0.;
 result = par[0]/(par[2]*TMath::Sqrt(2*TMath::Pi())) * TMath::Exp((-TMath::Power((x-par[1]),2))/(2*TMath::Power(par[2],2))) + par[3]*x + par[4];
 return result;
}

Double_t mgauss(Double_t *xarg, Double_t *par)
{
 Double_t x = xarg[0] , result = 0.;
 result = par[0]/(par[2]*TMath::Sqrt(2*TMath::Pi())) * TMath::Exp((-0.5*TMath::Power((x-par[1]),2))/(TMath::Power(par[2],2)));
 return result;
}
   

void dop (TH1I* hspect){
    
     TSpectrum* Analyser1dim = new TSpectrum (20);
     
     Int_t NFound = Analyser1dim->Search (hspect, 10, "", 0.5);
     
    Int_t NPeaks = 0;
     
     Double_t* xpeaks = Analyser1dim->GetPositionX ();
     Double_t* ypeaks = Analyser1dim->GetPositionY ();
     
     Double_t fitrange = 50;
     
    Double_t wyniki[2];
     
     
     for (Int_t p = 0; p < NFound; p++)
      {
        Double_t xp = xpeaks[p];
        Int_t   bin = hspect->GetXaxis()->FindBin (xp);
        Double_t yp = hspect->GetBinContent (bin);
          
          
       if ((xp-2*fitrange)<0) { continue;}
          
     TF1* dopasowanie1 = new TF1("dopasowanie1",gausswithlinearbkg,xp-fitrange,xp+fitrange,5);
          
        //  TF1* dopasowanie = new TF1("dopasowanie",mgauss,xp-fitrange,xp+fitrange,3);
          
          dopasowanie1->SetParameters(1000, xp, 50, 0.0000000001, 0.00001);
          dopasowanie1->SetParNames("Amplituda", "Srednia", "Sigma", "A", "B");
          
      hspect->Fit("dopasowanie1","MIQR","",xp-fitrange,xp+fitrange); //MIQR
          
    Double_t sigma = dopasowanie1->GetParameter(2);
          
    TF1* dopasowanie = new TF1("dopasowanie",gausswithlinearbkg,xp-5*sigma,xp+5*sigma,5);
          
          dopasowanie->SetParameters(1000, xp, sigma, 0.0000000001, 0.00001);
          dopasowanie->SetParNames("Amplituda", "Srednia", "Sigma", "A", "B");
          
           hspect->Fit("dopasowanie","MIQR","",xp-4*sigma,xp+4*sigma); //MIQR
          
          Double_t srodek = dopasowanie->GetParameter(1);
          Double_t err = dopasowanie->GetParError(1);
          
          
      hspect->Draw("");
        
          
          if(NPeaks<2){
          wyniki[NPeaks]=srodek;
          }
          
          
        NPeaks++;
          
      }
     
    if(NPeaks==2){
    for (Int_t p = 0; p < NPeaks; p++){
        printf ( "%f ", wyniki[p]);
    }
        //return wyniki;
    }
    
    
      //printf ( "Found %d useful peaks to fit\n", NPeaks);
     


    
    
}

TH1I* getHistoFromTxtFile (TString fileName) {
	
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
	
	while (myStr >> value_temp) value.push_back(value_temp);
	myStr.close (); 
   
	//wpisywanie z wektorów do histogramu 
	
	Int_t nChannels = value.size();
	if (verbatim) std::cout << "\nLoaded file " << fileName << "; nChannels = " << nChannels << ". \n";
	
	TH1I* Spec_uncal = new TH1I ("Spec_uncal", "myhist", nChannels, 0.5, nChannels+0.5);
	for (int i=0; i<nChannels; i++) Spec_uncal->SetBinContent(i+1, value[i]);
    
	//editing the uncalibrated spectrum 
	
	Spec_uncal->SetTitle("Uncalibrated spectrum '" + fileName + "'");
	
	return Spec_uncal;
	
}

void getHistosFromTxtFileList (TString fileList) {
	
	//setup
	
	Bool_t verbatim = true; //choose whether you want to see all debug prompts
	TString fileNameRoot = "output.root";
	if (verbatim) std::cout << "hello" << std::endl;
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
		TH1I* histo = getHistoFromTxtFile(TString(file));
        
        TString hname = TString(file);
        hname.ReplaceAll(".txt", "");
        hname.ReplaceAll("100", "hist100");

        histo->SetName(hname);
        
        dop(histo);
        
        //Double_t result = dop(histo);
        
        //for (Int_t p = 0; p < 2; p++){
       //     printf ( "%f ", result[p]);
       // }
        
        
        histo->Write();
		
	}
	
	myStr.close (); 
	fout->Close();
	
}

int callibrator (TString listName) {
    
    TStopwatch t;
    t.Start();

	std::cout << "Hello world! This is the callibrator!" << std::endl;
	
	getHistosFromTxtFileList(listName);
    
    t.Stop();
    t.Print();
    
	return 0;
	
}
