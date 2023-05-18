// callibrator.C

#include <iostream>
#include <fstream>
#include <vector>

#include "TH1I.h"
#include "TFile.h"

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
	if (verbatim) std::cout << "Loaded file " << fileName << "; nChannels = " << nChannels << ". \n"; 
	
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
		histo->SetName(TString(file));
		histo->Write();
		
	}
	
	myStr.close (); 
	fout->Close();
	
}

int convertTextToHisto (TString listName) {

	std::cout << "Hello world! This is the callibrator!" << std::endl;
	
	getHistosFromTxtFileList(listName);

	return 0;
	
}
