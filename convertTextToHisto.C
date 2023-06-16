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
	std::vector < double > controlValues(16);
	controlValues[0] = 0.0;
	controlValues[1] = 0.0;
	controlValues[2] = 0.0;
	controlValues[3] = 0.0;
	controlValues[4] = 0.0;
	controlValues[5] = 0.0;
	controlValues[6] = 0.0;
	controlValues[7] = 0.0;
	controlValues[8] = 0.0;
	controlValues[9] = 0.0;
	controlValues[10] = 0.0;
	controlValues[11] = 0.0;
	controlValues[12] = 0.0;
	controlValues[13] = 0.0;
	controlValues[14] = 0.0;
	controlValues[15] = 0.0;

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

	Histogram.Fit("dopasowanie1","MIRQ","",Posit[0]-fitrange, Posit[0]+fitrange); //MIQR
	Histogram.Fit("dopasowanie2","MIRQ+","",Posit[1]-fitrange,Posit[1]+fitrange); //MIQR

	Double_t amplituda1 = dopasowanie1->GetParameter(0);
	Double_t amplituda2 = dopasowanie2->GetParameter(0);

	Double_t srednia1 = dopasowanie1->GetParameter(1);
	Double_t srednia2 = dopasowanie2->GetParameter(1);

	Double_t bladSredniej1 = dopasowanie1->GetParError(1);
	Double_t bladSredniej2 = dopasowanie2->GetParError(1);

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
	controlValues[3] = chi2NDF2;
	controlValues[4] = amplituda1;
	controlValues[5] = amplituda2;
	controlValues[6] = srednia1;
	controlValues[7] = srednia2;
	controlValues[8] = bladSredniej1;
	controlValues[9] = bladSredniej2;

	std::cout << "Found position: " << dopasowanie1->GetParameter(1) << " +/- " << dopasowanie1->GetParError(1) << "\n";
	std::cout << "Found position: " << dopasowanie2->GetParameter(1) << " +/- " << dopasowanie2->GetParError(1) << "\n";

	// Histogram.Write();

	// obliczanie kalibracji

	Double_t y1 = 279.01;
	Double_t y2 = 547.5;

	Double_t x1 = dopasowanie1->GetParameter(1);
	Double_t x2 = dopasowanie2->GetParameter(1);

	Double_t calibFunctionSlope = (y1 - y2)/(x1 - x2);
	Double_t calibFunctionYIntercept  = (-1.0) * calibFunctionSlope * x1 + y1;

	//---------------------------------------------------------------------------
	//dopasowanie trzeciego maksimum
	//---------------------------------------------------------------------------
	// Double_t E3_init = 367.62109;
	// Double_t sigma3_init = 1.8229187;
	// Double_t E3_init = 510.58772;
	// Double_t sigma3_init = 3.1642358;
	Double_t E3_init = 442.50045;
	Double_t sigma3_init = 1.9985245;
	Double_t chan_E3_init = (E3_init-calibFunctionYIntercept)/calibFunctionSlope;
	Double_t chan_sigma3_init = (sigma3_init-calibFunctionYIntercept)/calibFunctionSlope;

	Double_t f1 = 1;
	Double_t f2 = 1;
	TF1* dopasowanie3 = new TF1("dopasowanie3", gausswithlinearbkg, chan_E3_init - f1*chan_sigma3_init, chan_E3_init + f2*chan_sigma3_init, 5);

	dopasowanie3->SetParameters(1000, chan_E3_init, chan_sigma3_init, 0.0000000001, 0.00001);
	dopasowanie3->SetParNames("Amplituda", "Srednia", "Sigma", "A", "B");

	Histogram.Fit("dopasowanie3","MIRQ+","",chan_E3_init - f1*chan_sigma3_init, chan_E3_init + f2*chan_sigma3_init);

	Double_t amplituda3 = dopasowanie3->GetParameter(0);

	Double_t srednia3 = dopasowanie3->GetParameter(1);

	Double_t bladSredniej3 = dopasowanie3->GetParError(1);

	Double_t sigma3 = dopasowanie3->GetParameter(2);

	Double_t chi23 = dopasowanie3->GetChisquare();

	Double_t NDF3 = dopasowanie3->GetNDF();

	Double_t chi2NDF3 = chi23 / NDF3;

	controlValues[10] = sigma3;
	controlValues[11] = chi2NDF3;
	controlValues[12] = amplituda3;
	controlValues[13] = srednia3;
	controlValues[14] = bladSredniej3;
	controlValues[15] = calibFunctionSlope*srednia3 + calibFunctionYIntercept;

	std::cout << "Found position: " << dopasowanie3->GetParameter(1) << " +/- " << dopasowanie3->GetParError(1) << "\n";

	Histogram.Write();
	//---------------------------------------------------------------------------

	TString name = Histogram.GetName();
    name.ReplaceAll(".txt", ".cal");

	ofstream calibOut(name);
	calibOut << calibFunctionYIntercept << "\n";
	calibOut << calibFunctionSlope << "\n";
	//calibOut << srednia3 << "\n";
    calibOut.close();

		name.ReplaceAll(".cal",".3pik");
	ofstream pik3out(name);
		pik3out << srednia3 << "\n";
	pik3out.close();

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
	std::vector < double > sigma3vec;
	std::vector < double > chi2NDF1vec;
	std::vector < double > chi2NDF2vec;
	std::vector < double > chi2NDF3vec;
	std::vector < double > amplituda1vec;
	std::vector < double > amplituda2vec;
	std::vector < double > amplituda3vec;
	std::vector < double > srednia1vec;
	std::vector < double > srednia2vec;
	std::vector < double > srednia3vec;
	std::vector < double > bladSredniej1vec;
	std::vector < double > bladSredniej2vec;
	std::vector < double > bladSredniej3vec;
	std::vector < double > energia3vec;
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
		amplituda1vec.push_back(controlValues[4]);
		amplituda2vec.push_back(controlValues[5]);
		srednia1vec.push_back(controlValues[6]);
		srednia2vec.push_back(controlValues[7]);
		bladSredniej1vec.push_back(controlValues[8]);
		bladSredniej2vec.push_back(controlValues[9]);

		sigma3vec.push_back(controlValues[10]);
		chi2NDF3vec.push_back(controlValues[11]);
		amplituda3vec.push_back(controlValues[12]);
		srednia3vec.push_back(controlValues[13]);
		bladSredniej3vec.push_back(controlValues[14]);

		energia3vec.push_back(controlValues[15]);
	}

	Double_t NHistos = histos.size();

	TH1F* sigma1 = new TH1F("sigma1", "sigma1", NHistos, 0.5, NHistos+0.5);
	TH1F* sigma2 = new TH1F("sigma2", "sigma2", NHistos, 0.5, NHistos+0.5);
	TH1F* chi2NDF1 = new TH1F("chi2NDF1", "chi2NDF1", NHistos, 0.5, NHistos+0.5);
	TH1F* chi2NDF2 = new TH1F("chi2NDF2", "chi2NDF2", NHistos, 0.5, NHistos+0.5);
	TH1F* amplituda1 = new TH1F("amplituda1", "amplituda1", NHistos, 0.5, NHistos+0.5);
	TH1F* amplituda2 = new TH1F("amplituda2", "amplituda2", NHistos, 0.5, NHistos+0.5);
	TH1F* srednia1 = new TH1F("srednia1", "srednia1", NHistos, 0.5, NHistos+0.5);
	TH1F* srednia2 = new TH1F("srednia2", "srednia2", NHistos, 0.5, NHistos+0.5);

	TH1F* sigma3 = new TH1F("sigma3", "sigma3", NHistos, 0.5, NHistos+0.5);
	TH1F* chi2NDF3 = new TH1F("chi2NDF3", "chi2NDF3", NHistos, 0.5, NHistos+0.5);
	TH1F* amplituda3 = new TH1F("amplituda3", "amplituda3", NHistos, 0.5, NHistos+0.5);
	TH1F* srednia3 = new TH1F("srednia3", "srednia3", NHistos, 0.5, NHistos+0.5);

	TH1F* energia3 = new TH1F("energia3", "energia3", NHistos, 0.5, NHistos+0.5);

	TH1F* energia3_hist = new TH1F("energia3_hist", "energia3_hist", 100, 440, 445);

	for (int i=0; i<NHistos; i++) {

		sigma1->SetBinContent(i+1, sigma1vec[i]);
		sigma2->SetBinContent(i+1, sigma2vec[i]);
		chi2NDF1->SetBinContent(i+1, chi2NDF1vec[i]);
		chi2NDF2->SetBinContent(i+1, chi2NDF2vec[i]);
		amplituda1->SetBinContent(i+1, amplituda1vec[i]);
		amplituda2->SetBinContent(i+1, amplituda2vec[i]);
		srednia1->SetBinContent(i+1, srednia1vec[i]);
		srednia2->SetBinContent(i+1, srednia2vec[i]);
		srednia1->SetBinError(i+1, bladSredniej1vec[i]);
		srednia2->SetBinError(i+1, bladSredniej2vec[i]);

		sigma3->SetBinContent(i+1, sigma3vec[i]);
		chi2NDF3->SetBinContent(i+1, chi2NDF3vec[i]);
		amplituda3->SetBinContent(i+1, amplituda3vec[i]);
		srednia3->SetBinContent(i+1, srednia3vec[i]);
		srednia3->SetBinError(i+1, bladSredniej3vec[i]);

		energia3->SetBinContent(i+1,energia3vec[i]);
		energia3_hist->Fill(energia3vec[i]);


		sigma1->GetXaxis()->SetBinLabel(i+1, names[i]);
		sigma2->GetXaxis()->SetBinLabel(i+1, names[i]);
		chi2NDF1->GetXaxis()->SetBinLabel(i+1, names[i]);
		chi2NDF2->GetXaxis()->SetBinLabel(i+1, names[i]);
		amplituda1->GetXaxis()->SetBinLabel(i+1, names[i]);
		amplituda2->GetXaxis()->SetBinLabel(i+1, names[i]);
		srednia1->GetXaxis()->SetBinLabel(i+1, names[i]);
		srednia2->GetXaxis()->SetBinLabel(i+1, names[i]);

		sigma3->GetXaxis()->SetBinLabel(i+1, names[i]);
		chi2NDF3->GetXaxis()->SetBinLabel(i+1, names[i]);
		amplituda3->GetXaxis()->SetBinLabel(i+1, names[i]);
		srednia3->GetXaxis()->SetBinLabel(i+1, names[i]);

		energia3->GetXaxis()->SetBinLabel(i+1, names[i]);

	}

	myStr.close ();
	sigma1->Write();
	sigma2->Write();
	chi2NDF1->Write();
	chi2NDF2->Write();
	amplituda1->Write();
	amplituda2->Write();
	srednia1->Write();
	srednia2->Write();

	sigma3->Write();
	chi2NDF3->Write();
	amplituda3->Write();
	srednia3->Write();

	energia3->Write();
	energia3_hist->Write();

	fout->Close();
	return histos;

}

int convertTextToHisto () {

	TStopwatch t;
	t.Start();

	TString listName = "list.list";
	std::cout << "Selected listName '" << listName << "'\n";

	std::vector<TH1I> histos = getHistosFromTxtFileList(listName);
	Int_t nHistos = histos.size();

	Double_t time = t.RealTime();
	std::cout << "Done analyzing " << nHistos << " spectra in " << time << " seconds. Av = " << time/nHistos << " seconds/spec. \n";
	return 0;

}
