
#include <iostream>
#include <fstream>
#include <vector>


//wypisuje nazwe pliku .cal ze sciezka
void sciezka(TString s){
    
    
    s.ReplaceAll("Run","/Run");
    s.ReplaceAll("no","/MCA/");
    s.ReplaceAll("Ge","/Ge");
    s.ReplaceAll(".txt",".cal");
    
    
    TString filename = "home/analecz/tv/data/128Cs_HIL/";
    filename.Append(s);
    
    std::cout<<filename<<std::endl;
    //wypisuje: home/analecz/tv/data/128Cs_HIL/100um/Run11/MCA/0005/Ge01.cal
    

}

//Wypisuje wszystkie parametry (numer runu, detektora itd.)
//Wersja 1, char, sprintf
void wersja1(TString s){
    
    s.ReplaceAll(".txt","");
    Int_t slen=s.Length();
    
    char um[10], run[10], no[10], ge[10];
    
    if(slen==20){sprintf(um,"%c%c%c",s[0],s[1],s[2]);};
    if(slen==19){sprintf(um,"%c%c",s[0],s[1]);};
    if(slen==18){sprintf(um,"%c",s[0]);};
    
    sprintf(run,"%c%c",s[slen-12],s[slen-11]);
    sprintf(no,"%c%c%c%c",s[slen-8],s[slen-7],s[slen-6],s[slen-5]);
    sprintf(ge,"%c%c",s[slen-2],s[slen-1]);
    
    
    std::cout<<"Odleglosc: "<<um<<" um"<<std::endl;
    std::cout<<"Numer runu: "<<run<<std::endl;
    std::cout<<"Numer: "<<no<<std::endl;
    std::cout<<"Numer detektora: "<<ge<<std::endl;
    
   
    
}

//Wypisuje wszystkie parametry (numer runu, detektora itd.), w wersji TString lub Int_t
//Wersja 2, TString, Append
void wersja2(TString s){
    
    s.ReplaceAll(".txt","");
    Int_t slen=s.Length();
    
    TString um, run, no, ge;
    
    um=s[0];
    
    TString u = "u";
    
    TString s1 = s[1];
    TString s2 = s[2];
    TString s3 = s[3];
    
    if(s3.CompareTo(u)==0){um.Append(s1); um.Append(s2);};
    if(s2.CompareTo(u)==0){um.Append(s1);};
    if(s1.CompareTo(u)==0){};
    Int_t um_i = um.Atoi();
    
    
    run=s[slen-12];
    run.Append(s[slen-11]);
    Int_t run_i = run.Atoi();
    
    no=s[slen-8];
    no.Append(s[slen-7]);
    no.Append(s[slen-6]);
    no.Append(s[slen-5]);
    Int_t no_i = no.Atoi();
    
    ge=s[slen-2];
    ge.Append(s[slen-1]);
    Int_t ge_i = ge.Atoi();
    
    
    //TString
    std::cout<<"Odleglosc: "<<um<<" um"<<std::endl;
    std::cout<<"Numer runu: "<<run<<std::endl;
    std::cout<<"Numer: "<<no<<std::endl;
    std::cout<<"Numer detektora: "<<ge<<std::endl;
    
    
    //Int_t
    std::cout<<"Odleglosc: "<<um_i<<" um"<<std::endl;
      std::cout<<"Numer runu: "<<run_i<<std::endl;
      std::cout<<"Numer: "<<no_i<<std::endl;
      std::cout<<"Numer detektora: "<<ge_i<<std::endl;
    
    
}





int mojprogram() {
    
    TString s = "100umRun11no0005Ge01.txt";
    std::cout<<s<<std::endl;
    

    //sciezka(s); //wypisuje nazwe pliku .cal ze sciezka
    
    //Wypisują wszystkie parametry (numer runu, detektora itd.):
    //wersja1(s); //char, sprintf
    wersja2(s); //TString, Append; wypisuje w wersji TString lub Int_t
    

    

    return 0;
    
}


