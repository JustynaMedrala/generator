#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <TRandom3.h>
#include <TF1.h>
#include <TMath.h>


using namespace std;

TH1D* hist_energy = nullptr;

vector<string> split(const string &str, char sep){
  istringstream ss(str);
  vector<string> v;
  string s;
  while(getline(ss, s, sep))
    v.push_back(s);
  return v;
}



void to_root()
{
  double Energy;
  TRandom3 gen;
  gen.SetSeed(0);
  TF1 *f = new TF1("f","1/(1+[0]*(1-cos(x)))^2*(1+[0]*(1-cos(x))+1/(1+[0]*(1-cos(x)))-(sin(x))^2)/2",0,TMath::Pi());

  ifstream inFile;
  inFile.open("energies.txt");
  string energy;
  
  TFile *hfile = 0;
  hfile = TFile::Open("energies.root","RECREATE");
  TTree *tree = new TTree("T","T");
  tree->Branch("Energy",&Energy,"Energy/D");
  tree->Branch("E0", &E0, "E0/D");

  double E0, Ek, Edep, theta;
  int ind = 0;
  hist_energy = new TH1D(" ", " ", 100, 0, 520);

   if(inFile.is_open())                                                               
   {
     while(getline(inFile, energy))                                                           
     {
       //getline(inFile, energy);
       vector<string> tmp = split(energy, ',');
       ind++;
       if(ind%10000==0) cout<<double(ind)/double(100000)*100<<"%"<<endl;
       for(int i=0; i<3; i++){
         E0 = stod(tmp[i]);         
         f->SetParameter(0, E0);
         theta = f->GetRandom(); 
         Ek = E0/(1+E0*(1-cos(theta)));
         Edep = E0-Ek;
         //cout<<E0<<", "<<Edep<<endl;
         hist_energy->Fill(Edep);
         Energy = Edep;
         tree->Fill();
       }  
    }
  }
  tree->Write();


  TCanvas *c1 = new TCanvas();
  hist_energy->GetXaxis()->SetTitle("E_{dep} [keV]");
  hist_energy->Draw();
  c1->SaveAs("../plots/energy_ops.png");
}
