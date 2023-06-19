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
TH1D* hist_e = nullptr;
TH1D* hist_E0 = nullptr;

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
  double Energy, E_0, Energy_att;
  TRandom3 gen;
  gen.SetSeed(0);
  TF1 *f = new TF1("f","1+[0]*(1-cos(x))+1/(1+[0]*(1-cos(x)))-(sin(x))^2",0,TMath::Pi());

  //double E0;

  ifstream inFile;
  inFile.open("energies.txt");
  string energy;
  
  TFile *hfile = 0;
  hfile = TFile::Open("energies.root","RECREATE");
  TTree *tree = new TTree("T","T");
  tree->Branch("Energy",&Energy,"Energy/D");
  tree->Branch("E_0", &E_0, "E_0/D");
  tree->Branch("Energy_att", &Energy_att, "Energy_att/D");

  int ind = 0;
  double E[1999], a[1999], eff[1999];
  string atten;
  fstream atFile;
  atFile.open("efficiency_vs_energy.txt");

  if(atFile.is_open()){
    while(getline(atFile, atten))
    {
      vector<string> tmp = split(atten, ',');
      E[ind] = stod(tmp[0]);
      a[ind] = stod(tmp[1]);
      eff[ind] = pow(1-((1-exp(-a[ind]*2))*(1 - sin(atan(2)/2.)*sin(atan(2)/2.)*2.)),2);
      ind++;
    }
  }


  double E0, Ek, Edep, theta;
  ind = 0;
  hist_energy = new TH1D("h_att", "h_att", 511, 0, 511);
  hist_e = new TH1D("h_Edep", "h_Edep", 511, 0, 511);
  hist_E0 = new TH1D("h_E0", "h_E0", 511, 0, 520);

   if(inFile.is_open())                                                               
   {
     while(getline(inFile, energy))                                                           
     {
       //getline(inFile, energy);
       vector<string> tmp = split(energy, ',');
       ind++;
       if(ind%10000==0) {cout<<double(ind)/double(100000)*100<<"%"<<endl;}
       for(int i=0; i<3; i++){
         E0 = stod(tmp[i])/511.;         
         f->SetParameter(0, E0);
         theta = f->GetRandom(); 
         Ek = E0/(1+E0*(1-cos(theta)));
         Edep = E0-Ek;
         //cout<<E0<<", "<<Edep<<", "<<theta<<endl;
         Energy = Edep*511.;
         E_0 = E0*511.;
         hist_E0->Fill(E_0);
         hist_e->Fill(Energy);
         if(E_0 > 1.) hist_energy->Fill(Energy, eff[int(E_0)-1]);
         else hist_energy->Fill(Energy, eff[0]);
         //cout<<eff[int(E_0)]<<endl;
         tree->Fill();
       }  
    }
  }
  tree->Write();
  hist_e->Write();
  hist_energy->Write();
  hist_E0->Write();

  TCanvas *c1 = new TCanvas();
  gStyle->SetOptStat(0);

  hist_e->SetTitle(" ");
  hist_e->GetXaxis()->SetTitle("E_{dep} [keV]");
  hist_e->SetLineColor(kRed);
  hist_e->Draw("hist");
  hist_energy->Draw("same hist");
  auto legend = new TLegend(0.65,0.7,0.9,0.9);
  legend->AddEntry(hist_e,"deposited energy","l");
  legend->AddEntry(hist_energy,"scaled distribution","l");
  legend->Draw();
  c1->SaveAs("../plots/energy_ops.pdf");

  hist_E0->GetXaxis()->SetTitle("E [keV]");
  hist_E0->SetTitle(" ");
  hist_E0->Draw();
  c1->SaveAs("../plots/E0_hist.pdf");
}
