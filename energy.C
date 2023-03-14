#include <TRandom3.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>

TH1F* hist_Edep =nullptr;
TH1F* hist_Edep_s =nullptr;
TH1F* hist_Edep_deex =nullptr;
TH1F* hist_Edep_deex_s =nullptr;
TH1F* hist_Ek =nullptr;
TH1F* hist_ops_energy =nullptr;

void energy(){
  TRandom3 gen;
  gen.SetSeed(0);
  TF1 *f = new TF1("f","1/(1+(1-cos(x)))^2*(1+(1-cos(x))+1/(1+(1-cos(x)))-(sin(x))^2)",0,TMath::Pi()); 
  double theta = f->GetRandom(); 
  double E0 = 1;
  double Ek = 0;
  double Edep = E0 - Ek;
  int N = int(1e6);
  double sigma = 0;
  double Energy = 0;

  hist_Edep = new TH1F("hist_Edep","hist_Edep",100,0,520);
  hist_Edep_s = new TH1F("hist_Edep_smear", "hist_Edep_smear", 100, 0, 1274);
  hist_Edep_deex = new TH1F("hist_Edep_deex", "hist_Edep_deex", 100, 0, 1274);
  hist_Edep_deex_s = new TH1F("hist_Edep_deex_smear", "hist_Edep_deex_smear", 100, 0, 1274);
  hist_Ek = new TH1F("hist_Ek","hist_Ek",100, 0,1.5);
  hist_ops_energy = new TH1F("hist_ops_energy", "hist_ops_energy", 100, 0, 1274);

  TFile *hfile = 0;
  hfile = TFile::Open("energies.root","READ");
  TTree *tree = hfile->Get<TTree>("T");
  tree->SetBranchAddress("Energy",&Energy);

  int nentries = tree->GetEntries();

  for(int i = 0; i < nentries; i++){
    tree->GetEntry(i);
    hist_ops_energy->Fill(Energy);
  }

  for(int i = 0; i<N; i++){
    theta = f->GetRandom();
    Ek = E0/(1+E0*(1-cos(theta)));
    Edep = (E0-Ek);
    hist_Edep->Fill(Edep*511);
    sigma = sqrt(Edep)*0.044*0.511;
    Edep = gen.Gaus(Edep*0.511, sigma)/0.511;
    hist_Edep_s->Fill(Edep*511);
    hist_Ek->Fill(Ek);
  }
  /////////////////////////////////
  for(int i = 0; i<N; i++){
    theta = f->GetRandom();
    E0 = 1274./511.;
    Ek = E0/(1+E0*(1-cos(theta)));
    Edep = (E0-Ek);
    hist_Edep_deex->Fill(Edep*511);
    sigma = sqrt(Edep)*0.044*0.511;
    Edep = gen.Gaus(Edep*0.511, sigma)/0.511;
    hist_Edep_deex_s->Fill(Edep*511);
  }

  TCanvas *c = new TCanvas();
  gStyle->SetOptStat(0);
  ////////////////////////////////////
  f->SetTitle("  ");
  f->GetXaxis()->SetTitle("#theta");
  f->GetYaxis()->SetTitle("f(#theta)");
  f->Draw();
  c->SaveAs("../plots/f(theta).png");
  ///////////////////////////////////
  auto *l_tres = new TLine(50,0,50,0.06);
  auto *l_min = new TLine(0.3*511, 0, 0.3*511, 0.06);
  auto *l_max = new TLine(0.8*511, 0, 0.8*511, 0.06);
  ///////////////////////////////////
  hist_ops_energy->SetTitle("   ");
  hist_ops_energy->GetXaxis()->SetTitle("E_{dep}");
  hist_ops_energy->SetLineColor(kGreen);
  hist_Edep_deex_s->SetLineColor(kViolet);
  l_tres->SetLineColor(kRed);
  hist_Edep_deex_s->Scale(1.0/hist_Edep_deex_s->GetEntries());
  hist_ops_energy->Scale(1.0/hist_ops_energy->GetEntries());
  hist_Edep_s->Scale(1.0/hist_Edep->GetEntries());
  hist_ops_energy->Draw("hist");
  hist_Edep_s->Draw("same hist");
  hist_Edep_deex_s->Draw("same hist");
  l_tres->Draw("same");
  l_min->Draw("same");
  l_max->Draw("same");

  auto legend = new TLegend(0.65,0.7,0.9,0.9);
  legend->AddEntry(hist_Edep_s,"energy with smear","l");
  legend->AddEntry(hist_ops_energy,"energy for o-Ps","l");
  legend->AddEntry(hist_Edep_deex_s, "deexcitation energy", "l");
  legend->AddEntry(l_tres, "detector threshold", "l");
  legend->Draw(); 

  c->SaveAs("../plots/hist_Edep_ws.png");
  ///////////////////////////////////
  hist_Edep_deex->SetLineColor(kViolet);
  hist_Edep->Scale(1.0/hist_Edep->GetEntries());
  hist_Edep_deex->Scale(1.0/hist_Edep_deex->GetEntries());
  hist_ops_energy->Draw("hist");
  hist_Edep->Draw("same hist");
  hist_Edep_deex->Draw("same hist");
  l_tres->Draw("same");
  l_min->Draw("same");
  l_max->Draw("same");

  auto legend1 = new TLegend(0.65,0.7,0.9,0.9);
  legend1->AddEntry(hist_Edep,"energy with smear","l");
  legend1->AddEntry(hist_ops_energy,"energy for o-Ps","l");
  legend->AddEntry(hist_Edep_deex, "deexcitation energy", "l");
  legend1->AddEntry(l_tres, "detector threshold", "l");
  legend1->Draw();

  c->SaveAs("../plots/hist_Edep_wos.png");
  //////////////////////////////////
  hist_Ek->SetTitle("   ");
  hist_Ek->GetXaxis()->SetTitle("E_{#gamma'}");
  hist_Ek->Draw();
  c->SaveAs("../plots/hist_Ek.png");
}
