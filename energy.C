#include <TStyle.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TFile.h>
#include <TLine.h>
#include <TF1.h>
#include <TMath.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <iostream>
#include <memory>

using namespace std;

TH1F* hist_Edep =nullptr;
TH1F* hist_Edep_s =nullptr;
TH1F* hist_Edep_deex =nullptr;
TH1F* hist_Edep_deex_s =nullptr;
TH1F* hist_Ek =nullptr;
TH1F* hist_ops_energy =nullptr;
TH1F* hist_Edep_color =nullptr;
TH1F* hist_Edep_color_s =nullptr;

void part_of_hist(TH1F* hist, TH1F* hist_part, double min_value, double max_value){
  int min = hist->FindBin(min_value);
  int max = hist->FindBin(max_value);
  cout<<min<<", "<<max<<endl;
  for(int i = min; i < max; i++){
    //cout<<hist->GetBinContent(i)<<endl;
    hist_part->AddBinContent(i, hist->GetBinContent(i));
    //cout<<hist_part->GetBinContent(i)<<endl;
  }
  hist_part->Draw("hist");
  //return hist_part;
}

double en_dep(TF1* f, double E0){
    double Ek, theta;
    f->SetParameter(0, E0);
    theta = f->GetRandom();
    Ek = E0/(1+E0*(1-cos(theta)));
    return (E0-Ek);
}



double helperTerm(double E0, double theta)
{
  return  1./(1.+ E0*(1- TMath::Cos(theta)));
}

double get_Ek(double E0, double theta)
{
  return E0/(1+E0*(1-cos(theta)));
}
void test_energy_function()
{
  auto KN2 = [&](double*x, double* p){ return helperTerm(p[0], x[0]) +1./helperTerm(p[0], x[0]) - TMath::Sin(x[0])*TMath::Sin(x[0]);};
  TH1F* hist_Edep = new TH1F("hist_Edep","hist_Edep",100,0,1274);
  hist_Edep->SetCanExtend(TH1::kAllAxes);
  TH1F* hist_theta = new TH1F("hist_theta","hist_theta",100,0,TMath::Pi());
  hist_theta->SetCanExtend(TH1::kAllAxes);
  /// KN = dsigma/dOmega
  TF1 *func0 = new TF1("KN","1/(1+[0]*(1-cos(x)))^2*(1+[0]*(1-cos(x))+1/(1+[0]*(1-cos(x)))-(sin(x))^2)",0,TMath::Pi());
  /// KN = dsigma/dTheta -> Jacobian included
  TF1 *func1 = new TF1("KN","1/(1+[0]*(1-cos(x))) + (1+[0]*(1-cos(x)))-(sin(x))^2",0,TMath::Pi());
  /// KN = dsigma/dTheta -> Jacobian included, the same as before but using predefined lambda
  TF1 *func = new TF1("KN",KN2,0,TMath::Pi(), 1);
  double E0 =1;
  func->SetParameter(0, E0);
  std::cout << func->GetXmax() << std::endl;
  std::cout << func->GetXmin() << std::endl;
  std::cout << func->GetMaximum() << std::endl;
  std::cout << func->GetMinimum() << std::endl;
  double Ek = 0;
  Ek = get_Ek(E0, func->GetXmax()); 
  std::cout <<E0 -Ek  << std::endl;
  Ek = get_Ek(E0, func->GetXmin()); 
  std::cout << E0-Ek  << std::endl;
  double theta =0;
  double Edep =0;
  int nevents = 100000;
  for (int i = 0; i < nevents; i++) {
    theta = func->GetRandom();    
    Edep = E0 - get_Ek(E0, theta);
    hist_Edep->Fill(Edep*511);
    assert(theta>=0);
    assert(theta<= TMath::Pi());
    hist_theta->Fill(theta);
  }
  hist_Edep->Draw();
}

void fill_hist_from_file(const char* filename , TH1* hist_ops_energy)
{
  double energy = 0;
  std::unique_ptr<TFile> hfile( TFile::Open(filename,"READ") );
  TTree *tree = hfile->Get<TTree>("T");
  tree->SetBranchAddress("E_0",&energy);

  int nentries = tree->GetEntries();

  for(int i = 0; i < nentries; i++){
    tree->GetEntry(i);
    hist_ops_energy->Fill(energy);
  }
}


void energy() {
  TRandom3 gen;
  gen.SetSeed(0);
  TF1 *f = new TF1("f","1/(1+[0]*(1-cos(x))) + (1+[0]*(1-cos(x)))-(sin(x))^2",0,TMath::Pi());
  double E0 = 1;
  double Ek = 0;
  double Edep = E0 - Ek;
  int N = int(1);
  double sigma = 0;
  int num_bins = 200;
  hist_Edep = new TH1F("hist_Edep","hist_Edep",num_bins,0,1274);
  hist_Edep_s = new TH1F("hist_Edep_smear", "hist_Edep_smear", num_bins, 0, 1274);
  hist_Edep_deex = new TH1F("hist_Edep_deex", "hist_Edep_deex", num_bins, 0, 1274);
  hist_Edep_deex_s = new TH1F("hist_Edep_deex_smear", "hist_Edep_deex_smear", num_bins, 0, 1274);
  hist_Ek = new TH1F("hist_Ek","hist_Ek",num_bins, 0,1.5);
  hist_ops_energy = new TH1F("hist_ops_energy", "hist_ops_energy", 200, 0, 520);
  hist_Edep_color = new TH1F("hist_Edep_color","hist_Edep_color",num_bins,0,1274);
  hist_Edep_color_s = new TH1F("hist_Edep_color_s","hist_Edep_color_s",num_bins,0,1274);


  //TFile *hfile = 0;
  //hfile = TFile::Open("energies.root","READ");
  //TTree *tree = hfile->Get<TTree>("T");
  //tree->SetBranchAddress("Energy",&Energy);

  //hist_ops_energy = (TH1F*)hfile->Get("h_att");
  fill_hist_from_file("energies.root", hist_ops_energy);
  cout<<"bkg: "<< hist_ops_energy->Integral(0.3*511, 0.8*511)/hist_ops_energy->Integral()<<endl; 

  for(int i = 0; i<N; i++){
    Edep = en_dep(f, E0);
    //cout<<Edep<<endl;
    hist_Edep->Fill(Edep*511);
    sigma = sqrt(Edep)*0.044*0.511;
    Edep = gen.Gaus(Edep*0.511, sigma)/0.511;
    hist_Edep_s->Fill(Edep*511);
    hist_Ek->Fill(Ek);
  }
  /////////////////////////////////
  //cout<<2<<endl;
  /////////////////////////////////
  //cout<<3<<endl;
  for(int i = 0; i<N; i++){
    E0 = 1274./511.;
    Edep = en_dep(f, E0);
    hist_Edep_deex->Fill(Edep*511);
    sigma = sqrt(Edep)*0.044*0.511;
    Edep = gen.Gaus(Edep*0.511, sigma)/0.511;
    hist_Edep_deex_s->Fill(Edep*511);
  }
  /////////////////////////////////
  hist_Edep_deex_s->Scale(1.0/hist_Edep_deex_s->GetEntries());
  //hist_ops_energy->Scale(1.0/hist_ops_energy->GetEntries());
  hist_Edep_s->Scale(1.0/hist_Edep->GetEntries());
  hist_Edep->Scale(1.0/hist_Edep->GetEntries());
  hist_Edep_deex->Scale(1.0/hist_Edep_deex->GetEntries());
  /////////////////////////////////
  part_of_hist(hist_Edep, hist_Edep_color, 0.3*511, 0.8*511);
  part_of_hist(hist_Edep_s, hist_Edep_color_s, 0.3*511, 0.8*511);
  //for(int i=0; i<100; i++)  cout<<hist_Edep_color->GetBinContent(i)<<endl;
  /////////////////////////////////
  TCanvas *c = new TCanvas(" ", " ", 900, 800);
  gStyle->SetOptStat(0);

  ////////////////////////////////////
  f->SetTitle("  ");
  f->GetXaxis()->SetTitle("#theta");
  f->GetYaxis()->SetTitle("f(#theta)");
  f->Draw();
  c->SaveAs("../plots/f(theta).pdf");
  ///////////////////////////////////
  auto *l_tres = new TLine(50,0,50,0.1);
  auto *l_min = new TLine(0.3*511, 0, 0.3*511, 0.1);
  auto *l_max = new TLine(0.8*511, 0, 0.8*511, 0.1);
  l_tres->SetLineWidth(2);
  l_min->SetLineWidth(2);
  l_max->SetLineWidth(2);
  ///////////////////////////////////
  hist_Edep_color_s->SetFillStyle(1001);
  hist_ops_energy->SetTitle("   ");
  hist_ops_energy->GetXaxis()->SetTitle("E [keV]");
  //hist_ops_energy->GetYaxis()->SetRangeUser(0, 0.14);
  hist_ops_energy->SetLineColor(kGreen);
  hist_Edep_color_s->SetFillColorAlpha(kBlack, 0.3);
  hist_Edep_deex_s->SetLineColor(kViolet);
  l_tres->SetLineColor(kRed);
  hist_ops_energy->Draw("hist");
  hist_Edep_s->Draw("same hist");
  //hist_ops_energy->Draw("same hist");
  hist_Edep_deex_s->Draw("same hist");
  hist_Edep_color_s->Draw("same hist");
  l_tres->Draw("same");
  l_min->Draw("same");
  l_max->Draw("same");

  auto legend = new TLegend(0.65,0.7,0.9,0.9);
  legend->AddEntry(hist_Edep_s,"energy with smear","l");
  legend->AddEntry(hist_ops_energy,"energy for o-Ps","l");
  legend->AddEntry(hist_Edep_deex_s, "deexcitation energy", "l");
  legend->AddEntry(l_tres, "detector threshold", "l");
  legend->Draw(); 
  //std::string outPlotFile = "../plots/hist_Edep_ws.pdf";
  std::string outPlotFile = "hist_Edep_ws.pdf";
  c->SaveAs(outPlotFile.c_str());
  ///////////////////////////////////
  hist_Edep_color->SetFillStyle(1001);
  hist_Edep_color->SetFillColorAlpha(kBlack, 0.3);
  hist_Edep_deex->SetLineColor(kViolet);
  //hist_ops_energy->GetYaxis()->SetRangeUser(0, 0.3);
  hist_ops_energy->Draw("hist");
  hist_Edep->Draw("same hist");
  hist_Edep_deex->Draw("same hist");
  hist_Edep_color->Draw("same hist");
  l_tres->Draw("same");
  l_min->Draw("same");
  l_max->Draw("same");

  auto legend1 = new TLegend(0.65,0.7,0.9,0.9);
  legend1->AddEntry(hist_Edep,"energy without smear","l");
  legend1->AddEntry(hist_ops_energy,"energy for o-Ps","l");
  legend->AddEntry(hist_Edep_deex, "deexcitation energy", "l");
  legend1->AddEntry(l_tres, "detector threshold", "l");
  legend1->Draw();

  //std::string outPlotFile = "../plots/hist_Edep_wos.pdf";
  std::string outPlotFile2 = "hist_Edep_wos.pdf";
  c->SaveAs(outPlotFile2.c_str());
  //////////////////////////////////
  hist_Ek->SetTitle("   ");
  hist_Ek->GetXaxis()->SetTitle("E_{#gamma'}");
  hist_Ek->Draw();
  //std::string outPlotFile3 = "../plots/hist_Ek.png";
  std::string outPlotFile3 = "hist_Ek.png";
  c->SaveAs(outPlotFile3.c_str());
  
  hist_Edep_s->SetTitle(" ");
  hist_Edep_s->SetFillColor(kRed);
  hist_Edep_s->SetLineColor(kRed);
  hist_Edep_s->Draw("hist");
  c->SaveAs("E_dep_s.pdf");

  hist_Edep->GetXaxis()->SetTitle("E_{dep} [keV]");
  hist_Edep->SetTitle(" ");
  hist_Edep->SetFillColor(kBlue);
  hist_Edep->SetLineColor(kBlue);
  hist_Edep->Draw("hist");
  c->SaveAs("E_dep.pdf");

  hist_ops_energy->GetYaxis()->SetTitle("Counts");
  hist_ops_energy->SetLineColor(kBlue);
  hist_ops_energy->Draw("hist");
  c->SaveAs("../plots/E_ops.pdf");
}
