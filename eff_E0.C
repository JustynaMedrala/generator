#include <TRandom3.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>

TH1F* hist_Edep =nullptr;
TH1F* hist_Edep_smear =nullptr;
TH1F* hist_ops_energy =nullptr;
TGraph* graph_eff =nullptr;
TGraph* graph_eff_smear =nullptr;
TGraph* graph_diff =nullptr;
TGraph* graph_E =nullptr;
TGraph* graph_bkg =nullptr;
TGraph* graph_pur =nullptr;

void eff_E0(){
  TRandom3 gen;
  gen.SetSeed(0);
  TF1 *f = new TF1("f","1/(1+[0]*(1-cos(x)))^2*(1+[0]*(1-cos(x))+1/(1+[0]*(1-cos(x)))-(sin(x))^2)/2",0,TMath::Pi());
  double E0, theta, Ek, Edep, eff, eff_smear, sigma, Energy, bkg = 0;
  int bin_min = 30;
  int bin_max = 80;
  int N = 21;
  int M = int(1e6);
  float m[N];
  double me = 511;
  double mu = 0;
 
  for(int i=0; i<N; i++){
    m[i] = float(i)/float(N-1);
  }


  TFile *hfile = 0;
  hfile = TFile::Open("energies.root","READ");
  TTree *tree = hfile->Get<TTree>("T");
  tree->SetBranchAddress("Energy",&Energy);

  TCanvas *c1 = new TCanvas();  
  graph_eff = new TGraph(N);
  graph_eff_smear = new TGraph(N);
  graph_diff = new TGraph(N);
  graph_E = new TGraph(N);
  graph_bkg = new TGraph(N);
  graph_pur = new TGraph(N);

  int nentries = tree->GetEntries();
  hist_ops_energy = new TH1F("hist_ops_energy", "hist_ops_energy",100,0,511);
  for(int i = 0; i < nentries; i++){
    tree->GetEntry(i);
    hist_ops_energy->Fill(Energy);
  }

  
  for(int i = 0; i < N; i++){
    mu = 2*me*m[i];
    E0 = me*(1-mu*mu/me/me/4.);//sqrt(mu*(4*me-mu))-mu;
    E0 = E0/511.;
    cout<<E0<<endl;
    f->SetParameter(0, E0);
    hist_Edep = new TH1F("hist_Edep","hist_Edep",100,0,511);
    hist_Edep_smear = new TH1F("hist_Edep_with_smear", "hist_Edep_with_smear",100,0,511);
    for(int j = 0; j < M; j++){
      theta = f->GetRandom();
      Ek = E0/(1+E0*(1-cos(theta)));
      Edep = E0 - Ek;
      hist_Edep->Fill(Edep*511);     
 
      sigma = sqrt(Edep)*0.044*0.511;
      Edep = gen.Gaus(Edep*0.511, sigma)/0.511;

      hist_Edep_smear->Fill(Edep*511);
    }
    
    eff = hist_Edep->Integral(bin_min, bin_max)/hist_Edep->Integral();
    eff_smear = hist_Edep_smear->Integral(bin_min, bin_max)/hist_Edep_smear->Integral();
    bkg = hist_ops_energy->Integral(bin_min, bin_max)/hist_Edep_smear->Integral();
    graph_E->SetPoint(i, m[i], E0);
    graph_eff->SetPoint(i, m[i], eff*100);
    graph_eff_smear->SetPoint(i, m[i], eff_smear*100);
    graph_diff->SetPoint(i, m[i], (eff-eff_smear)*100);
    graph_bkg->SetPoint(i, m[i], bkg*100);
    graph_pur->SetPoint(i, m[i], eff_smear/(eff_smear+bkg));
    hist_Edep->Delete();
    hist_Edep_smear->Delete();
  }
  //Drawing
  ////////////////////////////////////////////////////////////////////////
  graph_eff->SetTitle("Efficiency vs m_{u}");
  graph_eff->GetXaxis()->SetTitle("#frac{m_{u}}{m_{e}}");
  graph_eff->GetYaxis()->SetTitle("Efficiency [%]");
  graph_eff->SetMarkerStyle(20);
  graph_eff->SetMarkerSize(1.);
  graph_eff_smear->SetMarkerStyle(20);
  graph_eff_smear->SetMarkerColor(kRed);
  graph_bkg->SetMarkerStyle(20);
  graph_bkg->SetMarkerColor(kGreen);
  graph_eff->Draw("AP");
  graph_eff_smear->Draw("same P");
  graph_bkg->Draw("same P");

  auto legend = new TLegend(0.65,0.2,0.9,0.4);
  legend->AddEntry(graph_eff,"energy without smear","p");
  legend->AddEntry(graph_eff_smear,"energy with smear","p");
  legend->AddEntry(graph_bkg, "background (o-Ps)", "p");
  legend->Draw();

  
  c1->SaveAs("../plots/Efficiency.png");

  ////////////////////////////////////////////////////////////////////////
  graph_diff->SetTitle("Difference in efficiencies vs m_{u}");
  graph_diff->GetXaxis()->SetTitle("#frac{m_{u}}{m_{e}}");
  graph_diff->GetYaxis()->SetTitle("Difference in efficiencies [%]");
  graph_diff->SetMarkerStyle(20);
  graph_diff->SetMarkerSize(1.);
  graph_diff->Draw("AP");
  
  c1->SaveAs("../plots/Difference_efficiency.png");
  

  ////////////////////////////////////////////////////////////////////////
  graph_E->SetTitle("Energy vs m_{u}");
  graph_E->GetXaxis()->SetTitle("#frac{m_{u}}{m_{e}}");
  graph_E->GetYaxis()->SetTitle("Energy #frac{E_{#gamma}}{E_{e}}");
  graph_E->SetMarkerStyle(20);
  graph_E->SetMarkerSize(1.);
  graph_E->Draw("AP");
  
  c1->SaveAs("../plots/Energy_vs_mu.png");

  ////////////////////////////////////////////////////////////////////////
  graph_pur->SetTitle("Purity vs m_{u}");
  graph_pur->GetXaxis()->SetTitle("#frac{m_{u}}{m_{e}}");
  graph_pur->GetYaxis()->SetTitle("Purity #frac{N_{sig}}{N_{sig}+N_{bkg}}");
  graph_pur->SetMarkerStyle(20);
  graph_pur->SetMarkerSize(1.);
  graph_pur->Draw("AP");

  c1->SaveAs("../plots/Purity_vs_mu.png");
  
}
