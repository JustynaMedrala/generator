#include <TRandom3.h>
#include <TF1.h>
#include <TMath.h>
#include <string>
#include <iostream>

TH1F* hist_Edep =nullptr;
TH1F* hist_Edep_smear =nullptr;
TGraph* graph_eff =nullptr;
TGraph* graph_eff_smear =nullptr;
TGraph* graph_diff =nullptr;
TGraph* graph_E =nullptr;
TGraph* graph_br =nullptr;

void eff_mass(){
  TRandom3 gen;
  gen.SetSeed(0);
  TF1 *f = new TF1("f","1/(1+[0]*(1-cos(x)))^2*(1+[0]*(1-cos(x))+1/(1+[0]*(1-cos(x)))-(sin(x))^2)/2",0,TMath::Pi());
  double E0, theta, Ek, Edep, eff, eff_smear, sigma, mu = 0;
  double me = 511;
  int bin_min = 30;
  int bin_max = 80;
  int N = 21;
  int M = int(1e6);
  float m[N];

  for(int i=0; i<N; i++){
    m[i] = float(i)/float(N-1);
  }

  TCanvas *c1 = new TCanvas();  
  graph_eff = new TGraph(N);
  graph_eff_smear = new TGraph(N);
  graph_diff = new TGraph(N);
  graph_E = new TGraph(N);
  graph_br = new TGraph(N);

  for(int i = 0; i < N; i++){
    mu = 2*me*m[i];
    E0 = sqrt(mu*(4*me-mu))-mu;
    f->SetParameter(0, E0);
    hist_Edep = new TH1F("hist_Edep","hist_Edep",100,0,1);
    hist_Edep_smear = new TH1F("hist_Edep_with_smear", "hist_Edep_with_smear",100,0,1);
    for(int j = 0; j < M; j++){
      theta = f->GetRandom();
      Ek = E0/(1+E0*(1-cos(theta)));
      Edep = E0 - Ek;
      hist_Edep->Fill(Edep);     
 
      sigma = sqrt(Edep)*0.044*0.511;
      Edep = gen.Gaus(Edep*0.511, sigma)/0.511;

      hist_Edep_smear->Fill(Edep);
    }
    eff = hist_Edep->Integral(bin_min, bin_max)/hist_Edep->Integral();
    eff_smear = hist_Edep_smear->Integral(bin_min, bin_max)/hist_Edep_smear->Integral();
    //cout<<eff_smear-eff<<endl;
    graph_E->SetPoint(i, m[i], E0);
    graph_eff->SetPoint(i, m[i], eff*100);
    graph_eff_smear->SetPoint(i, m[i], eff_smear*100);
    graph_diff->SetPoint(i, m[i], (eff-eff_smear)*100);
    graph_br->SetPoint(i, m[i], 4*(1-m[i]*m[i]));
    if(i == N-1){
      string tit = "E0 =";
      const char *title = (tit+to_string(E0)).c_str();
      hist_Edep->SetTitle(title);
      hist_Edep->Draw();
      c1->SaveAs("../plots/hist_Edep_max_m.png");}
    hist_Edep->Delete();
    hist_Edep_smear->Delete();
  }
  //Drawing
  ////////////////////////////////////////////////////////////////////////
  graph_eff->SetTitle("Efficiency vs m_{u}");
  graph_eff->GetXaxis()->SetTitle("#frac{m_{u}}{2m_{e}}");
  graph_eff->GetYaxis()->SetTitle("Efficiency [%]");
  graph_eff->SetMarkerStyle(20);
  graph_eff->SetMarkerSize(1.);
  graph_eff_smear->SetMarkerStyle(20);
  graph_eff_smear->SetMarkerColor(kRed);
  graph_eff->Draw("AP");
  graph_eff_smear->Draw("same P");

  auto legend = new TLegend(0.65,0.2,0.9,0.4);
  legend->AddEntry(graph_eff,"energy without smear","p");
  legend->AddEntry(graph_eff_smear,"energy with smear","p");
  legend->Draw();
  
  c1->SaveAs("../plots/Efficiency.png");

  ////////////////////////////////////////////////////////////////////////
  graph_diff->SetTitle("Difference in efficiencies vs m_{u}");
  graph_diff->GetXaxis()->SetTitle("#frac{m_{u}}{2m_{e}}");
  graph_diff->GetYaxis()->SetTitle("Difference in efficiencies [%]");
  graph_diff->SetMarkerStyle(20);
  graph_diff->SetMarkerSize(1.);
  graph_diff->Draw("AP");
  
  c1->SaveAs("../plots/Difference_efficiency.png");
  
  ////////////////////////////////////////////////////////////////////////
  graph_E->SetTitle("Energy vs m_{u}");
  graph_E->GetXaxis()->SetTitle("#frac{m_{u}}{2m_{e}}");
  graph_E->GetYaxis()->SetTitle("Energy E_{#gamma} [keV]");
  graph_E->SetMarkerStyle(20);
  graph_E->SetMarkerSize(1.);
  graph_E->Draw("AP");
  
  c1->SaveAs("../plots/Energy_vs_mu.png");
  ////////////////////////////////////////////////////////////////////////////
  graph_br->SetTitle("Branching ratio vs m_{u}");
  graph_br->GetXaxis()->SetTitle("#frac{m_{u}}{2m_{e}}");
  graph_br->GetYaxis()->SetTitle("Branching ratio");
  graph_br->SetMarkerStyle(20);
  graph_br->SetMarkerSize(1.);
  graph_br->Draw("AP");

  c1->SaveAs("../plotsBranching_ratio_vs_mu.png");



}
