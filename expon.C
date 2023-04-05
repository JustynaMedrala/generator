#include <TRandom3.h>
#include <TH1D.h>
#include <iostream>
#include <vector>
#include <algorithm>

TH1D* hist_prompt =nullptr;
TH1D* hist_hit=nullptr;
void expon()
{
  TRandom3 gen;
  gen.SetSeed(0);
  double A = 1e6;
  double t_w = 20e-6;
  double tau = 142e-9;
  int N = int(1e6);
  double prompt, hit;
  vector<pair<double, double>> hits;
  hist_prompt= new TH1D("hist_prompt","hist_prompt",200, 0,t_w*1.1);
  hist_hit= new TH1D("hist_hit", "hist_hit", 200, 0, t_w*1.1);
  for(int j = 0; j < 5; j++){
    prompt = gen.Uniform(0, 1)*t_w;
    hist_prompt->Fill(prompt);
    for (int i = 0; i < N; i++) {
      //cout<<i<<endl;
      hit = gen.Exp(tau)+prompt;
      hist_hit->Fill(hit);
      //cout<<hit<<endl;
    }
  }
  TCanvas *c1 = new TCanvas();
  gStyle->SetOptStat(0);
  hist_hit->SetFillColorAlpha(kRed, 0.3);
  hist_hit->SetLineColor(kRed);
  hist_hit->SetTitle(" ");
  hist_hit->Scale(1./hist_hit->GetEntries());
  hist_prompt->Scale(1./hist_prompt->GetEntries());
  hist_hit->GetXaxis()->SetTitle("t [s]");
  hist_hit->GetYaxis()->SetRangeUser(0., 0.3);
  //hist_prompt->Draw("hist");
  hist_hit->Draw("hist");
  hist_prompt->Draw("hist same");

  auto legend = new TLegend(0.7,0.8,0.9,0.9);
  legend->AddEntry(hist_prompt,"prompts","l");
  legend->AddEntry(hist_hit,"hits","l");
  legend->Draw();
  c1->SaveAs("../plots/generated.pdf");
}  
