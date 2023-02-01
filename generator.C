#include <TRandom3.h>
#include <TH1D.h>
#include <iostream>
#include <vector>
#include <algorithm>

TH1D* hist_prompt =nullptr;
TGraph* hist_A=nullptr;
void generator()
{
  TRandom3 gen;
  gen.SetSeed(0);
  double A = 0.1;
  double t_w = 20e-6;
  double tau = 140e-9;
  int N = int(A*t_w);
  double prompt, hit;
  int num_wrong = 0;
  vector<pair<double, double>> hits;
  hist_prompt= new TH1D("hist_prompt","hist_prompt",100, 0,t_w*1.1);
  hist_A= new TGraph();
  TCanvas *c1 = new TCanvas();
  for(int j = 0; j < 10; j++){
    N = int(A*t_w*1e6);
    for (int i = 0; i < N; i++) {
      prompt = gen.Uniform(0, 1)*t_w;
      hit = gen.Exp(tau)+prompt;
      hits.push_back(make_pair(prompt, hit));
    }
    sort(hits.begin(), hits.end());
    num_wrong = 0;
    for(int i = 0; i < N-1; i++){
      if(hits[i].second - hits[i+1].first > 0) num_wrong++;
      hist_prompt->Fill(hits[i+1].first - hits[i].first);
    }
    hist_prompt->SetTitle(" ");
    hist_prompt->GetXaxis()->SetTitle("#Deltat [s]");
    hist_prompt->Draw();
    if(A == 1) c1->SaveAs("delta_t.png");
    cout<<"A: "<<A<<" Hits that appeared after next prompt: "<<num_wrong<<" ("<<num_wrong/double(N)<<")"<<endl;
    hist_A->SetPoint(j, A, num_wrong);
    A = int(A+1);
  }    
  hist_A->SetTitle(" ");
  hist_A->GetYaxis()->SetTitle("wrong order");
  hist_A->GetXaxis()->SetTitle("A [MBq]");
  hist_A->SetMarkerStyle(20);
  //hist_A->SetMarkerSize(2);
  hist_A->Draw("AP");
  c1->SaveAs("aktywnosc.png");
}
