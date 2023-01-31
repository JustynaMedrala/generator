#include <TRandom3.h>
#include <TH1D.h>
#include <iostream>
#include <vector>
#include <algorithm>

TH1D* hist_prompt =nullptr;
TH1D* hist_hit=nullptr;
void generator()
{
  TRandom3 gen;
  gen.SetSeed(0);
  double A = 1e6;
  double t_w = 20e-6;
  double tau = 140e-9;
  int N = 1000; //int(A*t_w);
  double prompt, hit;
  vector<pair<double, double>> hits;
  hist_prompt= new TH1D("hist_prompt","hist_prompt",100, 0,t_w*1.1);
  hist_hit= new TH1D("hist_hit", "hist_hit", 100, 0, t_w*1.1);
  for (int i = 0; i < N; i++) {
    prompt = gen.Uniform(0, 1)*t_w;
    hit = gen.Exp(tau)+prompt;
    hits.push_back(make_pair(prompt, hit));
    hist_prompt->Fill(prompt);
    hist_hit->Fill(hit);   
  }
  hist_hit->SetLineColor(kRed);
  hist_hit->SetTitle(" ");
  hist_hit->GetXaxis()->SetTitle("t [s]");
  hist_hit->Draw();
  hist_prompt->Draw("same");
  for(int i = 0; i < N; i++) cout<<hits[i].first<<", "<<hits[i].second<<endl;
  cout<<"Po sortowaniu: "<<endl;
  sort(hits.begin(), hits.end());
  for(int i = 0; i < N; i++) cout<<hits[i].first<<", "<<hits[i].second<<endl;
  int num_wrong = 0;
  for(int i = 0; i < N-1; i++){
   if(hits[i].second - hits[i+1].first > 0) num_wrong++;
  }
  cout<<"Hits that apperade after next prompt: "<<num_wrong<<" ("<<num_wrong/double(N)<<")"<<endl;

}
