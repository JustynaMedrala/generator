#include <TRandom3.h>
#include <TH1D.h>
#include <iostream>
#include <vector>
#include <algorithm>

TGraph* graph_A=nullptr;
TGraph* graph_t_shift=nullptr;
TGraph* graph_t_acc=nullptr;


void plot(TCanvas *c, TGraph* graph, const char *title, const char *x_label, const char *y_label, const char *name){
  gPad->SetLeftMargin(0.13);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitle(x_label);
  graph->GetYaxis()->SetTitle(y_label);
  graph->SetMarkerStyle(20);
  graph->Draw("AP");
  c->SaveAs(name);
}

vector<vector<double>> generator(TRandom3 gen, int N, double t_w, double tau){
  vector<double> event = {0, 0};
  vector<vector<double>> hits;
  double prompt, hit= 0;
  for(int i = 0; i < N; i++) {
    prompt = gen.Uniform(0, 1)*t_w;
    hit = gen.Exp(tau) + prompt;
    event[0] = prompt; 
    event[1] = hit; 
    hits.push_back(event);
  }
  return hits;
}

double num_random(vector<vector<double>> hits, double t_shift, double t_acc, int N){
  double time_start, time_stop;
  int N_rand = 0;
  int N_gamma = 0;
  for(int i = 0; i < N; i++){
  time_start = hits[i][0]+t_shift;
  time_stop = time_start+t_acc;
  for(int k = i; k < N; k++){
    if(hits[i][1] < time_stop && hits[i][1] > time_start) N_rand++;
    N_gamma++;
  }
  }
 return double(N_rand)/double(N_gamma);
}

void eff_random()
{
  TRandom3 gen;
  gen.SetSeed(0);
  double A = 1e6;
  double t_w = 200e-6;
  double tau = 142e-9;
  int N = 0;
  double prompt, hit = 0;
  double t_shift = 100e-9;
  double t_acc = 5e-9;
  double time_start, time_stop = 0;
  double M = 100;

  graph_A = new TGraph(M);
  graph_t_shift = new TGraph(M);
  graph_t_acc = new TGraph(M);

  for(int j = 1; j < M+1; j++){
    A = j*0.1e6;
    N = int(A*t_w);
    vector<vector<double>> hits = generator(gen, N, t_w, tau);
    sort(hits.begin(), hits.end());
    graph_A->SetPoint(j, A/1e6, num_random(hits, t_shift, t_acc, N));
  }
  A = 1e6;
  N = int(A*t_w);
  for(int j = 1; j < M+1; j++){
    t_shift = j*5e-9;
    vector<vector<double>> hits = generator(gen, N, t_w, tau);
    sort(hits.begin(), hits.end());
    graph_t_shift->SetPoint(j, t_shift/1e-9, num_random(hits, t_shift, t_acc, N));
  }
  t_shift = 100e-9;
  for(int j = 1; j < M+1; j++){
    t_acc = j*1e-9;
    vector<vector<double>> hits = generator(gen, N, t_w, tau);
    sort(hits.begin(), hits.end());
    graph_t_acc->SetPoint(j, t_acc/1e-9, num_random(hits, t_shift, t_acc, N));
  }


  double w = 800;
  double h = 600;
  TCanvas *c = new TCanvas("c", "c", w, h);
  plot(c, graph_A, " ", "A [MBq]", "#frac{N_{rand}}{N_{#gamma}}", "/mnt/home/jmedrala/plots/A_vs_randoms.png");
  plot(c, graph_t_shift, " ", "t_{shift} [ns]", "#frac{N_{rand}}{N_{#gamma}}", "/mnt/home/jmedrala/plots/t_shift_vs_randoms.png");
  plot(c, graph_t_acc, " ", "t_{acc} [ns]", "#frac{N_{rand}}{N_{#gamma}}", "/mnt/home/jmedrala/plots/t_acc_vs_randoms.png");

}
