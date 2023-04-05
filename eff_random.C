#include <TRandom3.h>
#include <TH1D.h>
#include <iostream>
#include <vector>
#include <algorithm>

TGraph* graph_A=nullptr;
TGraph* graph_t_shift=nullptr;
TGraph* graph_t_acc=nullptr;
TGraph* graph_A_integral=nullptr;
TGraph* graph_t_shift_integral=nullptr;
TGraph* graph_t_acc_integral=nullptr;


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
  //TH1D *hist_hit = new TH1D(" ", " ", 100, 0, 1.1*t_w);
  //TH1D *hist_prompt = new TH1D(" ", " ", 100, 0, 1.1*t_w);
  for(int i = 0; i < N; i++) {
    prompt = gen.Uniform(0, 1)*t_w;
    hit = gen.Exp(tau) + prompt;
    event[0] = prompt; 
    event[1] = hit; 
    //hist_hit->Fill(hit);
    //hist_prompt->Fill(prompt);
    hits.push_back(event);
  }
  /*TCanvas *c = new TCanvas();
  hist_hit->SetLineColor(kRed);
  hist_hit->Draw();
  hist_prompt->Draw("same");
  c->SaveAs("../plots/prompt_and_hit.png");*/
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
    if(hits[k][1] < time_stop && hits[k][1] > time_start) { N_rand++;}
  }
  N_gamma++;
  }
 return double(N_rand)/double(N_gamma);
}

double integral(TRandom3 gen, double tau, int N, double t_shift, double t_acc){
  double time_start = t_shift;
  double time_stop = t_shift + t_acc;
  double u,t = 0;
  double a = 0;
  double b = 5*t_shift;
  TH1F *hist_prompt = new TH1F(" ", " ", 100, a, b);

  for(int i = 0; i<N; i++){
    u = gen.Uniform(0,1);
    t = -tau*log(exp(-a/tau)-u*(exp(-a/tau)-exp(-b/tau)));
    hist_prompt->Fill(t);
  }
  int start = hist_prompt->FindBin(time_start);
  int stop = hist_prompt->FindBin(time_stop);
  double eff = double(hist_prompt->Integral(start, stop))/double(hist_prompt->Integral());
  delete hist_prompt;
  return eff;
}

void eff_random()
{
  TRandom3 gen;
  gen.SetSeed(0);
  double A = 1e6;
  double t_w = 2e-3;
  double tau = 10e-9;
  int N = 0;
  double prompt, hit = 0;
  double t_shift = 100e-9;
  double t_acc = 50e-9;
  double time_start, time_stop = 0;
  double M = 100;
  double tau_prompt = 4;
  int N_int = int(1e6);

  graph_A = new TGraph(M);
  graph_t_shift = new TGraph(M);
  graph_t_acc = new TGraph(M);
  graph_A_integral = new TGraph(M);
  graph_t_shift_integral = new TGraph(M);
  graph_t_acc_integral = new TGraph(M);
  
  //double i = integral(gen, t_w, N_int, t_shift, t_acc); 

  for(int j = 1; j < M+1; j++){
    A = j*0.1e6;
    N = int(A*t_w);
    cout<<double(j/M*100)<<", "<<N<<endl;
    vector<vector<double>> hits = generator(gen, N, t_w, tau);
    sort(hits.begin(), hits.end());
    //for(int i = 0; i<N; i++) cout<<hits[i][0]*1e7<<", "<<hits[i][1]*1e7<<endl;
    graph_A->SetPoint(j-1, A/1e6, num_random(hits, t_shift, t_acc, N));
    graph_A_integral->SetPoint(j-1, A/1e6, integral(gen, tau_prompt, N_int, t_shift*1e2, t_acc*1e2)*100);
  }
  A = 1e6;
  N = int(A*t_w);
  for(int j = 1; j < M+1; j++){
    t_shift = j*10e-9;
    cout<<double(j/M*100)<<", "<<N<<endl;
    vector<vector<double>> hits = generator(gen, N, t_w, tau);
    sort(hits.begin(), hits.end());
    double num = num_random(hits, t_shift, t_acc, N);
    graph_t_shift->SetPoint(j-1, t_shift/1e-9, num);
    graph_t_shift_integral->SetPoint(j-1, t_shift/1e-9, integral(gen, tau_prompt, N_int, t_shift, t_acc)*100);
  }
  t_shift = 100e-9;
  for(int j = 1; j < M+1; j++){
    t_acc = j*10e-9;
    cout<<double(j/M*100)<<", "<<N<<endl;
    vector<vector<double>> hits = generator(gen, N, t_w, tau);
    sort(hits.begin(), hits.end());
    graph_t_acc->SetPoint(j-1, t_acc/1e-9, num_random(hits, t_shift, t_acc, N));
    graph_t_acc_integral->SetPoint(j-1, t_acc/1e-9, integral(gen, tau_prompt, N_int, t_shift, t_acc)*100);  
  }


  double w = 800;
  double h = 600;
  TCanvas *c = new TCanvas("c", "c", w, h);
  plot(c, graph_A, " ", "A [MBq]", "#frac{N_{rand}}{N_{#gamma}}", "/mnt/home/jmedrala/plots/A_vs_randoms.png");
  plot(c, graph_t_shift, " ", "t_{shift} [ns]", "#frac{N_{rand}}{N_{#gamma}}", "/mnt/home/jmedrala/plots/t_shift_vs_randoms.png");
  plot(c, graph_t_acc, " ", "t_{acc} [ns]", "#frac{N_{rand}}{N_{#gamma}}", "/mnt/home/jmedrala/plots/t_acc_vs_randoms.png");
  plot(c, graph_A_integral, " ", "A [MBq]", "Efficiency [%]", "/mnt/home/jmedrala/plots/A_vs_eff.png");
  plot(c, graph_t_shift_integral, " ", "t_{shift} [ns]", "Efficiency [%]", "/mnt/home/jmedrala/plots/t_shift_vs_eff.png");
  plot(c, graph_t_acc_integral, " ", "t_{acc} [ns]", "Efficiency [%]", "/mnt/home/jmedrala/plots/t_acc_vs_eff.png");
  
}
