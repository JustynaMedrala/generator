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
TGraph* graph_integral_vs_rand_01MBq=nullptr;
TGraph* graph_integral_vs_rand_1MBq=nullptr;
TGraph* graph_integral_vs_rand_10MBq=nullptr;
TGraph* graph_integral_vs_rand_5MBq=nullptr;
TGraph* graph_tau_integral=nullptr;
TGraph* graph_tau=nullptr;


void plot(TCanvas *c, TGraph* graph, const char *title, const char *x_label, const char *y_label, const char *name){
  gPad->SetLeftMargin(0.13);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitle(x_label);
  graph->GetYaxis()->SetTitle(y_label);
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kBlue);
  graph->Draw("AP");
  c->SaveAs(name);
}

vector<vector<double>> generator(TRandom3 gen, double A, double t_w, double tau){
  vector<double> event = {0, 0};
  vector<vector<double>> hits;
  double prompt, hit= 0;
  int N = int(A*t_w);
  
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
    if(hits[k][1] < time_stop && hits[k][1] > time_start) { 
      N_rand++;
      break;
    }
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
  double b = tau;
  int ind = 0;
  int tot = 0;
  TH1F *hist_prompt = new TH1F(" ", " ", 100, a, b);

  for(int i = 0; i<N; i++){
    u = gen.Uniform(0,1);
    t = -tau*log(exp(-a/tau)-u*(exp(-a/tau)-exp(-b/tau)));
    if(t < time_stop && t > time_start) ind++;
    tot++;
    hist_prompt->Fill(t);
  }
  int start = hist_prompt->FindBin(time_start);
  int stop = hist_prompt->FindBin(time_stop);
  double eff = double(hist_prompt->Integral(start, stop))/double(hist_prompt->Integral());
  /*double w = 900;
  double h = 700;
  TCanvas *c = new TCanvas("c", "c", w, h);
  gPad->SetLeftMargin(0.13);
  gStyle->SetOptStat(0);
//  hist_prompt->GetYaxis()->SetRangeUser(5000, 17000);
  hist_prompt->SetLineColor(kViolet);
  hist_prompt->GetYaxis()->SetTitle("number of events");
  hist_prompt->GetXaxis()->SetTitle("t [s]");
  hist_prompt->Draw();
  auto *l_min = new TLine(1,5000,1,15000);
  auto *l_max = new TLine(1.2, 5000,1.2, 15000);
  auto *l_shift = new TArrow(0, 7000, 1, 7000, 0.01, "<>");
  auto *l_acc = new TArrow(1, 7000, 1.2, 7000, 0.01, "<>");
  TLatex *text = new TLatex(0.45, 7200, "#scale[0.55]{t_{shift}}"); 
  TLatex *text1 = new TLatex(1.02, 7200, "#scale[0.55]{t_{acc}}");  
  l_min->SetLineWidth(2);
  l_max->SetLineWidth(2);
  start= hist_prompt->FindBin(1);
  stop= hist_prompt->FindBin(1.2);
  TH1F *hist_part = new TH1F(" ", " ", 100, a, b);
  for(int i = start; i < stop; i++){
    hist_part->AddBinContent(i, hist_prompt->GetBinContent(i));}
  hist_part->SetFillColorAlpha(kBlack, 0.3);
  hist_part->SetLineColor(kBlack);
  hist_part->Draw("same");
  l_min->Draw("same");
  l_max->Draw("same");
  l_shift->Draw();
  l_acc->Draw();
  text->Draw();
  text1->Draw();
  c->SaveAs("../plots/ekspon.pdf");*/
  delete hist_prompt;
  //delete c;
  return double(ind)/double(tot);
}

void eff_random()
{
  TRandom3 gen;
  gen.SetSeed(0);
  double A = 1e6;
  double t_w = 1e-3;
  double tau = 10e-9;
  int N = 0;
  double prompt, hit = 0;
  double t_shift = 200e-9;
  double t_acc = 50e-9;
  double time_start, time_stop = 0;
  double M = 100;
  double tau_prompt = 4;
  int N_int = int(1e9);
  double num_rand = 0;
  double ef_integral = 0;
  double ef_det_ortho = 0.085; 
  double ef_det_para = 1-exp(-0.0959816284917518*2);
  double ef_ortho = 0.423327*ef_det_ortho*pow((1-ef_det_ortho),2);  
  double ef_para = 0.744442*ef_det_para*(1-ef_det_para);
  double norm_rand = ef_ortho + ef_para;

  ofstream t_shift_eff;
  t_shift_eff.open ("t_shift_efficiency.txt");

  graph_A = new TGraph(M);
  graph_t_shift = new TGraph(M);
  graph_t_acc = new TGraph(M);
  graph_A_integral = new TGraph(M);
  graph_t_shift_integral = new TGraph(M);
  graph_t_acc_integral = new TGraph(M);
  graph_integral_vs_rand_01MBq = new TGraph(M);
  graph_integral_vs_rand_1MBq = new TGraph(M);
  graph_integral_vs_rand_10MBq = new TGraph(M);
  graph_integral_vs_rand_5MBq = new TGraph(M);
  graph_tau_integral = new TGraph(M);
  graph_tau = new TGraph(M);

  //double i = integral(gen, t_w, N_int, t_shift, t_acc); 

  /*for(int j = 0; j < M; j++){
    A = 0.1e6+j*0.1e6;
    N = int(A*t_w);
    if(int(j*100/M)%20==0) cout<<"Progress: "<<double(j/M*100)<<"%, N: "<<N<<endl;
    vector<vector<double>> hits = generator(gen, A, t_w, tau);
    sort(hits.begin(), hits.end());
    num_rand = num_random(hits, t_shift, t_acc, N);
    graph_A->SetPoint(j, A/1e6, norm_rand*num_rand);
    ef_integral = integral(gen, tau_prompt, N_int, t_shift, t_acc)*100;
    //cout<<A<<", "<<ef_integral<<endl;
    graph_A_integral->SetPoint(j, A/1e6, ef_integral);
  }
  A = 1e6;
  N = int(A*t_w);
  for(int j = 1; j < M+1; j++){
    t_shift = j*20e-9;
    if(int(j*100/M)%20==0) cout<<"Progress: "<<double(j/M*100)<<"%, N: "<<N<<endl;
    vector<vector<double>> hits = generator(gen, A, t_w, tau);
    sort(hits.begin(), hits.end());
    num_rand = num_random(hits, t_shift, t_acc, N);
    ef_integral = integral(gen, tau_prompt, N_int, t_shift, t_acc)*100;
    graph_t_shift->SetPoint(j-1, t_shift/1e-9, norm_rand*num_rand);
    graph_t_shift_integral->SetPoint(j-1, t_shift/1e-9, ef_integral);
    t_shift_eff<<t_shift<<","<<ef_integral<<endl;
  }
  t_shift = 200e-9;
  for(int j = 1; j < M+1; j++){
    t_acc = j*2e-9;
    if(int(j*100/M)%20==0) cout<<"Progress: "<<double(j/M*100)<<"%, N: "<<N<<endl;
    vector<vector<double>> hits = generator(gen, A, t_w, tau);
    sort(hits.begin(), hits.end());
    num_rand = num_random(hits, t_shift, t_acc, N);
    ef_integral = integral(gen, tau_prompt, N_int, t_shift, t_acc)*100;
    graph_t_acc->SetPoint(j-1, t_acc/1e-9, norm_rand*num_rand);
    graph_t_acc_integral->SetPoint(j-1, t_acc/1e-9, ef_integral);  
  }
  
  t_shift = 200e-9; 

  A = 0.1e6;
  N = int(A*t_w);
  for(int j = 1; j < M+1; j++){
    t_acc = j*2e-9;
    if(int(j*100/M)%20==0) cout<<"Progress: "<<double(j/M*100)<<"%, N: "<<N<<endl;
    vector<vector<double>> hits = generator(gen, A, t_w, tau);
    sort(hits.begin(), hits.end());
    num_rand = num_random(hits, t_shift, t_acc, N);
    ef_integral = integral(gen, tau_prompt, N_int, t_shift, t_acc)*100;
    if(j == 25) cout<<"t_acc: "<<t_acc/1e-9<<", "<<ef_integral<<endl;
    graph_integral_vs_rand_01MBq->SetPoint(j-1, (1-norm_rand*num_rand)*100, ef_integral);
  }

  A = 1e6;
  N = int(A*t_w);
  for(int j = 1; j < M+1; j++){
    t_acc = j*2e-9;
    if(int(j*100/M)%20==0) cout<<"Progress: "<<double(j/M*100)<<"%, N: "<<N<<endl;
    vector<vector<double>> hits = generator(gen, A, t_w, tau);
    sort(hits.begin(), hits.end());
    num_rand = num_random(hits, t_shift, t_acc, N);
    ef_integral = integral(gen, tau_prompt, N_int, t_shift, t_acc)*100;
    if(j == 25) cout<<"t_acc: "<<t_acc/1e-9<<", "<<ef_integral<<endl;    
    graph_integral_vs_rand_1MBq->SetPoint(j-1, (1-norm_rand*num_rand)*100, ef_integral);
  }*/
  t_acc = 50e-9;
  t_shift = 200e-9;
  for(int j = 0; j < M; j++){
    tau_prompt = 4./(1-pow(float(j)/float(M), 4));
    if(int(j*100/M)%20==0) cout<<"Progress: "<<double(j/M*100)<<"%, N: "<<N<<endl;
    vector<vector<double>> hits = generator(gen, A, t_w, tau);
    sort(hits.begin(), hits.end());
    num_rand = num_random(hits, t_shift, t_acc, N);
    ef_integral = integral(gen, tau_prompt, N_int, t_shift, t_acc)*100;
    cout<<tau_prompt<<" <- tau_prompt "<<ef_integral<<endl;
    graph_tau->SetPoint(j-1, tau_prompt, norm_rand*num_rand);
    graph_tau_integral->SetPoint(j-1, tau_prompt, ef_integral);
  }

  /*tau_prompt = 4;

  A = 10e6;
  N = int(A*t_w);
  for(int j = 1; j < M+1; j++){
    t_acc = j*2e-9;
    if(int(j*100/M)%20==0) cout<<"Progress: "<<double(j/M*100)<<"%, N: "<<N<<endl;
    vector<vector<double>> hits = generator(gen, A, t_w, tau);
    sort(hits.begin(), hits.end());
    num_rand = num_random(hits, t_shift, t_acc, N);
    ef_integral = integral(gen, tau_prompt, N_int, t_shift, t_acc)*100;    
    if(j == 25) cout<<"t_acc: "<<t_acc/1e-9<<", "<<ef_integral<<endl;
    graph_integral_vs_rand_10MBq->SetPoint(j-1, (1-norm_rand*num_rand)*100, ef_integral);
  }

  A = 5e6;
  N = int(A*t_w);
  for(int j = 1; j < M+1; j++){
    t_acc = j*2e-9;
    if(int(j*100/M)%20==0) cout<<"Progress: "<<double(j/M*100)<<"%, N: "<<N<<endl;
    vector<vector<double>> hits = generator(gen, A, t_w, tau);
    sort(hits.begin(), hits.end());
    num_rand = num_random(hits, t_shift, t_acc, N);
    ef_integral = integral(gen, tau_prompt, N_int, t_shift, t_acc)*100;
    if(j == 25) cout<<"t_acc: "<<t_acc/1e-9<<", "<<ef_integral<<endl;
    graph_integral_vs_rand_5MBq->SetPoint(j-1, (1-num_rand*norm_rand)*100, ef_integral);
  }*/

  double w = 900;
  double h = 700;
  TCanvas *c = new TCanvas("c", "c", w, h);
  /*plot(c, graph_A, " ", "A [MBq]", "N_{rand} / N_{#gamma*}", "/mnt/home/jmedrala/plots/A_vs_randoms.pdf");
  plot(c, graph_t_shift, " ", "t_{shift} [ns]", "N_{rand} / N_{#gamma*}", "/mnt/home/jmedrala/plots/t_shift_vs_randoms.pdf");
  plot(c, graph_t_acc, " ", "t_{acc} [ns]", "N_{rand} / N_{#gamma*}", "/mnt/home/jmedrala/plots/t_acc_vs_randoms.pdf");
  plot(c, graph_integral_vs_rand_01MBq, "A = 0.1 [MBq]", "Purity [%]", "Efficiency [%]", "/mnt/home/jmedrala/plots/A_01MBq_purity_vs_randoms.pdf");
  plot(c, graph_integral_vs_rand_1MBq, "A = 1 [MBq]", "Purity [%]", "Efficiency [%]", "/mnt/home/jmedrala/plots/A_1MBq_purity_vs_randoms.pdf");
  plot(c, graph_integral_vs_rand_10MBq, "A = 10 [MBq]", "Purity [%]", "Efficiency [%]", "/mnt/home/jmedrala/plots/A_10MBq_purity_vs_randoms.pdf");
  plot(c, graph_integral_vs_rand_5MBq, "A = 5 [MBq]", "Purity [%]", "Efficiency [%]", "/mnt/home/jmedrala/plots/A_5MBq_purity_vs_randoms.pdf");
  plot(c, graph_A_integral, " ", "A [MBq]", "Efficiency [%]", "/mnt/home/jmedrala/plots/A_vs_eff.pdf");
  plot(c, graph_t_shift_integral, " ", "t_{shift} [ns]", "Efficiency [%]", "/mnt/home/jmedrala/plots/t_shift_vs_eff.pdf");
  plot(c, graph_t_acc_integral, " ", "t_{acc} [ns]", "Efficiency [%]", "/mnt/home/jmedrala/plots/t_acc_vs_eff.pdf");
  plot(c, graph_tau, " ", "#tau [s]", "N_{rand} / N_{#gamma*}", "/mnt/home/jmedrala/plots/tau_vs_randoms.pdf");*/
  plot(c, graph_tau_integral, " ", "#tau [s]", "Efficiency [%]", "/mnt/home/jmedrala/plots/tau_vs_eff.pdf");

 
}
