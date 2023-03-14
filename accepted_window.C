#include <TRandom3.h>
#include <TH1D.h>
#include <iostream>

TGraph* graph_sig =nullptr;

void accepted_window(){
  TRandom3 gen;
  gen.SetSeed(0);
  const int M = 5;
  double t_w[M] = {400, 500, 1000, 1420, 2000};
  int sig_acc[M];
  double sig[M];
  double tau2 = 142e-9;
  double tau1 = 4;
  double t1, u;
  double a = 0;
  double b = 5e-6;
  int N = int(1e6);

  for(int i = 0; i < M; i++){
    sig_acc[i] = 0;
  }

  TH1D *hist_tau = new TH1D(" ", " ", 100, a, 2*b);
  for(int i = 0; i < N; i++){
    u = gen.Uniform(0,1);
    t1 = -tau1*log(exp(-a/tau1)-u*(exp(-a/tau1)-exp(-b/tau1)));
    hist_tau->Fill(t1);
    for(int j = 0; j < M; j++){
      if(t1>tau2&&t1<t_w[j]*1e-9) sig_acc[j]++;
    }
  }

  for(int i = 0; i < M; i++){
    sig[i] = double(sig_acc[i])/double(N)*100;
    cout<<sig[i]<<endl;
  }

  TCanvas* c1 = new TCanvas();

  graph_sig = new TGraph(M, t_w, sig);

  graph_sig->SetTitle("Accepted events in time window");
  graph_sig->GetXaxis()->SetTitle("t_{w} [ns]");
  graph_sig->GetYaxis()->SetTitle("#frac{N_{acc}}{N_{all}} [%]");
  graph_sig->SetMarkerStyle(20);
  graph_sig->SetMarkerSize(1.);
  graph_sig->GetYaxis()->SetRangeUser(-5,105);
  graph_sig->Draw("AP");

  c1->SaveAs("../plots/N_sig.png");

  hist_tau->Draw();

}
