#include <TRandom3.h>
#include <TH1D.h>
#include <iostream>

TGraph* graph_sig =nullptr;
TGraph* graph_bkg =nullptr;

void turn_accepted(){
  TRandom3 gen;
  gen.SetSeed(0);
  const int M = 5;
  double t_w[M] = {40, 50, 100, 142, 200};
  int sig_acc[M], bkg_acc[M];
  double sig[M], bkg[M];
  double tau2 = 142e-9;
  double tau1 = 4;
  double t1, t2;
  double a = 0;
  double b = 1e-6;
  double u;
  int N = int(1e6);

  for(int i = 0; i < M; i++){
    sig_acc[i] = 0;
    bkg_acc[i] = 0;
    //cout<<bkg_acc[i]<<endl;
  }


  for(int i = 0; i < N; i++){
    u = gen.Uniform(0,1);
    t1 = -tau1*log(exp(-a/tau1)-u*(exp(-a/tau1)-exp(-b/tau1)));
    t2 = gen.Exp(tau2);
    //cout<<t1<<", "<<(t1>40e-9)<<", "<<t2<<", "<<(t2>40e-9)<<endl;
    for(int j = 0; j < M; j++){
      //cout<<t_w[j]<<endl;
      if(t1>t_w[j]*1e-9) sig_acc[j]++;
      if(t2>t_w[j]*1e-9) bkg_acc[j]++;
    }
  }

  for(int i = 0; i < M; i++){
    sig[i] = double(sig_acc[i])/double(N)*100;
    bkg[i] = double(bkg_acc[i])/double(N)*100;
    cout<<sig[i]<<", "<<bkg[i]<<endl;
  }

  TCanvas* c1 = new TCanvas();

  graph_sig = new TGraph(M, t_w, sig);
  graph_bkg = new TGraph(M, t_w, bkg);

  graph_sig->SetTitle(" ");
  graph_sig->GetXaxis()->SetTitle("t_{acc} [ns]");
  graph_sig->GetYaxis()->SetTitle("#frac{N_{acc}}{N_{all}} [%]");
  graph_sig->SetMarkerStyle(20);
  graph_sig->SetMarkerSize(1.);
  graph_sig->GetYaxis()->SetRangeUser(-5, 105);
  graph_bkg->SetMarkerStyle(20);
  graph_bkg->SetMarkerSize(1.);
  graph_bkg->SetMarkerColor(kRed);
  graph_sig->Draw("AP");
  graph_bkg->Draw("same P");

  auto legend = new TLegend(0.65,0.75,0.9,0.9);
  legend->AddEntry(graph_sig,"signal","p");
  legend->AddEntry(graph_bkg,"background","p");
  legend->Draw();

  c1->SaveAs("../plots/turn_N_acc.png");

}
