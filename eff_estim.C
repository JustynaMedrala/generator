#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <TRandom3.h>
#include <TF1.h>
#include <TMath.h>

TGraph* hist_e_mu = nullptr;
TGraph* hist_e_eff = nullptr;
TGraph* hist_N_L = nullptr;
TGraph* hist_sigma = nullptr;
TGraph* hist_BR = nullptr;

vector<string> split(const string &str, char sep){
  istringstream ss(str);
  vector<string> v;
  string s;
  while(getline(ss, s, sep))
    v.push_back(s);
  return v;
}

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


void eff_estim(){
  fstream inFile;
  inFile.open("mu_vs_energy.txt"); 
  string atten;
  double E, a, a_temp;
  int M = 36;
  int ind = 0; 
  double A[4] = {0.1, 1, 5, 10}; 
  double L = 0.64*0.25;
  double T, N_ops, eff;
  int N = 100;
  double T_year = 365*24*60*60;
  double br, m, eff_tw, eff_tw_2 = 0;
  double me = 511;
  double N_o = 200.;
  double BRv = 0;

  m = 2*me*0;
  br = 3.5e-8*(1-m*m);
  eff_tw = 0.060164;
  eff_tw_2 = 0.402847;
  
  TCanvas *c = new TCanvas(" ", " ",  700, 600);

  hist_e_mu = new TGraph(M);
  hist_e_eff = new TGraph(M);


  if(inFile.is_open()){
    while(getline(inFile, atten))
    {
      vector<string> tmp = split(atten, ',');
      E = stod(tmp[0]);
      a_temp = stod(tmp[1]);
      eff = (1-exp(-a_temp*2))*(1 - sin(atan(2)/2.)*sin(atan(2)/2.)*2.);
      if(E > 0.4 && E > 0.6){a = a_temp;}
      hist_e_mu->SetPoint(ind, E, (1-exp(-a_temp*2))*100);    
      hist_e_eff->SetPoint(ind, E, eff*100);
      ind++;                                          
    }
  }
  eff = (1-exp(-a*2))*(1 - sin(atan(2)/2.)*sin(atan(2)/2.)*2.);
  eff = eff*eff_tw*eff_tw_2;


  for(int j = 0; j < 4; j++){
    L = 0.64*0.25*A[j]*1e6;
    hist_N_L = new TGraph(N); 
    hist_BR = new TGraph(N);
    hist_sigma = new TGraph(N); 
    for(int i = 0; i < N; i++){
      T = 0.1*(i+1)*T_year;
      N_ops = L*eff*T*br;
      BRv = N_o/(L*eff*T);
      hist_N_L->SetPoint(i, T/T_year, N_ops);
      hist_BR->SetPoint(i, T/T_year, BRv);
      hist_sigma->SetPoint(i, T/T_year, sqrt(N_ops)*1.959964);
      }
    string file_name = "../plots/N_ops_"+to_string(int(A[j]))+".pdf";
    string title = "A = "+to_string(A[j]).substr(0,4)+" [MBq]";
    plot(c, hist_N_L, &title[0], "T [year]", "N_{o-Ps}", &file_name[0]);
    file_name = "../plots/sigma_ops_"+to_string(int(A[j]))+".pdf";
    plot(c, hist_sigma, &title[0], "T [year]", "#sigma_{o-Ps}", &file_name[0]);
    file_name = "../plots/BR_ops_"+to_string(int(A[j]))+".pdf";
    gPad->SetLogy();
    plot(c, hist_BR, &title[0], "T [year]", "Branching ratio", &file_name[0]);
    delete hist_N_L; 
  }
  gPad->SetLogx();
  plot(c, hist_e_eff, "Attenuation in water and geometrical efficiency", "E [MeV]", "Efficiency [%]", "../plots/Attenuation_eff.pdf");
  gPad->SetLogy();
  plot(c, hist_e_mu, "Attenuation in water", "E [MeV]", "Efficiency [%]", "../plots/Attenuation.pdf");

  return 0;
}
