/*Purpose:
  1. generating plots:
     - number of o-ps -> gamma U decay vs time of observation
*/

#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <TRandom3.h>
#include <TF1.h>
#include <TMath.h>
#include <sstream>

TGraph* hist_e_mu = nullptr;
TGraph* hist_e_eff = nullptr; 
TGraph* hist_e_mu_part = nullptr;
TGraph* hist_e_eff_part = nullptr; 
TGraphErrors* hist_N_L = nullptr;
TGraph* hist_N_L_low = nullptr;
TGraph* hist_N_L_high = nullptr;
TGraph* hist_sigma = nullptr;
TGraph* hist_BR = nullptr;

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 1)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return std::move(out).str();
}


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
  gPad->SetBottomMargin(0.13);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitle(x_label);
  graph->GetYaxis()->SetTitle(y_label);
  graph->GetXaxis()->SetTitleOffset(1.2);
  graph->GetYaxis()->SetTitleOffset(1.3);
  graph->GetXaxis()->SetTitleSize(0.045);
  graph->GetYaxis()->SetTitleSize(0.045); 
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kBlue);
  graph->Draw("AP");

  TLatex *text = new TLatex();
  text->SetTextSize(0.05);
  text->SetTextAlign(13);
  text->SetTextColorAlpha(kBlue, 0.35);
  text->DrawLatexNDC(0.2, 0.8, "PRELIMINARY");


  TLatex *text1 = new TLatex();
  text1->SetTextSize(0.05);
  text1->SetTextAlign(13);
  text1->SetTextColorAlpha(kRed, 0.35);
  text1->DrawLatexNDC(0.5, 0.2, "toy Monte Carlo");
 
  c->SaveAs(name);
}

void plot2(TCanvas *c, TGraph* graph,TGraph* graph2, const char *title, const char *x_label, const char *y_label, const char *name){
  gPad->SetLeftMargin(0.13);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitle(x_label);
  graph->GetYaxis()->SetTitle(y_label);
  graph->GetXaxis()->SetTitleOffset(1.2);
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kBlue);
  graph2->SetMarkerColor(kRed);
  graph2->SetMarkerStyle(20);
  graph->Draw("AP");
  graph2->Draw("same P");
  c->SaveAs(name);
}

void plot3(TCanvas *c, TGraph* graph,TGraph* graph2, TGraph* graph3, const char *title, const char *x_label, const char *y_label, const char *name){
  gPad->SetLeftMargin(0.13);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitle(x_label);
  graph->GetYaxis()->SetTitle(y_label);
  graph->GetXaxis()->SetTitleOffset(1.2);
  //graph->SetMarkerStyle(20);
  graph->SetLineColor(kBlue);
  graph->SetLineWidth(2);
  //graph2->SetMarkerColor(kRed);
  graph2->SetLineStyle(9);
  graph3->SetLineStyle(9);
  graph2->SetLineWidth(2);
  graph3->SetLineWidth(2);  
  graph->Draw("APL");
  graph2->Draw("same PL");
  graph3->Draw("same PL");
  c->SaveAs(name);
}


void eff_estim(){
  int M = 36;
  int ind = 0; 
  double A[4] = {0.1, 1, 5, 10}; 
  double L = 0.36*0.75;
  double T, N_ops, eff;
  int N = 50;
  double T_year = 365*24*60*60;
  double br, m, E0, eff_det = 0;
  double me = 511;
  double N_o = 200.;
  double BRv = 0;
  double m_u[4] = {0, 0.25, 0.5, 0.7};
  double eff_tw = 0.0124222;
  double eff_tw_2[4] = {0.745171, 0.737113, 0.611849, 0.00770294};
  double a[4] = {0.0959816284917518, 0.09861948616546291, 0.10794745523227846, 0.12498662494925405};
  string plot_path, plot_name, plot_title;
  plot_path="../../plots/";
  //E0: 511, 479.062, 383.25, 260.61
  
  for(int k = 0; k < 4; k++){
  m = 2*me*m_u[k];
  E0 = me*(1-m*m/me/me/4);
  br = 3.5e-8*(1-pow(m_u[k],4));

  TCanvas *c = new TCanvas(" ", " ",  700, 600);

  eff = (1-exp(-a[k]*2))*(1 - sin(atan(1)/2.)*sin(atan(1)/2.)*2.); 
  eff = eff_tw*eff*eff_tw_2[k];

  //cout<<"Average: "<<avg_att/double(ind_att)<<endl;
  //cout<<"Eff detekcji: "<<1-exp(-a*2)<<", "<<a<<endl;
  //cout<<"Eff geom: "<<(1 - sin(atan(2)/2.)*sin(atan(2)/2.)*2.)<<endl;
  for(int j = 0; j < 4; j++){
    L = 0.36*0.75*A[j]*1e6;
    hist_N_L = new TGraphErrors(N); 
    hist_N_L_low = new TGraph(N);      
    hist_N_L_high = new TGraph(N);
    hist_BR = new TGraph(N);
    hist_sigma = new TGraph(N); 
    for(int i = 0; i < N; i++){
      T = 0.2*(i+1)*T_year;
      N_ops = L*eff*T*br;
      if(T==3*T_year) N_o = N_ops;
      hist_N_L->SetPoint(i, T/T_year, N_ops);
      hist_N_L->SetPointError(i, 0, sqrt(N_ops));
      hist_N_L_low->SetPoint(i, T/T_year, N_ops-sqrt(N_ops));
      hist_N_L_high->SetPoint(i, T/T_year, N_ops+sqrt(N_ops));
      hist_sigma->SetPoint(i, T/T_year, sqrt(N_ops)/double(N_ops));
      }
    cout<<"N_obs: "<<N_o<<", aktywność: "<<A[j]<<endl;
    for(int i = 0; i < N; i++){
      T = 0.2*(i+1)*T_year;
      BRv = N_o/(L*eff*T);
      hist_BR->SetPoint(i, T/T_year, BRv);
      }
    gPad->SetLogy(0);
    plot_name = plot_path+"N_ops_"+to_string_with_precision<float>(A[j])+"_"+to_string_with_precision<float>(m)+".pdf";
    plot_title = "m_{U} = "+to_string_with_precision<float>(m)+" [keV], A = "+to_string_with_precision<float>(A[j])+" [MBq]";
    plot(c, hist_N_L, &plot_title[0], "T [year]", "N_{o-Ps}", &plot_name[0]);
    plot_name = plot_path+"sigma_ops_"+to_string_with_precision<float>(A[j])+"_"+to_string_with_precision<float>(m)+".pdf";
    plot(c, hist_sigma, &plot_title[0], "T [year]", "#sqrt{N_{o-Ps}}/N_{o-Ps}", &plot_name[0]);
    plot_name = plot_path+"BR_ops_"+to_string_with_precision<float>(A[j])+"_"+to_string_with_precision<float>(m)+".pdf";
    gPad->SetLogy();
    plot(c, hist_BR, &plot_title[0], "T [year]", "Branching ratio", &plot_name[0]);
    delete hist_N_L; 
    delete hist_BR;
    delete hist_sigma;
  }
  }
  return 0;
}
