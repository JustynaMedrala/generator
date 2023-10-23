/* Purpose:
   1. generation of plots:
	-
        - energy of monochromatic photon vs U boson mass
	- lifetime vs U boson mass
	- branching ratio vs U boson mass
 *
*/
#include <TRandom3.h>
#include <TF1.h>
#include <TMath.h>
#include <string>
#include <iostream>
#include <fstream>

TH1F* hist_Edep =nullptr;
TH1F* hist_Edep_smear =nullptr;
TGraph* graph_eff =nullptr;
TGraph* graph_eff_smear =nullptr;
TGraph* graph_diff =nullptr;
TGraph* graph_E =nullptr;
TGraph* graph_br =nullptr;
TGraph* graph_time =nullptr;

void plot(TCanvas *c, int num, bool leg_opt, TGraph* graph1, const char *label1, TGraph* graph2, const char *label2, const char *title, const char *x_label, const char *y_label, const char *name){
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.13);
  graph1->SetTitle(title);
  graph1->GetXaxis()->SetTitle(x_label);
  graph1->GetXaxis()->SetTitleOffset(1.3);
  graph1->GetYaxis()->SetTitle(y_label);
  graph1->SetMarkerStyle(20);
  graph1->Draw("AP");
  if(num > 1){
    graph2->SetMarkerColor(kRed); 
    graph2->SetMarkerStyle(20);
    graph2->Draw("same P");}
  if(leg_opt == 1){
    auto legend = new TLegend(0.65,0.2,0.9,0.4);
    legend->AddEntry(graph1,label1,"p");
    legend->AddEntry(graph2,label2,"p");
    legend->Draw();
  }
  c->SaveAs(name);
}

void Umass_life(){
  TRandom3 gen;
  gen.SetSeed(0);
  TF1 *f = new TF1("f","(1+[0]*(1-cos(x))+1/(1+[0]*(1-cos(x)))-(sin(x))^2)",0,TMath::Pi());
  double E0, theta, Ek, Edep, eff, eff_smear, sigma, mu = 0;
  double me = 511;
  int bin_min = 30;
  int bin_max = 80;
  int N = 21;
  int M = int(1e6);
  float m[N];
  string plot_path, plot_name, data_path;
  plot_path = "../../plots/";
  data_path = "../data/";

  ofstream branching_ratio, E_gamma, efficiency;
  branching_ratio.open (data_path+"branching_ratio.txt");
  E_gamma.open(data_path+"E_gamma.txt");
  efficiency.open(data_path+"efficiency_vs_mU.txt");

  for(int i=0; i<N-1; i++){
    m[i] = float(i)/float(N-1);
  }

  TCanvas *c = new TCanvas(" ", " ", 700, 600); 
  graph_eff = new TGraph(N);
  graph_eff_smear = new TGraph(N);
  graph_diff = new TGraph(N);
  graph_E = new TGraph(N);
  graph_br = new TGraph(N);
  graph_time = new TGraph(N);

  for(int i = 0; i < N; i++){
    mu = 2*me*m[i];
    E0 = me*(1-mu*mu/me/me/4.);
    f->SetParameter(0, E0/511.);
    hist_Edep = new TH1F("hist_Edep","hist_Edep",100,0,1);
    hist_Edep_smear = new TH1F("hist_Edep_with_smear", "hist_Edep_with_smear",100,0,1);
    for(int j = 0; j < M; j++){
      theta = f->GetRandom();
      Ek = E0/(1+E0/me*(1-cos(theta)));
      Edep = (E0 - Ek)/511.;
      hist_Edep->Fill(Edep);
 
      sigma = 0.044*sqrt(Edep*0.511);//sigma = 0.044*sqrt(Edep [MeV])
      Edep = gen.Gaus(Edep*0.511, sigma)/0.511;

      hist_Edep_smear->Fill(Edep);
    }
    eff = hist_Edep->Integral(bin_min, bin_max)/hist_Edep->Integral();
    eff_smear = hist_Edep_smear->Integral(bin_min, bin_max)/hist_Edep_smear->Integral();
    branching_ratio<<m[i]<<", "<<3.5e-8*(1-m[i]*m[i]*m[i]*m[i])<<endl;
    E_gamma<<m[i]<<", "<<E0<<endl;
    efficiency<<m[i]<<", "<<eff_smear<<endl;
    graph_E->SetPoint(i, m[i], E0);
    graph_eff->SetPoint(i, m[i], eff*100);
    graph_eff_smear->SetPoint(i, m[i], eff_smear*100);
    graph_diff->SetPoint(i, m[i], (eff-eff_smear)*100);
    graph_br->SetPoint(i, m[i], 3.5e-8*(1-m[i]*m[i]*m[i]*m[i]));
    graph_time->SetPoint(i, m[i], 4./(1-m[i]*m[i]*m[i]*m[i]));
    if(i == 0){
      string tit = "E0 =";
      const char *title = (tit+to_string(E0)).c_str();
      hist_Edep_smear->SetTitle(title);
      hist_Edep_smear->Draw();
      plot_name = plot_path+"hist_Edep_max_m.png";
      c->SaveAs(&plot_name[0]);}
    hist_Edep->Delete();
    hist_Edep_smear->Delete();
  }

  branching_ratio.close();
  E_gamma.close();
  efficiency.close();
  //Drawing
  ////////////////////////////////////////////////////////////////////////
  plot_name = plot_path+"Efficiency.pdf";
  plot(c, 2, 1, graph_eff_smear,"energy with smear", graph_eff, "energy without smear", "  ", "#frac{m_{U}}{2m_{e}}", "Efficiency [%]", &plot_name[0]);    

  ////////////////////////////////////////////////////////////////////////
  plot_name = plot_path+"Difference_efficiency.pdf";
  plot(c, 1, 0, graph_diff,nullptr, nullptr, nullptr, "Difference in efficiencies vs m_{U}", "#frac{m_{U}}{2m_{e}}", "Difference in efficiencies [%]", &plot_name[0]);
  
  ////////////////////////////////////////////////////////////////////////
  plot_name = plot_path+"Energy_vs_mu.pdf";
  plot(c, 1, 0, graph_E,nullptr, nullptr, nullptr, "Energy vs m_{U}", "#frac{m_{U}}{2m_{e}}", "Energy E_{#gamma} [keV]", &plot_name[0]);
  
  ////////////////////////////////////////////////////////////////////////////
  plot_name = plot_path+"Branching_ratio_vs_mu.pdf";
  plot(c, 1, 0, graph_br,nullptr, nullptr, nullptr, "Branching ratio vs m_{U}", "#frac{m_{U}}{2m_{e}}", "Branching ratio", &plot_name[0]);  

  ////////////////////////////////////////////////////////////////////////////
  plot_name = plot_path+"Lifetime_vs_mu.pdf";
  plot(c, 1, 0, graph_time,nullptr, nullptr, nullptr, "Lifetime vs m_{U}", "#frac{m_{U}}{2m_{e}}", "Lifetime [s]", &plot_name[0]);    

}
