#include <TRandom3.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>

TH1F* hist_Edep =nullptr;
TH1F* hist_Edep_smear =nullptr;
TH1F* hist_ops_energy =nullptr;
TH1F* hist_E0 =nullptr;
TGraph* graph_eff =nullptr;
TGraph* graph_eff_smear =nullptr;
TGraph* graph_diff =nullptr;
TGraph* graph_E =nullptr;
TGraph* graph_bkg =nullptr;
TGraph* graph_pur =nullptr;
TGraph* graph_pur_vs_min =nullptr;
TGraph* graph_eff_vs_min =nullptr;
TGraph* graph_bkg_vs_min =nullptr;

vector<string> split(const string &str, char sep){
  istringstream ss(str);
  vector<string> v;
  string s;
  while(getline(ss, s, sep))
    v.push_back(s);
  return v;
}

void fill_hist_from_file(const char* filename , TH1* hist_ops_energy)
{
  double energy = 0;
  std::unique_ptr<TFile> hfile( TFile::Open(filename,"READ") );
  TTree *tree = hfile->Get<TTree>("T");
  tree->SetBranchAddress("Energy_att",&energy);

  int nentries = tree->GetEntries();

  for(int i = 0; i < nentries; i++){
    tree->GetEntry(i);
    hist_ops_energy->Fill(energy);
  }
}


void E_dep(int M, double E0, TH1* hist_Edep, TH1* hist_Edep_smear){
    double theta, Ek, Edep, sigma;
    TF1 *f = new TF1("f","(1+[0]*(1-cos(x))+1/(1+[0]*(1-cos(x)))-(sin(x))^2)",0,TMath::Pi());
    f->SetParameter(0, E0);
    TRandom3 gen;
    gen.SetSeed(0);
    for(int j = 0; j < M; j++){
      theta = f->GetRandom();
      Ek = E0/(1+E0*(1-cos(theta)));
      Edep = E0 - Ek;
      hist_Edep->Fill(Edep*511);

      sigma = sqrt(Edep)*0.044*0.511;
      Edep = gen.Gaus(Edep*0.511, sigma)/0.511;

      hist_Edep_smear->Fill(Edep*511);
    }
}

double efficiency(int bin_min, int bin_max, TH1* hist){
  return hist->Integral(bin_min, bin_max)/hist->Integral();
}

void purity(TH1* hist_Edep, TH1* hist_Edep_smear, int N_bin, int bin_min, int bin_max, TGraph* graph_pur_vs_min, TGraph* graph_eff_vs_min, TGraph* graph_bkg_vs_min){
  double eff, bkg, e_min;  

  for(int j = 0; j<=N_bin; j++){
      eff = efficiency(bin_min+j, bin_max, hist_Edep_smear);
      bkg = efficiency(bin_min+j, bin_max, hist_ops_energy);
      //cout<<"Purity: "<<eff/(eff+bkg)<<" bkg: "<<bkg<<" eff: "<<eff<<endl; 
      e_min = hist_Edep_smear->GetBinCenter(j+bin_min)-hist_Edep_smear->GetBinWidth(j+bin_min)/2.;
      graph_pur_vs_min->SetPoint(j, e_min, eff/(eff+bkg)*100);
      graph_eff_vs_min->SetPoint(j, e_min, eff*100);
      graph_bkg_vs_min->SetPoint(j, e_min, bkg*100);
      }                            
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
  c->SaveAs(name);
}


void plot3(TCanvas *c, TGraph* graph, const char *lname_1, TGraph* graph2, const char *lname_2, TGraph* graph3, const char *lname_3, const char *title, const char *x_label, const char *y_label, const char *name){
  gPad->SetLeftMargin(0.13);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitle(x_label);
  graph->GetYaxis()->SetTitle(y_label);
  graph->GetYaxis()->SetRangeUser(0, 101);
  graph->GetXaxis()->SetTitleOffset(1.2);
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kBlue);
  graph2->SetMarkerColor(kRed);
  graph2->SetMarkerStyle(20);
  graph3->SetMarkerStyle(20);
  graph->Draw("AP");
  graph2->Draw("same P");
  graph3->Draw("same P");
  

  auto legend = new TLegend(0.65,0.67,0.98,0.98);
  legend->AddEntry(graph,lname_1,"p");
  legend->AddEntry(graph2,lname_2,"p");
  legend->AddEntry(graph3,lname_3, "p");
  legend->SetTextSize(0.04);
  legend->Draw();

  c->SaveAs(name);
}

void eff_E0(){
  double E0, theta, Ek, Edep, eff, eff_smear, sigma, Energy, bkg = 0;
  int num_bins = 511;
  int bin_min = 0.3*num_bins;
  int bin_max = 0.8*num_bins;
  int N = 21;
  int M = int(1e6);
  float m[N];
  double me = 511;
  double mu = 0;
  float e_min = 0;
  int N_bin = int(bin_max-bin_min+1)/5.11;
  int ind = 10; //0 - mu = 0 keV, 5 - mu = 255.5 keV, 10 - mu = 511 keV, 14 - mu = 715.4 keV
  string plot_path, plot_name, data_path, filename, plot_title;
  plot_path = "../../plots/";
  data_path = "../data/";


  for(int i=0; i<N; i++){
    m[i] = float(i)/float(N-1);
  }


  TFile *hfile = 0;
  filename = data_path+"energies.root";
  hfile = TFile::Open(&filename[0],"READ");
  TTree *tree = hfile->Get<TTree>("T");
  tree->SetBranchAddress("Energy",&Energy);

  TCanvas *c1 = new TCanvas();  
  graph_eff = new TGraph(N);
  graph_eff_smear = new TGraph(N);
  graph_diff = new TGraph(N);
  graph_E = new TGraph(N);
  graph_bkg = new TGraph(N);
  graph_pur = new TGraph(N);
  graph_pur_vs_min = new TGraph(N);
  graph_eff_vs_min = new TGraph(N);
  graph_bkg_vs_min = new TGraph(N);

  hist_ops_energy = new TH1F("hist_ops_energy", "hist_ops_energy", num_bins, 0, 511);
  hist_E0 = new TH1F("hist_energy", "hist_energy", 511, 0, 511);
  hist_ops_energy = (TH1F*)hfile->Get("h_att");
  //fill_hist_from_file("energies.root", hist_ops_energy);
  //fill_hist_from_file("energies.root", "E_0", hist_E0);

  hist_ops_energy->Draw();
  plot_name = plot_path+"hist_ops.png";
  c1->SaveAs(&plot_name[0]);

  for(int i = 0; i < N; i++){
    mu = 2*me*m[i];
    E0 = me*(1-mu*mu/me/me/4.);//sqrt(mu*(4*me-mu))-mu;
    E0 = E0/511.;
    //cout<<E0<<endl;
    hist_Edep = new TH1F("hist_Edep","hist_Edep",num_bins,0,511);
    hist_Edep_smear = new TH1F("hist_Edep_with_smear", "hist_Edep_with_smear",num_bins,0,511);
    E_dep(M, E0, hist_Edep, hist_Edep_smear);
    if(i == ind) purity(hist_Edep, hist_Edep_smear, N_bin, bin_min, bin_max, graph_pur_vs_min, graph_eff_vs_min, graph_bkg_vs_min);
   
    bkg = efficiency(bin_min, bin_max, hist_ops_energy);
    eff = efficiency(bin_min, bin_max, hist_Edep);
    eff_smear = efficiency(bin_min, bin_max, hist_Edep_smear);
    cout<<i<<", "<<E0*511.<<", "<<mu<<", "<<eff_smear<<", "<<bkg<<endl;
    graph_E->SetPoint(i, m[i], E0);
    graph_eff->SetPoint(i, m[i], eff*100);
    graph_eff_smear->SetPoint(i, m[i], eff_smear*100);
    graph_diff->SetPoint(i, m[i], (eff-eff_smear)*100);
    graph_bkg->SetPoint(i, m[i], bkg*100);
    graph_pur->SetPoint(i, m[i], eff_smear/(eff_smear+bkg));
    
    hist_Edep->Delete();
    hist_Edep_smear->Delete();
  }
  //Drawing
  ////////////////////////////////////////////////////////////////////////
  plot_name = plot_path + "Efficiency.pdf";
  plot3(c1, graph_eff, "energy without smear", graph_eff_smear, "energy with smear", graph_bkg, "background (o-Ps)", " ", "#frac{m_{U}}{m_{e}}", "Efficiency [%]", &plot_name[0]);

  ////////////////////////////////////////////////////////////////////////
  plot_name = plot_path + "Difference_efficiency.pdf";
  plot(c1, graph_diff, "Difference in efficiencies vs m_{U}", "#frac{m_{U}}{m_{e}}", "Difference in efficiencies [%]", &plot_name[0]);
  
  ////////////////////////////////////////////////////////////////////////
  plot_name = plot_path + "Energy_vs_mu.pdf";
  plot(c1, graph_E, "Energy vs m_{U}", "#frac{m_{U}}{m_{e}}", "Energy #frac{E_{#gamma}}{E_{e}}", &plot_name[0]);

  ////////////////////////////////////////////////////////////////////////
  plot_name = plot_path + "Purity_vs_mu.pdf";
  plot(c1, graph_pur, " ", "#frac{m_{U}}{m_{e}}", "Purity #frac{N_{sig}}{N_{sig}+N_{bkg}}", &plot_name[0]);
  
  ///////////////////////////////////////////////////////////////////////
  plot_name = plot_path + "Purity_vs_min.pdf";
  stringstream s;
  s<<2*me*m[ind];
  plot_title = "m_{U} = "+s.str()+" [keV]"; 
  plot3(c1, graph_pur_vs_min, "purity #frac{#varepsilon_{sig}}{#varepsilon_{sig}+#varepsilon_{bkg}}", graph_bkg_vs_min, "#varepsilon_{bkg}", graph_eff_vs_min, "#varepsilon_{sig}", &plot_title[0], "E_{min} [keV]", "[%]", &plot_name[0]);
}
