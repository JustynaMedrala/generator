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
 
  for(int i=0; i<N; i++){
    m[i] = float(i)/float(N-1);
  }


  TFile *hfile = 0;
  hfile = TFile::Open("energies.root","READ");
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
  c1->SaveAs("hist_ops.png");
  for(int i = 0; i < N; i++){
    mu = 2*me*m[i];
    E0 = me*(1-mu*mu/me/me/4.);//sqrt(mu*(4*me-mu))-mu;
    E0 = E0/511.;
    //cout<<E0<<endl;
    hist_Edep = new TH1F("hist_Edep","hist_Edep",num_bins,0,511);
    hist_Edep_smear = new TH1F("hist_Edep_with_smear", "hist_Edep_with_smear",num_bins,0,511);
    E_dep(M, E0, hist_Edep, hist_Edep_smear);
    if(i == 14) purity(hist_Edep, hist_Edep_smear, N_bin, bin_min, bin_max, graph_pur_vs_min, graph_eff_vs_min, graph_bkg_vs_min);
   
    bkg = efficiency(bin_min, bin_max, hist_ops_energy);
    eff = efficiency(bin_min, bin_max, hist_Edep);
    eff_smear = efficiency(bin_min, bin_max, hist_Edep_smear);
    cout<<E0*511.<<", "<<mu<<", "<<eff_smear<<", "<<bkg<<endl;
    graph_E->SetPoint(i, m[i], E0);
    graph_eff->SetPoint(i, m[i], eff*100);
    graph_eff_smear->SetPoint(i, m[i], eff_smear*100);
    graph_diff->SetPoint(i, m[i], (eff-eff_smear)*100);
    graph_bkg->SetPoint(i, m[i], bkg*100);
    graph_pur->SetPoint(i, m[i], eff_smear/(eff_smear+bkg));
    //hist_Edep_smear->Draw();
    //c1->SaveAs("hist_dep.pdf");
    hist_Edep->Delete();
    hist_Edep_smear->Delete();
  }
  //Drawing
  ////////////////////////////////////////////////////////////////////////
  graph_eff->SetTitle(" ");
  graph_eff->GetXaxis()->SetTitle("#frac{m_{U}}{m_{e}}");
  graph_eff->GetYaxis()->SetTitle("Efficiency [%]");
  graph_eff->SetMarkerStyle(20);
  graph_eff->SetMarkerSize(1.);
  graph_eff_smear->SetMarkerStyle(20);
  graph_eff_smear->SetMarkerColor(kRed);
  graph_bkg->SetMarkerStyle(20);
  graph_bkg->SetMarkerColor(kGreen);
  graph_eff->Draw("AP");
  graph_eff_smear->Draw("same P");
  graph_bkg->Draw("same P");

  auto legend = new TLegend(0.65,0.2,0.9,0.4);
  legend->AddEntry(graph_eff,"energy without smear","p");
  legend->AddEntry(graph_eff_smear,"energy with smear","p");
  legend->AddEntry(graph_bkg, "background (o-Ps)", "p");
  legend->Draw();

  
  c1->SaveAs("../plots/Efficiency.pdf");

  ////////////////////////////////////////////////////////////////////////
  graph_diff->SetTitle("Difference in efficiencies vs m_{u}");
  graph_diff->GetXaxis()->SetTitle("#frac{m_{U}}{m_{e}}");
  graph_diff->GetYaxis()->SetTitle("Difference in efficiencies [%]");
  graph_diff->SetMarkerStyle(20);
  graph_diff->SetMarkerSize(1.);
  graph_diff->Draw("AP");
  
  c1->SaveAs("../plots/Difference_efficiency.pdf");
  

  ////////////////////////////////////////////////////////////////////////
  graph_E->SetTitle("Energy vs m_{u}");
  graph_E->GetXaxis()->SetTitle("#frac{m_{U}}{m_{e}}");
  graph_E->GetYaxis()->SetTitle("Energy #frac{E_{#gamma}}{E_{e}}");
  graph_E->SetMarkerStyle(20);
  graph_E->SetMarkerSize(1.);
  graph_E->Draw("AP");
  
  c1->SaveAs("../plots/Energy_vs_mu.pdf");

  ////////////////////////////////////////////////////////////////////////
  gPad->SetLeftMargin(.13);
  graph_pur->SetTitle(" ");
  graph_pur->GetXaxis()->SetTitle("#frac{m_{U}}{m_{e}}");
  graph_pur->GetYaxis()->SetTitle("Purity #frac{N_{sig}}{N_{sig}+N_{bkg}}");
  graph_pur->GetYaxis()->SetTitleOffset(1.3);
  graph_pur->SetMarkerStyle(20);
  graph_pur->SetMarkerColor(kViolet);
  graph_pur->SetMarkerSize(1.);
  graph_pur->Draw("AP");

  c1->SaveAs("../plots/Purity_vs_mu.pdf");
  
  ///////////////////////////////////////////////////////////////////////
  gPad->SetLeftMargin(.13);
  graph_pur_vs_min->SetTitle("m_{U} = 715.4 keV");
  graph_pur_vs_min->GetXaxis()->SetTitle("E_{min} [keV]");
  graph_pur_vs_min->GetYaxis()->SetTitle("[%]");
  graph_pur_vs_min->GetYaxis()->SetTitleOffset(1.3);
  graph_pur_vs_min->GetYaxis()->SetRangeUser(0, 101);
  graph_pur_vs_min->SetMarkerStyle(20);
  graph_bkg_vs_min->SetMarkerStyle(20);
  graph_eff_vs_min->SetMarkerStyle(20);
  graph_pur_vs_min->SetMarkerColor(kViolet);
  graph_bkg_vs_min->SetMarkerColor(kBlack);
  graph_eff_vs_min->SetMarkerColor(kRed);
  graph_pur_vs_min->SetMarkerSize(1.);
  graph_pur_vs_min->Draw("AP");
  graph_eff_vs_min->Draw("same P");
  graph_bkg_vs_min->Draw("same P");

  auto legend1 = new TLegend(0.7,0.7,0.90,0.95);
  legend1->AddEntry(graph_eff_vs_min,"#varepsilon_{sig}","p");
  legend1->AddEntry(graph_bkg_vs_min,"#varepsilon_{bkg}","p");
  legend1->AddEntry(graph_pur_vs_min, "purity #frac{#varepsilon_{sig}}{#varepsilon_{sig}+#varepsilon_{bkg}}", "p");
  legend1->SetTextSize(0.035);
  legend1->Draw();

  c1->SaveAs("../plots/Purity_vs_min.pdf");

}
