/*Purpose: 
  1. production of plots:
     - attenuation multiplied by geometrical efficiency vs photon energy (graph and colz)
     - histogram of efficiencies
*/

#include <iostream>
#include <sstream>
#include <fstream>

TH2F* hist_e_eff = nullptr;
TH1F* hist_eff = nullptr;
TH1F* hist_E0 =nullptr;
TGraph* graph_eff =nullptr;

vector<string> split(const string &str, char sep){
  //spliting string to vector of strings separated by "sep"
  istringstream ss(str);
  vector<string> v;
  string s;
  while(getline(ss, s, sep))
    v.push_back(s);
  return v;
}

void fill_hist_from_file(const char* filename , const char* branch,  TH1* hist_ops_energy)
{
  //filling histogram hist_ops_energy from specific branch of file
  double energy = 0;
  std::unique_ptr<TFile> hfile( TFile::Open(filename,"READ") );
  TTree *tree = hfile->Get<TTree>("T");
  tree->SetBranchAddress(branch,&energy);

  int nentries = tree->GetEntries();

  for(int i = 0; i < nentries; i++){
    tree->GetEntry(i);
    hist_ops_energy->Fill(energy);
  }
}

void attenuation(){
  double E[39999], a[39999], eff[39999], att[39999];
  int ind = 0;
  int num_bins = 511;
  int canvas_w = 700;
  int canvas_h = 600;

  string atten, plot_path, plot_name;
  fstream inFile;
  ofstream outFile;
  inFile.open("../data/attenuation_more.txt");
  outFile.open("../data/efficiency_vs_energy.txt");
  plot_path = "../../plots/";

  hist_E0 = new TH1F("hist_energy", "hist_energy", num_bins, 0, 511);
  fill_hist_from_file("../data/energies.root", "E_0", hist_E0);
  hist_eff = new TH1F("hist_eff", "hist_eff", 30, 5, 50);
  hist_e_eff = new TH2F("hist_E_eff", "hist_E_eff", num_bins, 0, num_bins, 35, 5, 50);
 
  
  if(inFile.is_open()){
    while(getline(inFile, atten))
    {
      vector<string> tmp = split(atten, ',');
      E[ind] = stod(tmp[0])*1000.;
      a[ind] = stod(tmp[1]);
      //efficiency from attenuation
      att[ind] = (1-exp(-a[ind]*2))*100;
      //efficiency from attenuation multiplied by geometrical efficiency
      eff[ind] = (1-exp(-a[ind]*2))*(1 - sin(atan(1)/2.)*sin(atan(1)/2.)*2.)*100;
      outFile<<E[ind]<<","<<eff[ind]<<endl;
      ind++;                              
    }
  }
  for(int i = 0; i< num_bins; i++){
     hist_eff->Fill(eff[2*i+1], hist_E0->GetBinContent(i));
     hist_e_eff->Fill(hist_E0->GetBinCenter(i), eff[2*i+1], hist_E0->GetBinContent(i));
  }
  
  //DRAWING and SAVING:
  TCanvas* c = new TCanvas("canvas", "canvas", canvas_w, canvas_h);

  gStyle->SetOptStat(0);

  hist_e_eff->SetTitle(" ");
  hist_e_eff->GetXaxis()->SetTitle("E_{#gamma} [keV]");
  hist_e_eff->GetYaxis()->SetTitle("Efficiency [%]");
  hist_e_eff->Draw("colz");
  plot_name = plot_path+"Eff_vs_energy.pdf";
  c->SaveAs(&plot_name[0]);

  gPad->SetLogy();
  hist_eff->SetTitle(" ");
  hist_eff->GetXaxis()->SetTitle("Efficiency [%]");
  hist_eff->SetMarkerStyle(20);
  hist_eff->SetLineColor(kBlue);
  hist_eff->SetFillColor(kBlue);
  hist_eff->Draw("hist");
  plot_name = plot_path+"Efficiency_histogram.pdf";
  c->SaveAs(&plot_name[0]);

  gPad->SetLogx();
  graph_eff = new TGraph(ind, E, eff);
  graph_eff->GetXaxis()->SetTitle("E_{#gamma} [keV]");
  graph_eff->GetXaxis()->SetTitleOffset(1.2);
  graph_eff->GetYaxis()->SetTitle("Efficiency [%]");
  graph_eff->SetTitle(" ");
  graph_eff->SetLineWidth(2.);
  graph_eff->SetLineColor(kRed);
  graph_eff->Draw("AL");
  plot_name = plot_path+"Efficiency_geometrical_and_attenuation.pdf";
  c->SaveAs(&plot_name[0]);
}
