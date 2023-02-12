#include <TRandom3.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>

TH1F* hist_Edep =nullptr;
TH1F* hist_Ek =nullptr;

void energy(){
  TRandom3 gen;
  gen.SetSeed(0);
  TF1 *f = new TF1("f","1/(1+(1-cos(x)))^2*(1+(1-cos(x))+1/(1+(1-cos(x)))-(sin(x))^2)",0,TMath::Pi()); 
  double theta = f->GetRandom(); 
  cout<<theta<<endl;
  double E0 = 1;
  double Ek = 0;
  int N = int(1e6);

  hist_Edep = new TH1F("hist_Edep","hist_Edep",100,-1,1);
  hist_Ek = new TH1F("hist_Ek","hist_Ek",100, 0,1);

  for(int i = 0; i<N; i++){
    theta = f->GetRandom();
    //cout<<theta<<endl;
    Ek = E0/(1+E0*(1-cos(theta)));
    //cout<<Ek-E0<<endl;
    hist_Edep->Fill(E0-Ek);
    hist_Ek->Fill(Ek);
  }
  TCanvas *c = new TCanvas();
  ////////////////////////////////////
  f->SetTitle("  ");
  f->GetXaxis()->SetTitle("#theta");
  f->GetYaxis()->SetTitle("f(#theta)");
  f->Draw();
  c->SaveAs("f(theta).png");
  ///////////////////////////////////
  hist_Edep->SetTitle("   ");
  hist_Edep->GetXaxis()->SetTitle("E_{dep}");
  hist_Edep->Draw();
  c->SaveAs("hist_Edep.png");
  //////////////////////////////////
  hist_Ek->SetTitle("   ");
  hist_Ek->GetXaxis()->SetTitle("E_{#gamma'}");
  hist_Ek->Draw();
  c->SaveAs("hist_Ek.png");
}
