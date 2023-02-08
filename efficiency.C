#include <TRandom3.h>
#include <TH1D.h>
#include <iostream>
#include <vector>
#include <algorithm>

TH1D* hist_prompt =nullptr;
TGraph* hist_A=nullptr;

void efficiency()
{
  TRandom3 gen;
  gen.SetSeed(0);
  double A = 1e6;
  double t_w = 1;
  double tau = 140e-9;
  int N = int(A*t_w);
  double prompt, hit;
  double t_acc = 200e-9;
  int num_1, num_2 = 0;
  int num_wrong, num_prop = 0;
  double U1, U2 = 0;
  double ef_prompt = 5*0.1;
  double ef_hit = 5*0.1*0.1*0.1;
  vector<double> event = {0, 0, 0, 0};
  vector<vector<double>> hits;
  for(int i = 0; i < N; i++) {
    prompt = gen.Uniform(0, 1)*t_w;
    U1 = gen.Uniform(0, 1);
    U2 = gen.Uniform(0, 1);
    hit = gen.Exp(tau) + prompt;
    event[0] = prompt; 
    event[1] = hit; 
    event[2] = double(U1 <= ef_prompt);
    event[3] = double(U2 <= ef_hit);
    hits.push_back(event);
  }
  sort(hits.begin(), hits.end());
  for(int i = 0; i < N-1; i++){
    if(hits[i][1] - hits[i+1][0] > 0) num_wrong++;
  }
  for(int i = 0; i < N; i++){
    //cout<<hits[i][0]<<", "<<hits[i][0]<<", "<<hits[i][2]<<", "<<hits[i][3]<<endl;
    if(hits[i][1] - hits[i][0] < t_acc && hits[i][2] > 0.5 && hits[i][3] > 0.5) num_prop++;
    if(hits[i][1] - hits[i][0] > t_acc && hits[i][2] > 0.5 && hits[i][3] > 0.5) num_1++;
    if(hits[i][2] < 0.5){
      for(int wrong_i = i; wrong_i < N; wrong_i++){
        if(hits[wrong_i][2] > 0.5){
          if(hits[i][1] - hits[wrong_i][0] < t_acc && hits[wrong_i][2] > 0.5 && hits[i][3] > 0.5) num_2++;
          break;
        }
      }
    }
  }
  cout<<"N: "<<N<<", True (delta(t) < tw): "<<num_prop<<", delta(t) > tw: "<<num_1<<", missing prompt and delta(t) < tw: "<<num_2<<endl;
  cout<<"Wrong order: "<<num_wrong<<endl;
}
