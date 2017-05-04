#include <TMath.h>
#include <TH1.h>
#include <TFile.h>
#include <TF1.h>
#include "TString.h"

Double_t accept(Double_t *x, Double_t *par)
{
  //The acceptance function gets flat after 2200 MeV
  Double_t xx =x[0];
  Double_t f;
  if(xx<2200) f = par[0]*xx/2200.; 
  else f = par[0];
  return f;
}

void Plotter_energy()
{
  //const int num=1;
  TH1F *h2,*ha;
  TH1F *hclone;
  // TCanvas *c;

  TFile *f=new TFile("gm2ringsim_truth.root");
  TFile *f1=new TFile("acceptance.root");
  TString title="truthAnalyzer/track_e_all";
  //  title[1]="truthAnalyzer/gasgun_t";

  TString name="Gas Gun E";
  //name[1]="Positron birth time";
  
  //Acceptance histogram to be rebinned
  ha=(TH1F*)f1->Get("overallAcceptance");
  ha->Rebin(3.226);
  
  //ha->Draw();
  //cout<<" Trying to find "<<title<<endl;
  h2=(TH1F*)f->Get(title);
  h2->SetTitle(title);
  
  TH1D* hEnergy = new TH1D("hEnergy","Positron Energy Deposited",310,0.,3100.);
  TH1D* hEnergy_prob = new TH1D("hEnergy_prob","Prob Positron Energy Deposited",310,0.,3100.);
  for (Int_t i=1;i < 311; i++) 
    hEnergy->SetBinContent(i,ha->GetBinContent(i)*h2->GetBinContent(i));
  double total = hEnergy->Integral();
  for (Int_t i=1;i < 311; i++)
    hEnergy_prob->SetBinContent(i,hEnergy->GetBinContent(i)/total);
  cout<<"Ok!!!\n";
  TFile *f2 = new TFile("prob.root","RECREATE");
  hEnergy_prob->Write("Prob");
  hEnergy->Write("AccEnergy"); 
  //h2->Draw();
  //hEnergy->Draw();
 
  gStyle->SetOptStat(0);
  //c->cd(1);
  //ha->Draw();
  
  //c->cd(2);
  //h2->Draw();   
  
  //c->cd(3);
  //hEnergy->Draw();
  //c->cd(4);
  hEnergy_prob->Draw();
  
  //c->Print("allEnergy.png");
 
  
}
