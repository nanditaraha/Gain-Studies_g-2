/*
Fitting the gain function.
  */

void fit_sipm() {
  TRandom3 rd;
  //TFile *f = TFile::Open("noWiggle_constE_fac10_g80_m100.root");
  TFile *f = TFile::Open("noWiggle_fac10_g320_m100.root");
  TH1D *h1=(TH1D*)f->Get("Muons");// As this file has no lasers - just the basic fit
  //h1->SetName("Laser");
  h1->SetLineColor(kBlue);
  //h1->GetXaxis()->SetRangeUser(0,700);
  /* 
  for(int i=0; i<700; i++){
    if(h1->GetBinContent(i+1)!=0){
	h1->SetBinError(i+1,0.02/sqrt(2000));
	h1->SetBinContent((i+1),h1->GetBinContent(i+1)+(0.02/sqrt(2000))*rd.Gaus(0.,1.));
	//h1->SetBinError(i+1,sqrt(0.02*0.02/2000+h1->GetBinError(i+1)*h1->GetBinError(i+1)));
	//Error set as Quadrature of a fixed error due to 2000 pts and 2% rms and error of the plot itself
	  
	//h1->SetBinError(i+1,sqrt(0.02/2000));// this works best but all points then have a constant error bar ans so not very convincing...
    }
    else h1->SetBinError(i+1,0.0);
  }
  */
 

  
  //These constants work great..
  TF1 *f3 = new TF1("f3","1-100*[0]*[1]*18*(exp(-((x+1.2)/64.4))- exp(-((x+1.2)/18)))/(64.4-18)",0.,400);
  //f3->SetParNames("P_{0}","n_{0}","#tau_{r}");
  f3->SetParameters(4.267e-4,1.45);//1.47 is the average energy in GeV. 4.27e-4 is the energy drop per GeV used in my simulaiton
  //f3->FixParameter(0,4.3e-5);
  //f3->SetParameters(1.052e-5,100,18);
  //f3->SetLineColor(kBlue);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  h1->SetTitle("");
  h1->GetYaxis()->SetRangeUser(0.987,1.001);
  h1->GetXaxis()->SetRangeUser(0,250);
  h1->GetXaxis()->SetTitle("Time(#mus)");
  h1->GetYaxis()->SetTitle("#frac{G}{G_{0}}");
  h1->SetStats(1);
  h1->Fit("f3","R");
  //h1->Draw();
  //f3->Draw("same");
  
  TLatex latex;
  latex.SetTextSize(0.037);
  latex.SetTextColor(kRed);
  latex.DrawLatex(121.5,.993,"#int_{a}^{b} x^{2} dx");
  gPad->SetGridx();
  gPad->SetGridy();
}

