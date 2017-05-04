/*
Fitting the gain function.
  */

void fit_sipm_wiggle() {
  TRandom3 rd;
  //TFile *f = TFile::Open("wiggle_constE_fac10_g80_m100.root");
  //TFile *f = TFile::Open("wiggle_fac10_g320_m100.root");// without phase and fixed A = 0.4
  //TFile *f = TFile::Open("wiggle_fac10_g80_m100.root");// without phase and fixed A = 0.4
  //TFile *f = TFile::Open("ph_wiggle_fac10_g80_m100.root");//With phase...
  TFile *f = TFile::Open("wiggle_try_80_5.root");// with phase and A from energy 
  //TFile *f = TFile::Open("ph_lowOffset_wiggle_try_320_100.root");
  
  //TH1D *h1=(TH1D*)f->Get("Laser");// As this file has no lasers - just the basic fit
  TH1D *h1=(TH1D*)f->Get("Muons");// As this file has no lasers - just the basic fit
  h1->SetLineColor(kBlue);
  //h1->GetXaxis()->SetRangeUser(0,700);
  /* 
  for(int i=0; i<700; i++){
   if(h1->GetBinContent(i+1)!=0)
     {
       h1->SetBinError(i+1,0.02/sqrt(2000));
       
       
       //h1->SetBinContent((i+1),h1->GetBinContent(i+1)+(0.02/sqrt(20000))*rd.Gaus(0.,1.));
       //h1->SetBinError(i+1,sqrt(0.02*0.02/2000+h1->GetBinError(i+1)*h1->GetBinError(i+1)));
       //Error set as Quadrature of a fixed error due to 2000 pts and 2% rms and error of the plot itself
       
       //h1->SetBinError(i+1,sqrt(0.02/2000));// this works best but all points then have a constant error bar ans so not very convincing...
     }
   else h1->SetBinError(i+1,0.0);
  }
  
  */
  
 //Trying to generalize with a and b as zero and get the plot w/o wiggle 
  
  TF1 *S_N= new TF1("S_N","24.9828*(exp(-(x+1.2)/64.4)-exp(-(x+1.2)/18.))/64.4", 0,700);
  TF1 *S_a= new TF1("S_a","0.000416279*((64.4-18)*exp(-(x+1.2)/64.4)*cos(1.449*(x+1.2))-exp(-(x+1.2)/18.) + 64.4*18*1.439*exp(-(x+1.2)/64.4)*sin(1.449*(x+1.2)))/64.4", 0,700);
  TF1 *S_b= new TF1("S_b","0.000416279*((64.4-18)*exp(-(x+1.2)/64.4)*sin(1.449*(x+1.2))-exp(-(x+1.2)/18.) + 64.4*18*1.439*exp(-(x+1.2)/64.4)*cos(1.449*(x+1.2)))/64.4", 0,700);
  /*

  TF1 *S_N= new TF1("S_N","24.9828*(exp(-x/64.4)-exp(-x/18.))/64.4", 0,700);
  TF1 *S_a= new TF1("S_a","0.000416279*((64.4-18)*exp(-x/64.4)*cos(1.449*x)-exp(-x/18.) + 64.4*18*1.439*exp(-x/64.4)*sin(1.449*x))/64.4", 0,700);
  TF1 *S_b= new TF1("S_b","0.000416279*((64.4-18)*exp(-x/64.4)*sin(1.449*x)-exp(-x/18.) + 64.4*18*1.439*exp(-x/64.4)*cos(1.449*x))/64.4", 0,700);
  */

  TF1 *f2 = new TF1("f2","1-100*[0]*([1]*S_N+[2]*S_a+[3]*S_b)", 0,700);
  //Parameters are avg. Asy, N, a and b
  //f2->SetParameters(2,0.001154,0.0002295,-3.523e-5);
  f2->SetParameters(4.27e-4,1.45,0.139,-0.002633);
  f2->SetNpx(1000);
  
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  h1->SetTitle("");
  h1->GetYaxis()->SetRangeUser(0.987,1.001);
  h1->GetXaxis()->SetRangeUser(0,250);
  h1->GetXaxis()->SetTitle("Time(#mus)");
  h1->GetYaxis()->SetTitle("#frac{G}{G_{0}}");
  h1->SetStats(1);
  h1->Fit("f2","R");
  //h1->Draw();
  //f2->Draw("same");

  TLatex latex;
  latex.SetTextSize(0.037);
  latex.SetTextColor(kRed);
  //latex.DrawLatex(100.5,.993,"#frac{G}{G_{0}} = 1 - #frac{p_{0} n_{0} #tau_{r}}{#tau_{#mu}-#tau_{r}}#left([1+A cos(#omega_a t + #phi)]e^{-t/#tau_{#mu}} - e^{-t/#tau_{r}} #right)");
  //Franco's avg. function...
  //latex.DrawLatex(121.5,.993,"#frac{G}{G_{0}} = 1 - A (N S_{N}(t, #tau_{r}) + a_{f} S_{a}(t, #tau_{r})) + b_{f} S_{b}(t, #tau_{r})");
  gPad->SetGridx();
  gPad->SetGridy();
  
}

