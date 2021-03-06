/*
 I am trying to restructure this code to handle lasers of higher energy and drops
 An instance of gainDropInFill (drop here) calls operator() which is like a root function.
 
 Laser: 
 drop
Function gainSimulation() draws the function for 100 muons of a fill (just one fill)
  */
const int cycle=100;//need 2000 trying with less first...
const int fills=128;//need 16 trying with less first...
const int laser=1;
//const int nMuon=100;
const int nMuon=1;

int factor = 10;
//Use a scale of dropE for 150 MeV positrons
//double dropE = 6.4e-5;//for 150 MeV positrons - to produce Jinst drop for 100 muons
double dropE = factor*4.267e-6;//for 100 MeV positrons - 4.267e-5 works best
double tau_r=18;
int gap = 640; //This is the gap in time in us arising from the frequency of the laser
int offset = 5; //This is the offset we move the laser from fill to fill in us
int bins=700;
class gainDropInFill {
 public:
  gainDropInFill(vector<double> decayTimesIn, vector<double> dropPerDecayIn,
                double recoveryTimeIn)
      : decayTimes(decayTimesIn),
        dropPerDecay(dropPerDecayIn),
        recoveryTime(recoveryTimeIn) {
    //sort(decayTimes.begin(), decayTimes.end());
  }

  double operator()(double* x, double* p) {
    double t = x[0];
    double gain = 1;
    if (t > decayTimes[0]) {
      // drop gain for the first pulse
      gain = 1 - dropPerDecay[0];
      unsigned int i = 1;
      for (; i < decayTimes.size() && decayTimes[i] < t; ++i) {
        // recover gain from pulse i - 1 to pulse i
        gain = 1 -
               (1 - gain) *
                   exp(-(decayTimes[i] - decayTimes[i - 1]) / recoveryTime);
        // drop gain for pulse i
        gain = gain * (1 - dropPerDecay[i]);
      }
      // recover gain from the last pulse before time t to time t
      gain = 1 - (1 - gain) * exp(-(t - decayTimes[i - 1]) / recoveryTime);
      //cout<<"Decay time: "<<decayTimes[i]<<" Drop "<<dropPerDecay[i]<<"Gain "<<gain<<"\n";
    }

    return gain;
  }

 private:
  vector<double> decayTimes;
  // decay drop for one decay (when gain is fully recovered)
  vector<double> dropPerDecay;
  // bias reconvery time constant for one pulse
  double recoveryTime;
};

void energy_dist(int cycle=2000, int fill=16){
  //TFile *f = new TFile::Open("prob.root");
  TFile *f = TFile::Open("prob.root");
  TH1 *Prob = (TH1*)f->Get("Prob");
  TH1D *h = new TH1D("h", "Random energy dist", 3100,0,3100.);
  for(int c=0;c<cycle;c++){
    for(int n=0;n<fill;n++){
      for(int n=0; n<nMuon;n++){
	h->Fill(Prob->GetRandom());
	//cout<<"Energy"<<Prob->GetRandom()<<"\n";
      }
    }
  }
  h->Draw();
  Prob->SetLineColor(kRed);
  Prob->SetLineWidth(2);
  //Prob->Scale(h->Integral()/Prob->Integral());
  Prob->Scale(h->GetMaximum()/(2*Prob->GetMaximum()));
  Prob->Draw("same");
}

void gainSimulationTest(int cycle,int fill, int laser) {
  gRandom->SetSeed(0);
  int no_pts=0;
  double mean,sigma,data,err,data1;
  TGraph* g = new TGraph();
  TCanvas* c1 = new TCanvas();
  TF1* f;
  TF1* f1;
  TFile *file = TFile::Open("prob.root");
  TH1 *Prob = (TH1*)file->Get("Prob");
 
  TH1D* hf; TH1D* hf1; TH1D* hf2; TH1D *hl;
  TH2D* xhist1 = new TH2D("xhist1","SiPM gain Vs. time - Muons only", bins,0,700,50000,0.95,1.0);//final gain for all 2000 cycles
  TH2D* xhist2 = new TH2D("xhist2","SiPM gain Vs. time - Laser and Muons", bins,0,700,50000,0.95,1.0);//final gain for all 2000 cycles
  TH2D* xhist3 = new TH2D("xhist3","SiPM gain Vs. time - Laser only", bins,0,700,50000,0.95,1.0);//final gain for all 2000 cycles
  TH1D* hist1 = new TH1D("hist1","SiPM gain Vs. time - Laser only", bins,0,700);//final gain for all 2000 cycles
  TH1D* hist2 = new TH1D("hist2","SiPM gain Vs. time - All pulses", bins,0,700);//final gain for all 2000 cycles
  TRandom3 *random = new TRandom3;
  int nPulses = nMuon;
  TProfile *py;
  bool laser_on = false;
  
  vector<double> times_Laser;
  vector<double> dropEnergy_Laser;

  int size=0, sizeL=0;
  for(int c=0;c<cycle;c++){
    for(int n=0;n<fill;n++){
      vector<double> decayTimes(cycle*fill*(nMuon+laser),0);
      vector<double> dropPerDecay(cycle*fill*(nMuon+laser),0);
      vector<double> decayTimes_sort;
      vector<double> dropPerDecay_sort;
      vector<double> decayTimes_noLaser;
      vector<double> dropPerDecay_noLaser;
      vector<double> decayTimes_Laser;
      vector<double> dropPerDecay_Laser;
      pair <double,double> apair;
      vector< pair<double,double> > time_E;
      gRandom->SetSeed(0);
      for (int i = 0; i < nPulses; ++i){ 
	double time = gRandom->Exp(64);
	double drop = Prob->GetRandom()*dropE/100.;//took 4.26 ...
	//double drop = dropE;//took 4.26 ...
	decayTimes.push_back(time);
	dropPerDecay.push_back(drop);
	decayTimes_noLaser.push_back(time);
	dropPerDecay_noLaser.push_back(drop);
	//dropPerDecay.push_back(dropE*Prob->GetRandom()/150.);
      }
      for (int i = 0; i < laser; ++i){
	decayTimes.push_back(gap*i+offset*n);
	dropPerDecay.push_back(dropE*20);//2 GeV laser pulses; 
	//dropPerDecay.push_back(1e-2);
	//decayTimes_Laser.push_back(gap*i+offset*n);
	//dropPerDecay_Laser.push_back(dropE*10);
	laser_on=true;
      }      	
      //Here I pair the vectors of times with drop and save them in a vector pair called time_E and then sort this pair wrt times.      
      for(int j=0;j<decayTimes.size();j++)
	{
	  if(decayTimes[j]>0)
	    {
	      apair=make_pair(decayTimes[j],dropPerDecay[j]);
	      time_E.push_back(apair);
	    }
	}
      
      sort(time_E.begin(), time_E.end());

      for (vector<pair<double,double> >::iterator it2 = time_E.begin(); it2 != time_E.end(); ++it2) {
	apair = *it2;
	decayTimes_sort.push_back(apair.first);
	dropPerDecay_sort.push_back(apair.second);
      }

      //sort(decayTimes.begin(), decayTimes.end());
      sort(decayTimes_noLaser.begin(), decayTimes_noLaser.end());
      
      gainDropInFill drop(decayTimes_sort, dropPerDecay_sort, tau_r);//this is the function call to plot SiPM behaviour - Muons + laser....
      //gainDropInFill drop(decayTimes_Laser, dropPerDecay_Laser, tau_r);//this is the function call to plot SiPM behaviour - Muons + laser....
      gainDropInFill drop_noLaser(decayTimes_noLaser, dropPerDecay_noLaser, tau_r);//Muons only..
      //gainDropInFill drop_Laser(decayTimes_Laser, dropPerDecay_Laser, tau_r);//this is the function call to plot SiPM behaviour..

      f = new TF1(Form("%i pulses", nPulses+laser), drop, -10, 700, 0,
		    "gainDropInFill");
      f->SetTitle(Form("%i pulses; time [#mus]; gain", nPulses+laser));
      f1 = new TF1(Form("%i pulses", nPulses), drop_noLaser, -10, 700, 0,
		  "gainDropInFill");
      f1->SetTitle(Form("%i pulses; time [#mus]; gain", nPulses));
      
      //PLot graph only for 1 cycle 1 fill to test
      
      TGraph* g = new TGraph();
      g->SetName(Form("g%i", nPulses));
      for (unsigned int i = 0; i < decayTimes_sort.size(); ++i) {
	g->SetPoint(i, decayTimes_sort[i], f->Eval(decayTimes_sort[i]));
	cout<<"Samples:"<<i<<"="<<decayTimes_sort[i]<<" Gain="<<f->Eval(decayTimes_sort[i])<<" E drop="<<dropPerDecay_sort[i]<<"\n";
      }
      g->SetMarkerStyle(kFullStar);
      g->SetMarkerSize(2);
      g->SetMarkerColor(kBlue);
      f->SetNpx(1000);
      f->Draw();
      g->Draw("p");
      cout<<"Min val:"<<f->GetMinimum()<<"\n";
      for (unsigned int i = 0; i < decayTimes_sort.size(); ++i){//size should be nPulses+laser = Muon + Laser

	xhist2->Fill(decayTimes_sort[i], f->Eval(decayTimes_sort[i]));
	hist2->SetBinContent(hist2->FindBin(decayTimes_sort[i])-1,f1->Eval(decayTimes_sort[i])); 
      }

      for (unsigned int i = 0; i < decayTimes_noLaser.size(); ++i){ //Moun only 
	xhist1->Fill(decayTimes_noLaser[i],f1->Eval(decayTimes_noLaser[i]));
      }
    
      f->SetNpx(1000);
      f1->SetNpx(1000);
      // f2->SetNpx(1000);
      //f->Draw();
            
      //printf("=============End of fill %d===============\n",n+1);
      //c->Print(Form("gainDrop%iPulses.pdf", nPulses));
    }
    printf("=============End of cycle %d===============\n",c+1);    
    
  }
  /*
  hf1 = new TH1D("Muons only","Gain versus time - Muons only",7e2,0,700);
  hf = new TH1D("Laser","Gain versus time - Laser",7e2,0,700);
  hf2 = new TH1D("All","Gain versus time - All pulses ",7e2,0,700);
    
  for ( Int_t i=0; i<bins; i++){   //Muons only
    xhist1->GetXaxis()->SetRange(i,(i+1));
    double time = xhist1->GetMean(1);
    double gain = xhist1->GetMean(2);
    //double gain_err = xhist1->GetRMS(2)/sqrt(xhist1->Integral());
    double gain_err;// = xhist1->GetRMS(2)/sqrt(cycle);
    if(xhist1->Integral()>0)
      gain_err= xhist1->GetRMS(2)/sqrt(xhist1->Integral()/(fill*nMuon));
    else gain_err= 0.0;
    hf1->SetBinContent(hf1->FindBin(time)-1, gain); // Bin starts from 1 and arrays from 0
    hf1->SetBinError(hf1->FindBin(time)-1,gain_err);
    //printf("No. of points %f, RMS %f\n", xhist1->Integral(),xhist1->GetRMS(2));
  }
  
  for ( Int_t i=0; i<bins; i++){   //Muons + Laser  
    xhist2->GetXaxis()->SetRange(i,(i+1));
    double time = xhist2->GetMean(1);
    double gain = xhist2->GetMean(2);
    double gain_err;
    if(xhist2->Integral()>0)
      gain_err= xhist2->GetRMS(2)/sqrt(xhist2->Integral()/(fill*nMuon));
    else gain_err= 0.0;
    //double gain_err = xhist2->GetRMS(2)/sqrt(cycle);
    hf2->SetBinContent(hf2->FindBin(time)-1, gain); // Bin starts from 1 and arrays from 0  
    hf2->SetBinError(hf2->FindBin(time)-1,gain_err);
    //if(xhist2->Integral()>0)
    //printf("Bin %d: No. of points %f, RMS %f, err %f\n", i+1,xhist2->Integral(),xhist2->GetRMS(2),gain_err);
  }  
  
  for ( Int_t i=0; i<(int)(700/offset); i++){   //Laser only     
    xhist2->GetXaxis()->SetRangeUser(i*offset-0.5,i*offset+0.5);
    double time = offset*i; 
    //double gain = xhist2->GetMean(2);                                                          
    //double gain_err = xhist2->GetRMS(2)/sqrt(xhist2->Integral());
    for(int j=0;j<laser;j++)
      if(time>gap*j && time<(gap*j+fill*offset)){
      double gain = xhist2->GetMean(2);
      //double gain_err = xhist2->GetRMS(2)/sqrt(xhist2->Integral());
      double gain_err;
      if(xhist2->Integral()>0)
	gain_err= xhist2->GetRMS(2)/sqrt(xhist2->Integral()/(fill*nMuon));
      else gain_err= 0.0;
      //= xhist2->GetRMS(2)/sqrt(cycle);
      hf->SetBinContent(hf->FindBin(time)-1, gain); 
      hf->SetBinError(hf->FindBin(time)-1,gain_err); 
    } 
  }
  hf1->GetYaxis()->SetRangeUser(0.988,1.002);
  hf1->GetXaxis()->SetRangeUser(0,300);
  hf1->GetYaxis()->SetTitle("#frac{G}{G_{0}}   ");
  hf1->GetXaxis()->SetTitle("Time (#mu s)");
  hf1->SetStats(0);
  gPad->SetGridx();  
  gPad->SetGridy(); 
  hf1->GetYaxis()->SetRangeUser(0.987,1.002);
  hf1->GetXaxis()->SetRangeUser(0,300);
  hf1->Draw("e");
  
  hf->SetMarkerStyle(kFullStar);
  hf->SetMarkerSize(1);
  hf->SetMarkerColor(kRed);
  hf->Draw("same");

  hf2->GetYaxis()->SetRangeUser(0.988,1.002);
  hf2->GetXaxis()->SetRangeUser(0,300);
  hf2->GetYaxis()->SetTitle("#frac{G}{G_{0}}   ");
  hf2->GetXaxis()->SetTitle("Time (#mu s)");
  hf2->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hf2->GetYaxis()->SetRangeUser(0.987,1.002);
  hf2->GetXaxis()->SetRangeUser(0,300);
  hf2->SetLineColor(kGreen);
  hf2->Draw("same,e");
  */
  /*
  TLatex latex;
  latex.SetTextSize(0.037);
  latex.DrawLatex(121.5,.995,Form("Laser %d #mus, %d #mus step",gap,offset));
  //latex.DrawLatex(121.5,.995,"No Laser - 100 mouns each fill");
  latex.DrawLatex(121.5,.994,Form("%d cycles, %d fills each, %d laser shot",cycle,fill,laser));
  latex.DrawLatex(121.5,.992,Form("p_{0} = P_{0} #tau_{#mu} = %0.1f #times 10^{-5}",dropE*1e5));
   
  latex.SetTextColor(kRed);
  latex.DrawLatex(121.5,.993,"#frac{G}{G_{0}} = 1 - #frac{p_{0} n_{0} #tau_{r}}{#tau_{#mu}-#tau_{r}}#left(e^{-t/#tau_{#mu}} - e^{-t/#tau_{r}} #right)");
  TF1 *f3 = new TF1("f3","1-[0]*[1]*64.4*[2]*(exp(-(x/64.4))- exp(-(x/[2])))/(64.4-[2])",0.,160.);
  f3->SetParameters(1e-5,nMuon,tau_r);
  f3->Draw("same");
  
  TFile *file1 = TFile::Open(Form("noWiggle_fac%d_g%d_m%d.root",factor, gap, nMuon), "RECREATE");
  //TFile *file1 = TFile::Open(Form("try_%d_%d.root",gap,offset), "RECREATE");
  hf->Write("Laser");
  hf1->Write("Muons");
  hf2->Write("All");
  */
}

