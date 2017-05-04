void my_accept()
{

  int bin =1000;
  //int bin =310;
  //PHASE def...
  //Convert mrad to rad....
  double p0 =    -0.255134/1000.;
  double p1 =      65.3034/1000.;
  double p2 =     -705.492/1000.;
  double p3 =      5267.21/1000.;
  double p4 =     -23986.5/1000.;
  double p5 =      68348.1/1000.;
  double p6 =      -121761/1000.;
  double p7 =       131393/1000.;
  double p8 =       -78343/1000.;
  double p9 =      19774.1/1000.;

  TF1* phase = new TF1("phase","([0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 + [5]*x^5 + [6]*x^6 + [7]*x^7 + [8]*x^8 + [9]*x^9)",0.,1.);
  phase->SetParNames("p0","p1","p2","p3","p4","p5","p6","p7","p8","p9");
  phase->SetParameters(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9);
  phase->SetNpx(10000000);
  phase->SetLineColor(2);
  //phase->Draw();
  //TH1D* hp = new TH1D("hp","Phi(y)",bin,0.,1.);
  TH1D* hp = new TH1D("hp","Phi(y)",bin,0.,1.);
  hp->Eval(phase);
  double para[10]={p0,p1,p2,p3,p4,p5,p6,p7,p8,p9};
  double average = phase->Mean(0,1,para,0.00001);

  //ASYMMETRY...
  TF1* Asy = new TF1("Asy","(-8.*x^2+x+1.)/(4.*x^2-5.*x-5.)",0.,1.);
  Asy->SetNpx(10000000);
  TH1D* ha = new TH1D("ha","A(y)",bin,0.,1.);
  ha->Eval(Asy);
  //ha->Draw();

  //Nfunction...
  //TF1* N = new TF1("N","(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0.,1.);
  TF1* N = new TF1("N","(x-1.)*(4.*x*x-5*x-5.)",0.,1.);
 
  TF1* N_Asy_p_sin = new TF1("Asy_p_s","-Asy*sin(phase)*x",0.,1.);
  N_Asy_p_sin->SetNpx(10000000);
  //N_Asy_p_sin->Draw();
  //double average = N_Asy_p_sin->Mean(0,1,para,0.00001);
  //printf("%f\n",average);
  TH1D* hb_1 = new TH1D("hb_1","a1(y)",bin,0.,1.);
  hb_1->Eval(N_Asy_p_sin);

  TF1* N_Asy_p_cos = new TF1("Asy_p_c","N*Asy*cos(phase)*x",0.,1.);
  //N_Asy_p_cos->SetNpx(10000000);
  //N_Asy_p_cos->Draw();
  //double average1 = N_Asy_p_cos->Mean(0,1,para,0.00001);
  
  
  TF1 *Acc = new TF1("Acc","(-0.01334+2.86*x-10.94*x^2+24.106*x^3-18.4*x^4-6.57*x^5+16.58*x^6-6.87*x^7)*(1-exp((x-0.998)/2.55625e-03))", 0,1);

  TF1 *f2 = new TF1("f2","N*Acc",0,1);
  TF1 *f3 = new TF1("f3","N*Asy*cos(phase)*Acc",0,1);
  TF1 *f4 = new TF1("f4","-N*Asy*sin(phase)*Acc",0,1);

  printf("Integral = %f \n",f4->Integral(0,1,0.00001));
  
  f4->Draw();
  TLatex latex;
  latex.SetTextSize(0.037);
  latex.SetTextColor(kRed);
  latex.DrawLatex(0.1,-0.002,"#int_{0}^{1} b(y) Acc(y) dy = -0.002633");
  gPad->SetGridx();
  gPad->SetGridy();
}


