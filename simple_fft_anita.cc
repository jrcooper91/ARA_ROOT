// #include <useful_things.h>
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

//int simple_fft_anita(Int_t event_no){
void simple_fft_anita(){
  // TString file="/home/uzair/Documents/Utah/ELS/20161212180856.root";
   TString file="/home/kuap2/Documents/physics/anita/root/20171001174651hpol.root";

   Int_t event_no=0;
  Int_t samprate_in_ghz=2;
  TString tit="300 Mhz, -10 dBm, Ice ";

  Int_t size = 32768;
  //  Int_t band = 100000000
  Double_t band = samprate_in_ghz*100000000;
  Double_t bandwidth = band/2;

  //  Int_t event_no = 2; //number of events to average over. 
  TString dirname = "";
  TString infile = dirname+file;
  TCanvas *c1 = new TCanvas(tit, tit, 800, 400);
  TPad *c1_1 = new TPad("c1_1", "c1_1",0.01,0.01,0.49,0.99);
  TPad *c1_2 = new TPad("c1_2", "c1_2",0.51,0.01,0.99,0.99);
  c1_1->Draw();
  c1_2->Draw();
  c1_1->cd();

  
 
  TFile *f = TFile::Open(infile);
  TTree *t = (TTree*)f->Get("tree");
  Int_t points = t->GetBranch("CH1")->GetTotBytes()/t->GetEntries()/4;
  //finds largest power of two smaller than points to speed up fft calculation.
  size = pow(2,floor(log(points)/log(2)));
  Int_t pad = (points-size)/2;
  //cout<<size<<endl;
cout<<"GetEntries "<<t->GetEntries()<<endl;
  //the time series trace
 Double_t timebase = (Double_t)size/(band*.00001);
  //Float_t data[points];
  const Int_t npoints=points;
  const Int_t nsize=size;
  Float_t data[npoints];
  // Double_t indices[nsize];

  t->SetBranchAddress("CH3", data);

  Double_t scale = timebase/(Double_t)size;
  
  TH1D *inhist = new TH1D("voltage", tit+"-time", size, 0, timebase);
  TH1D *inhist_v = new TH1D("v", tit+"-time", size, 0, timebase);
  // TH1D *inhist = new TH1D("voltage", tit+"-time", size, 5.5, 6.5);
  //event_no = t->GetEntries();
  //  t->GetEntry(event_no);
  cout<<"timebase "<<timebase<<endl;
  t->GetEntry(event_no);
  for(Int_t i=0;i<=size;i++){
    inhist->SetBinContent(i, (data[i+pad]));
    // cout<<data[i+pad]<<endl;
    //       inhist->SetBinContent(i, (data[i]+inhist->GetBinContent(i)));
  }
  
  
   inhist->SetXTitle("time (#mus)");
   inhist->SetYTitle("voltage (V)");
   inhist->GetYaxis()->SetTitleOffset(1.5);
  
  inhist_v->Scale(1./size);
  cout<<"time: "<<inhist_v->Integral()<<endl;
  inhist->Draw();
   //fhe fft bit 
   c1_2->cd();
  

   TH1 *out = 0;//dummy th1 for the transform
   out = inhist->FFT(out, "mag");
   Int_t nbins = out->GetXaxis()->GetNbins();

   //output histogram
   TH1D *outhist = new TH1D("freq",tit+"-fft",nbins,0,band/100000);
   for (Int_t i=0;i<=nbins;i++) {
     Double_t y = out->GetBinContent(i);
     y = (10.*log10(pow(y, 2.)/50.))+30.;//normalized v->dbW->dbm
     //y=10.*log10(y)+30;
     outhist->SetBinContent(i, y-(10.*log10(band)));//dmb/hz
     Double_t x = y-(10.*log10(band));
     // cout<<x<<endl;
    //	  outhist->SetBinContent(i, y);
     }

   Double_t rebinval = 1024;
      outhist->Rebin(size/rebinval);//rebin to make it easier to read
      outhist->Scale(rebinval/size);
   outhist->GetXaxis()->SetRangeUser(0,(band/200000));
   outhist->GetYaxis()->SetTitle("dBm/Hz");
   //   outhist->GetYaxis()->SetTitle("normalized amplitude");
   outhist->GetXaxis()->SetTitle("freq (MHz)");
   outhist->GetYaxis()->SetTitleOffset(1.5);
   outhist->Draw("HIST");
   cout<<"freq: "<<outhist->Integral()/pow(size, 2)<<endl;
   //cout<<outhist->Integral()<<endl;

   //   out->Draw();
   //f->Close();
   //   cin>>event_no;
   c1->Update();
 

}

