
Float_t *getRMS_Vpp(Float_t *yVals, Float_t Vpphi);

void readtree(){
  TFile * nfile= new TFile("20171001143118vpol.root");

  Int_t samprate_in_ghz=2;
  Int_t size = 32768;
  //  Int_t band = 100000000
  Float_t band = samprate_in_ghz*100000000;
  Float_t timebase = (Float_t)size/(band*.00001);

  cout<<timebase<<endl;
  
  //TChain * tree= new TChain("tree");
  TTree * t= (TTree*)nfile->Get("tree");
  Float_t CH1[32768];
  Float_t vpp;
  Float_t rms;
  Float_t time_arr[32768];
  // Float_t t;//timestamp
  Float_t dummy=0;
  for (Int_t it=0;it<size;it++){
    dummy=(16.34/(Float_t)32768)*((Float_t)it);
    time_arr[it]=dummy;
    //cout<<dummy<<endl;
  }
  
   t->SetBranchAddress("CH1",CH1);
  // tree->SetBranchAddress("vpp",&vpp);
  // tree->SetBranchAddress("rms",&rms);

  Int_t totent=t->GetEntries();
  Int_t highbin=0;
  Float_t max=0;
  cout<<totent<<endl;
  
  for (Int_t ievt=0;ievt<totent;ievt++){
    t->GetEntry(ievt);
    //  max=TMath::MaxElement(32768,CH1);
    // cout<<max<<endl;
    
    if(ievt<totent){
      Float_t *rms_vpp1=getRMS_Vpp(CH1,0);
      highbin=rms_vpp1[2];
      //cout<<"The time is "<<time_arr[highbin]<<endl;
      //  cout<<time_arr[highbin]<<endl;

      Float_t *rms_vpp2=getRMS_Vpp(CH1,(2./3.0)*(rms_vpp1[0]));
      highbin=rms_vpp2[2];
      //cout<<"The time is "<<time_arr[highbin]<<endl;
      //cout<<time_arr[highbin]<<endl;
      
      Float_t *rms_vpp3=getRMS_Vpp(CH1,(1./3.0)*(rms_vpp2[0]));
      highbin=rms_vpp3[2];
      //cout<<"The time is "<<time_arr[highbin]<<endl;
      cout<<time_arr[highbin]<<endl;
    

    }
    //cout<<rms_vpp[0]<<" "<<rms_vpp[1]<<endl;
   
  }
  
}

//FFTtools::FFTtools() 
Float_t *getRMS_Vpp(Float_t *yVals, Float_t Vpphi)
//Float_t FFTgetWaveformSNR(TGraph *gr,Float_t &peakToPeak,Float_t &rms)
{
  Float_t peakToPeak;
  Float_t rms;

 //Int_t nBins = gr->GetN();
   Int_t nBins = 32768;
   //  Float_t *xVals = gr->GetX();
  //Float_t *yVals = gr->GetY();

   Float_t mean=0.;
   Float_t meanSq=0.;
   Int_t nRMS=32768;
   Float_t avgp2p=0;

   for(Int_t i=0;i<nRMS;i++){
     mean+=yVals[i];
     meanSq+=yVals[i]*yVals[i];
   }
   mean/=static_cast<Float_t>(nRMS);
   meanSq/=static_cast<Float_t>(nRMS);

   Int_t trending=3;
   Float_t p2p=0;
   Int_t firstBin=0;
   Float_t y;
   Float_t dumpk=-1000000;
   Int_t hibin=0;
   
   for(Int_t i=0;i<nBins;i++){
     y=yVals[i];
     if(i>0){
       if(y<yVals[i-1] && trending==0){
         if(TMath::Abs(y-yVals[firstBin]>p2p)){
           p2p=TMath::Abs(y-yVals[firstBin]);
         }
       }
       else if(y<yVals[i-1] && (trending==1 || trending==2)){
         trending=0;
         firstBin=i-1;
         if(TMath::Abs(y-yVals[firstBin]>p2p)){
           p2p=TMath::Abs(y-yVals[firstBin]);
         }
       }
       else if(y>yVals[i-1] && (trending==0 || trending==2)){
         trending=1;
         firstBin=i-1;
         if(TMath::Abs(y-yVals[firstBin]>p2p)){
           p2p=TMath::Abs(y-yVals[firstBin]);
         }
       }
       else if(y>yVals[i-1] && trending==1){
         if(TMath::Abs(y-yVals[firstBin]>p2p)){
           p2p=TMath::Abs(y-yVals[firstBin]);
         }
       }
       else if(y==yVals[i-1]){
         trending=2;
       }
       else if(trending==3){
         if(y<yVals[i-1]){
           trending=0;
           firstBin=0;
         }
         if(y>yVals[i-1]){
           trending=1;
           firstBin=0;
         }
       }
       else{
         std::cout << "trending cock up!" << std::endl;
         std::cout << "y " << y << " yVals[i] " << yVals[i] << " yVals[i-1] " << yVals[i-1] << std::endl;
      
       }
     }
     if(Vpphi>0){
     if(p2p>dumpk && p2p<Vpphi+0.001){
       hibin=i;
       dumpk=p2p;
     }
     }

     if(Vpphi==0){
     if(p2p>dumpk){
       hibin=i;
       dumpk=p2p;
     }
     }
   }

   p2p/=2.;

   rms=sqrt(meanSq-mean*mean);
   peakToPeak = p2p;
   //avgp2p+ = p2p/100;
   // cout<<avgp2p<<endl;
   //cout<<p2p<<endl;
   // cout<<ta<<endl;
   // cout<<rms<<endl;
   // cout<<'a'<<endl;
   //   cout<<"Compare "<<dumpk<<" "<<p2p<<" "<<hibin<<endl;
   Float_t *output=new Float_t[2];
   output[0]=p2p;
   output[1]=rms;
   output[2]=hibin;
   // output[0]=0;
   //output[1]=0;
   return output;
 
 }

//Get timestamp at 1/3 pvoltage//
//for (Int_t ievt=0;ievt<totent;ievt++){
// p2p1/3=
