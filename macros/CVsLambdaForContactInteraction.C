#include <fstream>
#include <iostream>
#include <string>
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "AtlasStyle.h"
#include "TPad.h"
#include "TGraphPainter.h"

using namespace std;

void GetLatex(bool MC,TLatex* &pLat1)
{

   float xp=0.6;
    pLat1->SetNDC();
    pLat1->SetTextColor(1);
    pLat1->SetTextSize(0.05);
    if (MC) pLat1->SetTextSize(0.08); //0.05
    pLat1->SetTextFont(72);
    float alt=0.56;
   pLat1->DrawLatex(xp,alt,"ATLAS");

      pLat1->SetTextFont(42);
      pLat1->SetTextSize(0.04);
    if (MC) pLat1->SetTextSize(0.07); //  /
       TString label="#font[42]{Internal}";
       float deltaX=0.17;// no cambiar!
       pLat1->DrawLatex(xp+deltaX,alt,label);
       if (!MC){
       pLat1->DrawLatex(xp,alt-0.06, "#sqrt{s} = 8 TeV, 20.3 fb^{-1}");
  //     pLat1->DrawLatex(xp,alt-0.09-0.09, "#sqrt{s} = 8 TeV");
       }
       else {
         pLat1->DrawLatex(xp+deltaX,alt-0.09, "MC Z+jets");

       }
}

void CVsLambdaForContactInteraction(double sfObs, double sfExp, bool overlay=false) 
{  
  double pi=3.14159;
  
  double constSF=4*pi*TMath::Sqrt(sfObs);
  double constExp=4*pi*TMath::Sqrt(sfExp);
  
  TF1 *funcExp = new TF1("f","[0]*x*x",2,5);
  funcExp->SetParameter(0,constExp);
  
  TF1 *funcObs = new TF1("f1","[0]*x*x",0,15);
  funcObs->SetParameter(0,constSF);

 TGraph *pExclusion = new TGraph;

 double x=0;
  int z=0;
  for (int i=1;i<=1500;i++)
    {
      x=i*0.01;
      pExclusion->SetPoint(z,x,constSF*x*x);
      ++z;      
    } 

  pExclusion->SetPoint(z,0,constSF*x*x);
  pExclusion->SetPoint(z+1,0,0);
  
  if(!overlay) {
    funcObs->SetLineWidth(3);
    funcExp->SetLineWidth(2);
    funcExp->SetLineColor(1);    
    funcExp->SetLineStyle(7);   
    pExclusion->SetFillColor(2);
    pExclusion->SetFillStyle(3004);
    funcExp->Draw("");  
  }
  else {
      funcObs->SetLineColor(kBlue);
      funcObs->SetLineWidth(3);
      funcExp->SetLineWidth(2);
      funcExp->SetLineColor(kBlue);    
      funcExp->SetLineStyle(9);
      pExclusion->SetFillColor(kBlue);
      pExclusion->SetFillStyle(3005);
      funcExp->Draw("same");  
  }

  funcObs->Draw("C,same");  
   pExclusion->Draw("F,same");

   TLegend *leg=new TLegend(0.4597315,0.3409091,0.8875839,0.4772727);
   leg->SetHeader("Bayesian (uniform prior)");
   leg->SetTextSize(0.028);
   leg->AddEntry(pExclusion,"95 % CI excluded region","f");
   leg->AddEntry(funcObs,"95 % CI observed limit","l");
   leg->AddEntry(funcExp,"95 % CI expected limit","l");
   leg->SetFillColor(0);
   leg->SetTextFont(42);    
   leg->Draw("same");
}

void CVsLambdaForContactInteraction(double sfObs, double sfExp, double sf1p, double sf2p, double sf1m, double sf2m, bool overlay=false, TString legHeader="OpTHyLiC") 
{  
  double pi=3.14159;
  
  double constSF=4*pi*TMath::Sqrt(sfObs);
  double constExp=4*pi*TMath::Sqrt(sfExp);
  double Exp1p=4*pi*TMath::Sqrt(sf1p);
  double Exp1m=4*pi*TMath::Sqrt(sf1m);
  double Exp2p=4*pi*TMath::Sqrt(sf2p);
  double Exp2m=4*pi*TMath::Sqrt(sf2m);
  
  TF1 *funcExp = new TF1("f","[0]*x*x",2,5);
  funcExp->SetParameter(0,constExp);
  
  TF1 *funcObs = new TF1("f1","[0]*x*x",0,15);
  funcObs->SetParameter(0,constSF);
 
  TGraphAsymmErrors *ph = new TGraphAsymmErrors;
  TGraphAsymmErrors *ph2 = new TGraphAsymmErrors;   

  TGraph *pExclusion = new TGraph;

  double x=0;
  int z=0;
  for (int i=1;i<=1500;i++)
    {
      x=i*0.01;
      ph->SetPoint(i-1,x,constExp*x*x);
      ph->SetPointError(i-1,0,0,x*x*(constExp-Exp1m),x*x*(Exp1p-constExp));
      ph2->SetPoint(i-1,x,constExp*x*x);
      ph2->SetPointError(i-1,0,0,x*x*(constExp-Exp2m),x*x*(Exp2p-constExp));
      pExclusion->SetPoint(z,x,constSF*x*x);
      ++z;      
    } 

  funcObs->SetLineWidth(3);
  funcExp->SetLineWidth(2);

  pExclusion->SetPoint(z,0,constSF*x*x);
  pExclusion->SetPoint(z+1,0,0);
  
  ph2->GetYaxis()->SetTitle("|C|");
  ph2->GetYaxis()->SetTitleSize(0.04);
  ph2->GetYaxis()->SetTitleOffset(1.5);
  ph2->GetXaxis()->SetTitle("#Lambda [TeV]");
  ph2->GetXaxis()->SetTitleOffset(1.2);

  ph2->GetYaxis()->SetRangeUser(0,160);
  ph2->GetXaxis()->SetRangeUser(2,5);

  if(!overlay) {
      pExclusion->SetFillColor(2);
      pExclusion->SetFillStyle(3004);
      ph->SetFillColor(kGreen);
      ph->SetFillStyle(1001);
      ph2->SetFillColor(kYellow);
      ph2->SetFillStyle(1001);  
      funcExp->SetLineColor(1);    
      ph2->Draw("A,E3");
      funcExp->SetLineStyle(7);
  }
  else {
    funcObs->SetLineColor(kBlue);
    pExclusion->SetFillColor(kBlue);
    pExclusion->SetLineStyle(9);
    pExclusion->SetFillStyle(3005);
    ph->SetFillStyle(3004);
    ph2->SetFillStyle(3005);  
    funcExp->SetLineColor(13);    
    ph2->Draw("E3,same");
    funcExp->SetLineStyle(7);
  }
  
  ph->Draw("E3,same");
  funcExp->Draw("same");  
  pExclusion->Draw("F,same");
  funcObs->Draw("C,same");  

  TLegend *leg;
  if(!overlay) 
    leg=new TLegend(0.4597315,0.1206294,0.8875839,0.3321678);
  else
    leg=new TLegend(0.5654362,0.3479021,0.9932886,0.5594406);
  leg->SetHeader(legHeader.Data());
  leg->SetTextSize(0.028);
  leg->AddEntry(pExclusion,"95 % CL excluded region","f");
  leg->AddEntry(funcObs,"95 % CL observed limit","l");
  leg->AddEntry(funcExp,"95 % CL expected limit","l");
  leg->AddEntry(ph,"#pm 1#sigma 95 % CL expected limit","f");
  leg->AddEntry(ph2,"#pm 2#sigma 95 % CL expected limit","f");
  leg->SetFillColor(0);
  //leg->SetFillColor(kWhite);
  leg->SetTextFont(42);    

  
  gPad->RedrawAxis();

  leg->Draw("same");
  //  TLatex *pLat2=new TLatex();
  //bool MC=false;
  //GetLatex(MC,pLat2);
       
  return;  
}

void CVsLambdaForContactInteraction() 
{
  TCanvas* c1=new TCanvas( "c1", "c1", 600, 600);
  c1->cd();
  c1->SetLeftMargin(0.12);
  c1->SetTicks();
  CVsLambdaForContactInteraction(1.45155,0.530867,0.814264,1.25033,0.37339,0.267759);  // Daniela's values (mclimit)
  //  CVsLambdaForContactInteraction(1.44608,0.529042,0.814834,1.24391,0.367955,0.296676,true);  // OTH mclimit+normal
  //CVsLambdaForContactInteraction(1.54063,0.548968,true); // OTH polyexpo+normal
  CVsLambdaForContactInteraction(1.6,0.45,true); // bayesian polyexpo+normal
  c1->SaveAs("CVsLambdaForContactInteractionHybridVsBayesian.pdf");

  TCanvas* c2=new TCanvas( "c2", "c2", 600, 600);
  c2->cd();
  c2->SetLeftMargin(0.12);
  c2->SetTicks();
  CVsLambdaForContactInteraction(1.45155,0.530867,0.814264,1.25033,0.37339,0.267759,0,"McLimit");  // Daniela's values (mclimit)
  CVsLambdaForContactInteraction(1.53374,0.455638,0.737753,1.28945,0.328312,0.244553,true,"Frequentist: asymptotic");  // Asymptotic
  c2->SaveAs("CVsLambdaForContactInteractionHybridVsAsymptotic.pdf");

  TCanvas* c3=new TCanvas( "c3", "c3", 600, 600);
  c3->cd();
  c3->SetLeftMargin(0.12);
  c3->SetTicks();
  CVsLambdaForContactInteraction(1.45155,0.530867,0.814264,1.25033,0.37339,0.267759,0,"McLimit");  // Daniela's values (mclimit)
//   CVsLambdaForContactInteraction(1.44608,0.529042,0.814834,1.24391,0.367955,0.296676,0,"OTH: mclimit+normal"); // OTH mclimit+normal
//   CVsLambdaForContactInteraction(1.65293,0.478975,0.737987,1.17081,0.343857,0.280377,true,"OTH: polyexpo+gammahyper");  // OTH: polyexpo+gammahyper
  CVsLambdaForContactInteraction(1.65293,0.478975,0.737987,1.17081,0.343857,0.26775,true,"OTH: polyexpo+gammahyper");  // change -2sig
  c3->SaveAs("CVsLambdaForContactInteractionMcLimitVsOTHpolyexpogammahyper.pdf");

}
