///////////////////////////////////////////////////////////
// Usage : 
//  > root -l load.C runSignificance.C
///////////////////////////////////////////////////////////

#if defined EXECUTABLE || defined __CLING__

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cstdlib>

#include <TStopwatch.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TArrow.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TMath.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TString.h>
#include <TH1.h>

#include "OpTHyLiC.h"

#endif


using namespace std;
using namespace OTH;
using namespace RooFit;

void setPlotStyle(TPad* c, RooPlot* frame, const std::string Xtitle, const std::string Ytitle)
{
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.17);
  
  frame->SetTitle("");
  frame->GetXaxis()->SetRangeUser(-70,50);
  frame->GetXaxis()->SetTitle(Xtitle.c_str());
  frame->GetYaxis()->SetTitle(Ytitle.c_str());
  frame->GetXaxis()->SetTitleSize(0.05);
  frame->GetXaxis()->SetTitleOffset(1.42);
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetXaxis()->SetLabelOffset(0.017);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetTitleOffset(0.74);
  frame->GetYaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTickLength(0);
}

void run(const OTH::SystType systtype, const OTH::StatType stattype, const int randomEngine, const std::string file1, const std::string file2, const std::string file3, const std::string file4, const std::string file5, const std::string outputFileName, const std::string text, const float scale, const float minY, const std::string textoptionnal1="", const std::string textoptionnal2="") {
 
  //OpTHyLiC oth(OTH::SystMclimit,OTH::StatNormal,OTH::STD_mt19937);
  OpTHyLiC oth(systtype,stattype,randomEngine);
  oth.addChannel("ch1",file1);
  if(file2!="")
    oth.addChannel("ch2",file2);
  if(file3!="")
    oth.addChannel("ch3",file3);
  if(file4!="")
    oth.addChannel("ch4",file4);
  if(file5!="")
    oth.addChannel("ch5",file5);
  
  oth.printSamples();
  
  const int nbExp=1e6;
  
  oth.setSigStrength(1);
  oth.generateDistrLLR(nbExp);
  const double pval=oth.pValueData();
  const double zval=sqrt(2.)*TMath::ErfInverse(1-2*pval);
  cout << endl << "Results: " << endl;
  cout << " -> p-value=" << pval << endl;
  cout << " -> significance=" << zval << endl;

  TH1F* hLLRb=(TH1F*) oth.getHisto(OpTHyLiC::hLLRb);
  hLLRb->Rebin(20);
  
  RooRealVar qmu("qmu","qmu",-70,60);
  RooDataHist* dataHistb = new RooDataHist("dataHistb","dataHistb",RooArgList(qmu),hLLRb);
  RooHistPdf* pdfHistb = new RooHistPdf("pdfHistb","pdfHistb",RooArgSet(qmu),*dataHistb);

  double llr=oth.computeLLRdata();

  TCanvas *c1 = new TCanvas("c1", "c1",600,600);
  c1->SetLogy();
  RooPlot* frame = qmu.frame();
  setPlotStyle(c1, frame, "q_{#mu}", "p(q_{ #mu}| #mu' = 0 )");
  pdfHistb->plotOn(frame,LineColor(14),Range(-70,50),Normalization(1,RooAbsReal::Raw));
  pdfHistb->plotOn(frame,LineColor(14),Range(-70,50),Normalization(scale,RooAbsReal::Raw),Range(-70,llr),FillColor(14),FillStyle(3005),DrawOption("F"),VLines());
  frame->GetYaxis()->SetRangeUser(minY,0.09);
  frame->Draw();
  
  TLatex latex0;
  latex0.SetTextSize(0.045);
  latex0.DrawLatexNDC(0.17,0.85,text.c_str());
  latex0.SetTextSize(0.045);
  if(textoptionnal1!="") 
    latex0.DrawLatexNDC(0.17,0.79,textoptionnal1.c_str());
  if(textoptionnal2!="")
  latex0.DrawLatexNDC(0.17,0.73,textoptionnal2.c_str());

  TLatex latex1;
  latex1.SetTextSize(0.05);
  ostringstream oss_pVal;
  oss_pVal << int(pval*1e5)/(float) 1e5;
  latex1.DrawLatexNDC(0.17,0.55,Form("p=%s",oss_pVal.str().c_str()));
  ostringstream oss_Z;
  oss_Z << int(zval*1e2)/(float) 1e2;
  latex1.DrawLatexNDC(0.17,0.48,Form("Z=%s",oss_Z.str().c_str()));

  TArrow* arr = new TArrow(llr,hLLRb->GetMaximum()/6,llr,0,0.02,"|>");
  arr->SetLineWidth(3);
  arr->SetLineColor(kRed);
  arr->SetFillColor(kRed);
  arr->Draw();

  TLatex latex2;
  latex2.SetTextSize(0.045);
  latex2.SetTextColor(kRed);
  latex2.DrawLatex(llr-10,hLLRb->GetMaximum()/5,"q_{ #mu}^{obs}");

  c1->SaveAs(outputFileName.c_str());
 
}

void runSignificance()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  
  run(OTH::SystMclimit,OTH::StatNormal,OTH::STD_mt19937,
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_2/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_3/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_2_METdown_100/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_2_METup_100/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_3_METup_40/CH.All.dat",
      "significanceContactInteractionMcLimitNormal.pdf",
      "OpTHyLiC",
      0.008,
      0.00002,
      "Contact interaction");
  
  run(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_mt19937,
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_2/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_3/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_2_METdown_100/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_2_METup_100/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_3_METup_40/CH.All.dat",
      "significanceContactInteractionPolyExpoGammaHyper.pdf",
      "OpTHyLiC (polyexpo + gamma hyper.)",
      0.0025,
      0.00002,
      "Contact interaction");
  
  run(OTH::SystPolyexpo,OTH::StatGammaUni,OTH::STD_mt19937,
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_2/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_3/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_2_METdown_100/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_2_METup_100/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_3_METup_40/CH.All.dat",
      "significanceContactInteractionPolyExpoGammaUni.pdf",
      "OpTHyLiC (polyexpo + gamma uni.)",
      0.0025,
      0.00002,
      "Contact interaction");
  
  run(OTH::SystPolyexpo,OTH::StatLogN,OTH::STD_mt19937,
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_2/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_3/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_2_METdown_100/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_2_METup_100/CH.All.dat",
      "Files4topContactInt_11dec2014/Ht+Lt_700_NbbJets_3_METup_40/CH.All.dat",
      "significanceContactInteractionPolyExpoLogN.pdf",
      "OpTHyLiC (polyexpo + lognormal.)",
      0.0025,
      0.00002,
      "Contact interaction");
  
  run(OTH::SystMclimit,OTH::StatNormal,OTH::STD_mt19937,
      "Files2UEDRPP_4nov2014/1000GeV/Opti_lowB_lowHt.txt",
      "Files2UEDRPP_4nov2014/1000GeV/Opti_highB_lowHt.txt",
      "Files2UEDRPP_4nov2014/1000GeV/Opti_lowB_highHt_lowMET.txt",
      "Files2UEDRPP_4nov2014/1000GeV/Opti_lowB_highHt_highMET.txt",
      "Files2UEDRPP_4nov2014/1000GeV/Opti_highB_highHt_highMET.txt",
      "significance2UEDRPPMkk1000McLimitNormal.pdf",
      "OpTHyLiC",
      0.006,
      0.00005,
      "2UED/RPP",
      "m_{KK}=1000 GeV");
  
  run(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_mt19937,
      "Files2UEDRPP_4nov2014/1000GeV/Opti_lowB_lowHt.txt",
      "Files2UEDRPP_4nov2014/1000GeV/Opti_highB_lowHt.txt",
      "Files2UEDRPP_4nov2014/1000GeV/Opti_lowB_highHt_lowMET.txt",
      "Files2UEDRPP_4nov2014/1000GeV/Opti_lowB_highHt_highMET.txt",
      "Files2UEDRPP_4nov2014/1000GeV/Opti_highB_highHt_highMET.txt",
      "significance2UEDRPPMkk1000PolyExpoGammaHyper.pdf",
      "OpTHyLiC (polyexpo + gamma hyper.)",
      0.0018,
      0.00001,
      "2UED/RPP",
      "m_{KK}=1000 GeV");
  
  run(OTH::SystPolyexpo,OTH::StatLogN,OTH::STD_mt19937,
      "Files2UEDRPP_4nov2014/1000GeV/Opti_lowB_lowHt.txt",
      "Files2UEDRPP_4nov2014/1000GeV/Opti_highB_lowHt.txt",
      "Files2UEDRPP_4nov2014/1000GeV/Opti_lowB_highHt_lowMET.txt",
      "Files2UEDRPP_4nov2014/1000GeV/Opti_lowB_highHt_highMET.txt",
      "Files2UEDRPP_4nov2014/1000GeV/Opti_highB_highHt_highMET.txt",
      "significance2UEDRPPMkk1000PolyExpoLogN.pdf",
      "OpTHyLiC (polyexpo + lognormal.)",
      0.0018,
      0.00001,
      "2UED/RPP",
      "m_{KK}=1000 GeV");
}
