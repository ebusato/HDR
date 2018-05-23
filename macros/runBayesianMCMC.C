///////////////////////////////////////////////////////////
// Usage : 
//  > source setupROOT.sh
//  > source setupTIFOSI.sh
//  > root -l runBayesianMCMC.C
//  or
//  > root -l 'runBayesianMCMC.C(true)'
///////////////////////////////////////////////////////////

void setGraphStyle(TGraph* g, int color, TString yTitle)
{
  g->SetMarkerSize(0.5);
  g->SetMarkerStyle(8);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->GetXaxis()->SetTitle("iteration");
  
  g->GetXaxis()->SetTitleOffset(1.35);
  g->GetXaxis()->SetLabelSize(0.085);
  g->GetYaxis()->SetTitle(yTitle);
  g->GetYaxis()->SetTitleOffset(0.46);
  if(TString(g->GetName()).Contains("Mean")) {
    g->GetYaxis()->SetLabelSize(0.085);
    g->GetYaxis()->SetTitleSize(0.09);
    g->GetXaxis()->SetLabelOffset(2);
    g->GetXaxis()->SetTitleSize(0.);
  }
  else {
    g->GetYaxis()->SetLabelSize(0.075);
    g->GetYaxis()->SetTitleSize(0.08);
    g->GetXaxis()->SetLabelOffset(0.01);
    g->GetXaxis()->SetTitleSize(0.08);
  }
  g->GetYaxis()->SetLabelOffset(0.006);
  g->GetYaxis()->CenterTitle();
  g->GetXaxis()->CenterTitle();
}

void createASCIIFile(const TString filename)
{
  ofstream myfile;
  myfile.open (filename);
  myfile << "+bg Bkg 10 3 \n";
  myfile << "+sig Sig 1 \n";
  myfile << "+data 10 \n";
  myfile.close();
}

TCanvas* makePlots(TString suffix)
{
  TFile* f1 = new TFile(Form("fileMeanStdDevVsIteration1_%s.root",suffix.Data()),"read");
  TFile* f2 = new TFile(Form("fileMeanStdDevVsIteration2_%s.root",suffix.Data()),"read");
  TFile* f3 = new TFile(Form("fileMeanStdDevVsIteration3_%s.root",suffix.Data()),"read");
  TFile* f4 = new TFile(Form("fileMeanStdDevVsIteration4_%s.root",suffix.Data()),"read");
  
  TGraph* gMean1 = (TGraph*) f1->Get("gMean");
  TGraph* gMean2 = (TGraph*) f2->Get("gMean");
  TGraph* gMean3 = (TGraph*) f3->Get("gMean");
  TGraph* gMean4 = (TGraph*) f4->Get("gMean");

  gMean1->GetYaxis()->SetRangeUser(0,38.5);
  setGraphStyle(gMean1,kBlack,"mean");
  setGraphStyle(gMean2,kRed,"mean");
  setGraphStyle(gMean3,kBlue,"mean");
  setGraphStyle(gMean4,kGreen+2,"mean");

  TGraph* gStdDev1 = (TGraph*) f1->Get("gStdDev");
  TGraph* gStdDev2 = (TGraph*) f2->Get("gStdDev");
  TGraph* gStdDev3 = (TGraph*) f3->Get("gStdDev");
  TGraph* gStdDev4 = (TGraph*) f4->Get("gStdDev");

  gStdDev1->GetYaxis()->SetRangeUser(0,13.5);
  setGraphStyle(gStdDev1,kBlack,"standard deviation");
  setGraphStyle(gStdDev2,kRed,"standard deviation");
  setGraphStyle(gStdDev3,kBlue,"standard deviation");
  setGraphStyle(gStdDev4,kGreen+2,"standard deviation");

  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);

  TCanvas* c = new TCanvas(Form("c%s",suffix.Data()),Form("c%s",suffix.Data()),900,500);

  TPad *pPad1=0,*pPad2=0;
  pPad1=new TPad(Form("p1%s",suffix.Data()),"",0.00,0.52,1.00,0.95);
  pPad1->SetBottomMargin(0.035);
  
  pPad2=new TPad(Form("p2%s",suffix.Data()),"",0.00,0.05,1.00,0.52);
  pPad2->SetTopMargin(0);
  pPad2->SetBottomMargin(0.21);

  pPad2->Draw();
  pPad1->Draw();

  pPad2->cd();
  gPad->SetLogx();
  gStdDev1->Draw("apl");
  gStdDev2->Draw("plsame");
  gStdDev3->Draw("plsame");
  gStdDev4->Draw("plsame");

  pPad1->cd();
  gPad->SetLogx();
  gMean1->Draw("apl");
  gMean2->Draw("plsame");
  gMean3->Draw("plsame");
  gMean4->Draw("plsame");

  return c;
}

void runBayesianMCMC(bool recreateMarkovChains=false)
{
  gSystem->Load("Models/Model_C");
  gSystem->Load("BayesianMCMC/BayesianMCMC_C");
  
  if(recreateMarkovChains) {
    model::Model mod(model::Model::normal,4);
    createASCIIFile("fileTemp.dat");
    mod.addChannel("ch1","fileTemp.dat");
    mod.printChannel("ch1",1);
    mod.makeModel();
    
    system("rm -f fileTemp.dat");

    BayesianMCMC bay(mod.getWorkspace());
    bay.setNumIters(1e5);
    //bay.setPriorPdf(BayesianMCMC::Exponential);
    bay.setSeed(1);
    bay.computeLimit(40);
    bay.saveMeanStdDevVsIterationToFile("fileMeanStdDevVsIteration1");
    bay.makePlots(true);
    bay.setSeed(2);
    bay.computeLimit(40);
    bay.saveMeanStdDevVsIterationToFile("fileMeanStdDevVsIteration2");
    bay.makePlots(true);
    bay.setSeed(7);
    bay.computeLimit(40);
    bay.saveMeanStdDevVsIterationToFile("fileMeanStdDevVsIteration3");
    bay.makePlots(true);
    bay.setSeed(4);
    bay.computeLimit(40);
    bay.saveMeanStdDevVsIterationToFile("fileMeanStdDevVsIteration4");
    bay.makePlots(true);
  }

  //TCanvas* cMu = makePlots("mu");
  TCanvas* cBkg = makePlots("BkgYield_ch1");

  //cMu->SaveAs("cMu.png");
  cBkg->SaveAs("cBkg.png");
}
