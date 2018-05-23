#include "atlasstyle-00-03-05/AtlasStyle.C"

void SetStyle()
{
   TStyle* atlasStyle = AtlasStyle();
   gROOT->SetStyle("ATLAS");
   gROOT->ForceStyle();
}

void runExclusionPlotRPPMcLimitVsOTH(const std::string file1, const std::string file2, const std::string outPdfFile, const std::string name1="", const std::string name2="")
{
  SetStyle();

  gSystem->Load("ExclusionPlot_C");

  LimitBrasilPlot limitBrasil1(name1,TheoXsec::RPPFullStat);
  limitBrasil1.readFile(file1);
  limitBrasil1.print();
  limitBrasil1.setXtitle("m_{KK} [TeV]");
  limitBrasil1.setYtitle("#sigma#left(pp#rightarrow A^{(1,1)}A^{(1,1)}#right)#times BR#left(A^{(1,1)}A^{(1,1)}#rightarrow t#bar{t}t#bar{t}#right) [pb]");
  limitBrasil1.setYtitleSize(0.045);
  limitBrasil1.setYtitleOffset(1.7);
  limitBrasil1.setYlabelSize(0.045);
  limitBrasil1.setXtitleSize(0.045);
  limitBrasil1.setXtitleOffset(1.18);
  limitBrasil1.setXlabelSize(0.045);
  limitBrasil1.makePlot(false,2e-3,0.3);
  limitBrasil1.getLegend()->SetY1(0.66);

  LimitBrasilPlot limitBrasil2(name2,TheoXsec::RPPFullStat);
  limitBrasil2.readFile(file2);
  limitBrasil2.makePlot(true);
  limitBrasil2.getLegend()->SetY1(0.66);
  
  gPad->SaveAs(outPdfFile.c_str());
}

void runExclusionPlotRPPOTHVsBayesian(const std::string file1, const std::string file2, const std::string file3)
{
  SetStyle();

  gSystem->Load("ExclusionPlot_C");

  LimitBrasilPlot limitBrasil1("OpTHyLiC",TheoXsec::RPPFullStat);
  limitBrasil1.readFile(file1);
  limitBrasil1.print();
  limitBrasil1.setXtitle("m_{KK} [TeV]");
  limitBrasil1.setYtitle("#sigma#left(pp#rightarrow A^{(1,1)}A^{(1,1)}#right)#times BR#left(A^{(1,1)}A^{(1,1)}#rightarrow t#bar{t}t#bar{t}#right) [pb]");
  limitBrasil1.setYtitleSize(0.045);
  limitBrasil1.setYtitleOffset(1.7);
  limitBrasil1.setYlabelSize(0.045);
  limitBrasil1.setXtitleSize(0.045);
  limitBrasil1.setXtitleOffset(1.18);
  limitBrasil1.setXlabelSize(0.045);
  limitBrasil1.makePlot(false,1e-3,0.7);
  limitBrasil1.plotTheoGraph(true);
  limitBrasil1.getLegend()->SetY1(0.66);
  TLegend* legBr = limitBrasil1.getLegend();
  
  LimitVsMass limitVsMass1("bayesian",Limit::obs,TheoXsec::RPPFullStat);
  limitVsMass1.readFile(file2);
  limitVsMass1.print();
  limitVsMass1.makePlot(false);
  limitVsMass1.plotTheoGraph();

  LimitVsMass limitVsMass2("bayesian",Limit::expMed,TheoXsec::RPPFullStat);
  limitVsMass2.readFile(file3);
  limitVsMass2.print();
  limitVsMass2.makePlot(false);
  limitVsMass2.plotTheoGraph();

  legBr->Draw();
  legBr->SetTextSize(0.025);
  legBr->SetX1(0.5100671);
  legBr->SetY1(0.7045455);
  legBr->SetX2(0.9177852);
  legBr->SetY2(0.9213287);

  TLegend* leg = new TLegend(0.5822148,0.5629371,0.9177852,0.6905594);
  leg->SetHeader("Bayesian (uniform prior)");
  leg->SetTextSize(0.025);
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->AddEntry(limitVsMass1.getGraph(),"95% CI observed limit","LEP");                                                                           
  leg->AddEntry(limitVsMass2.getGraph(),"95% CI expected limit","L");                                                                           
  leg->Draw();
  
  gPad->SaveAs("ExclusionPlot_RPPFullStat_OTHVsBayesian.pdf");
}

void runExclusionPlots()
{
  
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  c1->SetLeftMargin(0.16);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.12);

  runExclusionPlotRPPMcLimitVsOTH("resultsOTHlocal/mclimit_normal.txt","resultsOTHlocal/fromDorianMcLimitResultsForComparisonWithOTH_2nov2014.txt","ExclusionPlot_RPPFullStat_McLimitVsOTH.pdf","OpTHyLiC","McLimit");
  //runExclusionPlotRPPMcLimitVsOTH("resultsOTHlocal/mclimit_normal.txt","../../stat/Limit/ComparisonForSameSignAnalysis/RunAsymptotics/root-files/results_2UEDRPP/results_asymptotic.txt");
  
  
  /*
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  c2->SetLeftMargin(0.16);
  c2->SetRightMargin(0.05);
  c2->SetTopMargin(0.05);
  c2->SetBottomMargin(0.12);
  runExclusionPlotRPPOTHVsBayesian("resultsOTHlocal/mclimit_normal.txt","resultsTIFOSIlocal/bayesian_polyexpo_normal.txt","resultsTIFOSIlocal/bayesian_polyexpo_normal_ASIMOV.txt");
  */
  TCanvas* c3 = new TCanvas("c3","c3",600,600);
  c3->SetLeftMargin(0.16);
  c3->SetRightMargin(0.05);
  c3->SetTopMargin(0.05);
  c3->SetBottomMargin(0.12);
  runExclusionPlotRPPMcLimitVsOTH("resultsOTHlocal/mclimit_normal.txt","resultsAsymptoticlocal/results_asymptotic_1200GeVchangedbyhand.txt","ExclusionPlot_RPPFullStat_McLimitVsAsymptotic.pdf","OTH: mclimit + normal","Frequentist: asymptotic");

  TCanvas* c4 = new TCanvas("c4","c4",600,600);
  c4->SetLeftMargin(0.16);
  c4->SetRightMargin(0.05);
  c4->SetTopMargin(0.05);
  c4->SetBottomMargin(0.12);
  runExclusionPlotRPPMcLimitVsOTH("resultsOTHlocal/fromDorianMcLimitResultsForComparisonWithOTH_2nov2014.txt","resultsOTHlocal/polyexpo_gammahyper.txt","ExclusionPlot_RPPFullStat_McLimitVsOTHpolyexpo_gammahyper.pdf","OTH: mclimit + normal","OTH: polyexpo + gamma hyper");

  TCanvas* c5 = new TCanvas("c5","c5",600,600);
  c5->SetLeftMargin(0.16);
  c5->SetRightMargin(0.05);
  c5->SetTopMargin(0.05);
  c5->SetBottomMargin(0.12);
  runExclusionPlotRPPMcLimitVsOTH("resultsOTHlocal/fromDorianMcLimitResultsForComparisonWithOTH_2nov2014.txt","resultsOTHlocal/polyexpo_logN.txt","ExclusionPlot_RPPFullStat_McLimitVsOTHpolyexpo_logN.pdf","OTH: mclimit + normal","OTH: polyexpo + lognormal");
}
