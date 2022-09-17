#include <iostream>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <vector>
#include <sstream>

TGraphErrors* ObservedGraph(const char* type, std::vector<double> xArray, std::vector<double> yArray, std::vector<double> yErrArray)
{
    int typeID = -1;
    // 0: CLS
    // 1: CLB
    // 2: CLS+B

    if (std::string("CLS") == std::string(type)) typeID = 0;
    else if (std::string("CLB") == std::string(type)) typeID = 1;
    else if ((std::string("CLS+B") == std::string(type)) || (std::string("CLSPLUSB") == std::string(type))) typeID = 2;
    else {
        printf("invalid type!\n");
        exit(1);
    }

    TGraphErrors* graph = new TGraphErrors(static_cast<Int_t>(xArray.size()), &xArray.front(), &yArray.front(), 0, &yErrArray.front());
    TString pValueName = "CLs";
    if (typeID == 1) pValueName = "CLb";
    if (typeID == 2) pValueName = "CLs+b";
    TString name = pValueName + TString("_observed");
    TString title = TString("Observed ") + pValueName;
    graph->SetName(name);
    graph->SetTitle(title);
    graph->SetMarkerStyle(20);
    graph->SetLineWidth(2);

    if (typeID == 0) {
        graph->SetMarkerColor(kRed);
    }
    if (typeID == 1) {
        graph->SetMarkerColor(kBlue + 4);
        // graph->Draw("PL");
    }
    else if (typeID == 2) {
        graph->SetMarkerColor(kBlue);
        graph->SetLineStyle(3);
        // graph->Draw("PL");
    }

    return graph;
}

TMultiGraph* ExpectedGraph(std::vector<double> mu_values, std::vector<double> ExpCLsMedian, std::vector<double> ExpCLs1sigplus, std::vector<double> ExpCLs1sigminus, std::vector<double> ExpCLs2sigplus, std::vector<double> ExpCLs2sigminus)
{
    // create the graphs
    TGraph* g0 = new TGraph;
    TString pValueName = "CLs";
    g0->SetTitle(TString::Format("Expected %s - Median", pValueName.Data()));

    TGraphAsymmErrors* g1 = new TGraphAsymmErrors;
    g1->SetTitle(TString::Format("Expected %s #pm 1 #sigma", pValueName.Data()));

    TGraphAsymmErrors* g2 = new TGraphAsymmErrors;
    g2->SetTitle(TString::Format("Expected %s #pm 2 #sigma", pValueName.Data()));

    for (unsigned int i = 0; i < mu_values.size(); i++) {
        g0->SetPoint(i, mu_values.at(i), ExpCLsMedian.at(i));

        g1->SetPoint(i, mu_values.at(i), ExpCLsMedian.at(i));
        g1->SetPointEYlow(i, ExpCLs1sigminus.at(i)); // -1 sigma errorr
        g1->SetPointEYhigh(i, ExpCLs1sigplus.at(i)); // +1 sigma error

        g2->SetPoint(i, mu_values.at(i), ExpCLsMedian.at(i));
        g2->SetPointEYlow(i, ExpCLs2sigminus.at(i)); // -2 sigma errorr
        g2->SetPointEYhigh(i, ExpCLs2sigplus.at(i)); // +2 sigma error
    }

    TString name = "CLs" + TString("_expected");
    TString title = TString("CLs scan");
    TMultiGraph* graph = new TMultiGraph(name, title);

    // set the graphics options and add in multi graph
    // orderof adding is drawing order
    if (g2) {
        g2->SetFillColor(kYellow);
        graph->Add(g2, "3");
    }
    if (g1) {
        g1->SetFillColor(kGreen);
        graph->Add(g1, "3");
    }
    g0->SetLineStyle(2);
    g0->SetLineWidth(2);
    graph->Add(g0, "L");

    return graph;
}

void DrawPlots(TGraphErrors* GraphObservedCLs, TGraphErrors* GraphObservedCLb, TGraphErrors* GraphObservedCLsplusb, TMultiGraph* GraphExpectedCLs, double conf = 0.9) {
    TCanvas* c = new TCanvas("c", "c", 696, 472);

    // draw observed CLs first
    GraphObservedCLs->Draw("APL");
    TGraphErrors* gplot = GraphObservedCLs;
    gplot->GetHistogram()->SetTitle("CLs scan");

    // draw Expected
    GraphExpectedCLs->Draw();

    // draw line for conf level
  //  GraphExpectedCLs->Draw();

    // draw line for conf level
    if (gplot) {
        double alpha = 1.0 - conf;
        double x1 = gplot->GetXaxis()->GetXmin();
        double x2 = gplot->GetXaxis()->GetXmax();
        TLine* line = new TLine(x1, alpha, x2, alpha);
        line->SetLineColor(kRed);
        line->Draw();
        // put axis labels
        gplot->GetXaxis()->SetTitle("#mu");
        gplot->GetYaxis()->SetTitle("p value");
    }

    // draw observed CLb
    GraphObservedCLb->Draw("PL");

    // draw CLs+b
    GraphObservedCLsplusb->Draw("PL");

    // draw observed CLs again
    GraphObservedCLs->Draw("PL");

    // draw legend
    double y0 = 0.6;
    double verticalSize = 0.3;
    double y1 = y0 + verticalSize;
    TLegend* l = new TLegend(0.6, y0, 0.9, y1);
    l->AddEntry(GraphObservedCLs, "", "PEL");
    l->AddEntry(GraphObservedCLsplusb, "", "PEL");
    //l->AddEntry(GraphObservedCLs, "", "PEL");
    l->AddEntry(GraphObservedCLb, "", "PEL");
    // loop in reverse order (opposite to drawing one)
    int ngraphs = GraphExpectedCLs->GetListOfGraphs()->GetSize();
    for (int i = ngraphs - 1; i >= 0; --i) {
        TObject* obj = GraphExpectedCLs->GetListOfGraphs()->At(i);
        TString lopt = "F";
        if (i == ngraphs - 1) lopt = "L";
        if (obj)  l->AddEntry(obj, "", lopt);
    }
    l->Draw();

    // redraw the axis
    if (gPad) gPad->RedrawAxis();

    c->SaveAs("CLs.png");
}

void ReadTxt(bool IsItFreq, std::vector<double>* mu_values, std::vector<double>* ObsCLss,
    std::vector<double>* ObsCLsErrors, std::vector<double>* ObsCLbs, std::vector<double>* ObsCLbErrors,
    std::vector<double>* ObsCLsPlusbs, std::vector<double>* ObsCLsPlusbErrors, std::vector<double>* ExpCLsMedian,
    std::vector<double>* ExpCLs1sigplus, std::vector<double>* ExpCLs1sigminus, std::vector<double>* ExpCLs2sigplus,
    std::vector<double>* ExpCLs2sigminus) {

    for (int i = 0; i < 101; i++) {
        std::ostringstream streamObj;
        streamObj << std::fixed;
        streamObj << std::setprecision(1);

        std::string fname;
        double mu = i * 0.2;
        streamObj << mu;
        if (IsItFreq) fname = "CLs_freq_" + streamObj.str() + ".txt";
        else fname = "CLs_hyb_" + streamObj.str() + ".txt";

        FILE* fp;
        fp = fopen(fname.c_str(), "r");

        double temp_ObsCLs;
        double temp_ObsCLsError;
        double temp_ObsCLb;
        double temp_ObsCLbError;
        double temp_ObsCLsPlusb;
        double temp_ObsCLsPlusbError;

        double temp_ExpCLsMedian;
        double temp_ExpCLs1sigplus;
        double temp_ExpCLs1sigminus;
        double temp_ExpCLs2sigplus;
        double temp_ExpCLs2sigminus;

        fscanf(fp, "expected CLs median: %lf\n", &temp_ExpCLsMedian);
        fscanf(fp, "expected CLs +1sigma: %lf\n", &temp_ExpCLs1sigplus);
        fscanf(fp, "expected CLs -1sigma: %lf\n", &temp_ExpCLs1sigminus);
        fscanf(fp, "expected CLs +2sigma: %lf\n", &temp_ExpCLs2sigplus);
        fscanf(fp, "expected CLs -2sigma: %lf\n", &temp_ExpCLs2sigminus);
        fscanf(fp, "observed CLs central value: %lf\n", &temp_ObsCLs);
        fscanf(fp, "observed CLs error: %lf\n", &temp_ObsCLsError);
        fscanf(fp, "observed CLb central value: %lf\n", &temp_ObsCLb);
        fscanf(fp, "observed CLb error: %lf\n", &temp_ObsCLbError);
        fscanf(fp, "observed CLs+b central value: %lf\n", &temp_ObsCLsPlusb);
        fscanf(fp, "observed CLs+b error: %lf\n", &temp_ObsCLsPlusbError);

        fclose(fp);

        mu_values->push_back(mu);

        ObsCLss->push_back(temp_ObsCLs);
        ObsCLsErrors->push_back(temp_ObsCLsError);
        ObsCLbs->push_back(temp_ObsCLb);
        ObsCLbErrors->push_back(temp_ObsCLbError);
        ObsCLsPlusbs->push_back(temp_ObsCLsPlusb);
        ObsCLsPlusbErrors->push_back(temp_ObsCLsPlusbError);

        ExpCLsMedian->push_back(temp_ExpCLsMedian);
        ExpCLs1sigplus->push_back(temp_ExpCLs1sigplus);
        ExpCLs1sigminus->push_back(temp_ExpCLs1sigminus);
        ExpCLs2sigplus->push_back(temp_ExpCLs2sigplus);
        ExpCLs2sigminus->push_back(temp_ExpCLs2sigminus);
    }
}

void 

void ReadCLsTextFile(){

	std::vector<double> mu_values;
	std::vector<double> ObsCLss;
	std::vector<double> ObsCLsErrors;
	std::vector<double> ObsCLbs;
	std::vector<double> ObsCLbErrors;
	std::vector<double> ObsCLsPlusbs;
	std::vector<double> ObsCLsPlusbErrors;

    std::vector<double> ExpCLsMedian;
    std::vector<double> ExpCLs1sigplus;
    std::vector<double> ExpCLs1sigminus;
    std::vector<double> ExpCLs2sigplus;
    std::vector<double> ExpCLs2sigminus;

    // read txt files
    ReadTxt(true, &mu_values, &ObsCLss,
        &ObsCLsErrors, &ObsCLbs, &ObsCLbErrors,
        &ObsCLsPlusbs, &ObsCLsPlusbErrors, &ExpCLsMedian,
        &ExpCLs1sigplus, &ExpCLs1sigminus, &ExpCLs2sigplus,
        &ExpCLs2sigminus);

    // get graphs
    TGraphErrors* GraphObservedCLs = ObservedGraph("CLS", mu_values, ObsCLss, ObsCLsErrors);
    TGraphErrors* GraphObservedCLb = ObservedGraph("CLB", mu_values, ObsCLbs, ObsCLbErrors);
    TGraphErrors* GraphObservedCLsplusb = ObservedGraph("CLS+B", mu_values, ObsCLsPlusbs, ObsCLsPlusbErrors);
    TMultiGraph* GraphExpectedCLs = ExpectedGraph(mu_values, ExpCLsMedian, ExpCLs1sigplus, ExpCLs1sigminus, ExpCLs2sigplus, ExpCLs2sigminus);

    // draw
    DrawPlots(GraphObservedCLs, GraphObservedCLb, GraphObservedCLsplusb, GraphExpectedCLs);
}
