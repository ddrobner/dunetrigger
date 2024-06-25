#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"

#include <vector>

void mgraph_m(){
    // lol using vector here to be able to iterate nicely
    std::vector<int> thresholds = {100, 150, 200, 250, 300};

    std::vector<TFile> file_refs;
    std::vector<TGraph*> graph_refs;

    gStyle->SetPalette(55);

    for(int i = 0; i < static_cast<int>(thresholds.size()); i++){
        TGraph* cur_gr = nullptr;
        graph_refs.push_back(cur_gr);

        std::ostringstream cur_fname;
        cur_fname << thresholds.at(i) << "_sens.root";
        TFile cur_fi(cur_fname.str().c_str());
        
        std::ostringstream gr_name;
        gr_name << "TG_" << thresholds.at(i);
        cur_fi.GetObject(gr_name.str().c_str(), graph_refs.at(i));
        cur_fi.Close();
    }

    TCanvas* c1 = new TCanvas();
    c1->SetCanvasSize(2088,1416);
    TMultiGraph *mg = new TMultiGraph();
    TLegend *leg = new TLegend(0.85, 0.3, 0.7, 0.15);
    leg->SetFillColor(0);
    leg->SetHeader("Trigger Thresholds");

    // set the colors for each graph individually
    for(size_t i = 0; i < graph_refs.size(); i++){
        std::ostringstream legend_text;
        legend_text << thresholds.at(i) << " ADC";
        //graph_refs.at(i)->SetLineColor(colour);
        graph_refs.at(i)->SetLineWidth(8);
        mg->Add(graph_refs.at(i), "l");
        leg->AddEntry(graph_refs.at(i), legend_text.str().c_str(), "l");
    }

    mg->SetTitle("Trigger Efficiency;Electron Threshold; Efficiency");
    mg->Draw("a pmc plc");
    leg->Draw();
    mg->GetXaxis()->SetLimits(0,25000);
    mg->GetYaxis()->SetLimits(0,1);

    c1->SaveAs("test_mg.png");
}