////////////////////////////////////////////////////////////////////////
// Class:       PlotTPsWaveform
// Plugin Type: analyzer (Unknown Unknown)
// File:        PlotTPsWaveform_module.cc
//
// Generated at Wed Jun 26 10:48:52 2024 by ddrobner using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/Simulation/SimChannel.h"

// art TFileService
#include "art_root_io/TFileService.h"

// datatypes for data products I need
// as well as misc other typedef'd types
#include "EventStructs.hh"
#include "detdataformats/DetID.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RawData/RawDigit.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"

// TEMPORARY, let's do this and then try and use the TFileService
#include "TFile.h"
#include <larcore/Geometry/Geometry.h>

namespace duneana {
class PlotTPsWaveform;
};

class duneana::PlotTPsWaveform : public art::EDAnalyzer {
public:
  explicit PlotTPsWaveform(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PlotTPsWaveform(PlotTPsWaveform const &) = delete;
  PlotTPsWaveform(PlotTPsWaveform &&) = delete;
  PlotTPsWaveform &operator=(PlotTPsWaveform const &) = delete;
  PlotTPsWaveform &operator=(PlotTPsWaveform &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

private:
  // Declare member data here.
  // ART input tags
  art::InputTag rawdigit_tag;
  art::InputTag tp_tag;

  art::ServiceHandle<geo::Geometry> geom;
  geo::SigType_t coll_t = geo::_plane_sigtype::kCollection;

  //short pedestal = 8300;
  short pedestal = 900;

  template <typename T> T u_absdiff(T a, T b) {
    return (a > b) ? a - b : b - a;
  }
};

duneana::PlotTPsWaveform::PlotTPsWaveform(fhicl::ParameterSet const &p)
    : EDAnalyzer{p}, 
    rawdigit_tag(p.get<art::InputTag>("rawdigit_tag")),
    tp_tag(p.get<art::InputTag>("tp_tag"))
// ,
// More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this
  // module.
  consumes<std::vector<raw::RawDigit>>(rawdigit_tag);
  consumes<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>(tp_tag);
}

void duneana::PlotTPsWaveform::analyze(art::Event const &e) {
  // Implementation of required member function here.

  std::map<uint32_t, raw::RawDigit::ADCvector_t> rawdigit_channels;
  std::map<raw::ChannelID_t, std::vector<tpEvent>> tp_channels;

  // now we need to handle the raw digits :)
  auto rawdigit_handle =
      e.getValidHandle<std::vector<raw::RawDigit>>(rawdigit_tag);
  std::vector<raw::RawDigit> rd_vec = *rawdigit_handle;
  for (raw::RawDigit rd : rd_vec) {
    //if (geom->SignalType(geom->ChannelToROP(rd.Channel())) == coll_t) {
      raw::RawDigit::ADCvector_t adc_vect = rd.ADCs();
      uint32_t chan = static_cast<uint32_t>(rd.Channel());
      rawdigit_channels[chan] = adc_vect;
    //}
  }

  std::vector<dunedaq::trgdataformats::TriggerPrimitive> tp_handle = *(e.getValidHandle<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>(tp_tag));
  tp_channels = SensUtilityFunctions::process_tps_into_map(tp_handle);

  // let's find channels with a lot of TPs now to look for good candidates
  /*
  for(auto& rd_pair : rawdigit_channels){
    uint32_t chan = rd_pair.first;
    if(tp_channels[chan].size() > 0){
      std::cout << "Channel: " << chan << std::endl;
    }
  }
  */
  /*
  for(auto tp : tp_handle){
    std::cout << tp.channel << std::endl;
  }
  */
  std::vector<uint32_t> high_tp_channels = {2244};


  // std::unique_ptr<TFile> o_file(TFile::Open("test_plot_rd.root",
  // "RECREATE"));

  art::ServiceHandle<art::TFileService> tfs;

  // std::vector<TMultiGraph*> graphs;
  //  ok, now time to make multiple TGraphs for each channel
  for (uint32_t chan : high_tp_channels) {


   // std::vector<pulse_t> pulses = SensUtilityFunctions::channel_pulse_charges(
     //   electron_channels[chan], 2, chan);


    size_t channel_size = rawdigit_channels[chan].size();
    std::cout << channel_size << std::endl;
    std::vector<int32_t> rd_data(channel_size);
    for (size_t i = 0; i < channel_size; i++) {
      rd_data[i] = rawdigit_channels[chan].at(i) - pedestal;
    }
    std::vector<int32_t> rd_timestamps(channel_size);
    for (size_t i = 0; i < channel_size; i++) {
      rd_timestamps[i] = i;
    }
    TCanvas *c1 = new TCanvas();
    c1->SetCanvasSize(2088, 1416);
    auto tg_rds =
        new TGraph(channel_size, rd_timestamps.data(), rd_data.data());
    tg_rds->SetLineColor(2);
    tg_rds->SetLineWidth(3);

    /*
    auto tg_el = new TGraph();
    for (electronEvent e : electron_channels[chan]) {
      double time = static_cast<double>(e.tdc_time);
      double deps = static_cast<double>(e.num_electrons) / (41.742);
      tg_el->AddPoint(time, deps);
    }
    tg_el->SetMarkerSize(2.0);
    tg_el->SetMarkerStyle(20);
    tg_el->SetMarkerColor(7);
    */

    TMultiGraph *mg;
    mg = tfs->make<TMultiGraph>();
    mg->Add(tg_rds, "l");
    //mg->Add(tg_el, "p");

    size_t num_tps = tp_channels[chan].size();
    auto tg_tps = new TGraph();
    if (num_tps != 0) {

      for (tpEvent tp : tp_channels[chan]) {
        double time = static_cast<double>(tp.tdc_time);
        double peak = static_cast<double>(tp.adc_peak);
        tg_tps->AddPoint(time, peak);
        std::cout << "TP Peak: " << tp.adc_peak << std::endl;
      }
      tg_tps->SetMarkerSize(2.0);
      tg_tps->SetMarkerStyle(21);
      tg_tps->SetMarkerColor(1);

      //mg->Add(tg_tps, "p");
    }


/*
    TLegend *lg = new TLegend(0.1, 0.8, 0.2, 0.9);
    lg->AddEntry(tg_rds, "Raw Waveform");
    if (num_tps != 0) {
      lg->AddEntry(tg_tps, "Trigger Primitives");
    }
    lg->AddEntry(tg_el, "Electron Deposits");
  */

    std::ostringstream tgraph_title;
    tgraph_title << "Channel " << chan
                 << " RDs;Timestamp;ADC Counts";
    mg->SetTitle(tgraph_title.str().c_str());
    mg->Draw("a");
    mg->GetXaxis()->SetLimits(4450, 4550);
    //lg->Draw();
    auto elaxis_min = (mg->GetYaxis()->GetXmin()) * 0.041742;
    auto elaxis_max = (mg->GetYaxis()->GetXmax()) * 0.041742;
    TGaxis *elaxis =
        new TGaxis(mg->GetXaxis()->GetXmax(), mg->GetYaxis()->GetXmin(),
                   mg->GetXaxis()->GetXmax(), mg->GetYaxis()->GetXmax(),
                   elaxis_min, elaxis_max, 510, "+L");
    elaxis->SetTitleFont(42);
    elaxis->SetTitle("Num. Electrons (x1000)");
    // elaxis->SetTitleOffset(5);
    elaxis->Draw();
    elaxis->SetLabelFont(42);
    elaxis->SetLabelSize(0.03);
    // graphs.push_back(mg);
    std::ostringstream tgraph_img_name;
    tgraph_img_name << "TG_" << chan << "_sim.png";
    c1->SaveAs(tgraph_img_name.str().c_str());
    std::ostringstream tgraph_name;
    tgraph_name << "TG_" << chan;
    mg->Write(tgraph_name.str().c_str());
  }

  // o_file.Close();
}

DEFINE_ART_MODULE(duneana::PlotTPsWaveform)
