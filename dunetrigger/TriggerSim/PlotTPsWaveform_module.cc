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

// art TFileService
#include "art_root_io/TFileService.h"

// datatypes for data products I need
// as well as misc other typedef'd types
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "detdataformats/DetID.hpp"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "EventStructs.hh"

#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"

// TEMPORARY, let's do this and then try and use the TFileService
#include "TFile.h"


namespace duneana {
  class PlotTPsWaveform;
};


class duneana::PlotTPsWaveform : public art::EDAnalyzer {
public:
  explicit PlotTPsWaveform(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PlotTPsWaveform(PlotTPsWaveform const&) = delete;
  PlotTPsWaveform(PlotTPsWaveform&&) = delete;
  PlotTPsWaveform& operator=(PlotTPsWaveform const&) = delete;
  PlotTPsWaveform& operator=(PlotTPsWaveform&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  // ART input tags
  art::InputTag rawdigit_tag;
  art::InputTag tp_tag;

  std::map<uint32_t, std::vector<tpEvent>> tp_channels;
  std::map<uint32_t, raw::RawDigit::ADCvector_t> rawdigit_channels;

  short pedestal = 0;

};


duneana::PlotTPsWaveform::PlotTPsWaveform(fhicl::ParameterSet const& p)
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

void duneana::PlotTPsWaveform::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  auto fTPHandle = e.getValidHandle<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>(tp_tag);
  std::vector<dunedaq::trgdataformats::TriggerPrimitive> fTPs = *fTPHandle;

  // First let's process the TPs into a vector of tpEvent structs
  for(dunedaq::trgdataformats::TriggerPrimitive i : fTPs) {
    // not dividing by 32 here since we are comparing to raw digits
    uint64_t adc_time = i.time_start;
    uint32_t adc_int = i.adc_integral;
    uint16_t adc_peak = i.adc_peak;
    int32_t channel = i.channel;

    // cant do this in one line with POD types apparently?
    tpEvent this_tp;
    this_tp.tdc_time = static_cast<TDCTime_t>(adc_time/32);
    this_tp.adc_int = adc_int;
    this_tp.adc_peak = adc_peak;

    if(tp_channels.count(channel) == 0){
      std::vector<tpEvent> new_ch_vect;
      new_ch_vect.push_back(this_tp);
      tp_channels[channel] = new_ch_vect;
    } else{
      tp_channels[channel].push_back(this_tp);
    }
  }

  // now we need to handle the raw digits :)
  auto rawdigit_handle = e.getValidHandle<std::vector<raw::RawDigit>>(rawdigit_tag);
  std::vector<raw::RawDigit> rd_vec = *rawdigit_handle;
  for(raw::RawDigit rd : rd_vec){
    raw::RawDigit::ADCvector_t adc_vect = rd.ADCs();
    uint32_t chan = static_cast<uint32_t>(rd.Channel());
    if(rawdigit_channels.count(chan) != 0){
      std::cout << "Channels are not unique !!" << std::endl;
    }
    rawdigit_channels[chan] = adc_vect;
  }

  // let's find channels with a lot of TPs now to look for good candidates
  for(auto& tp_pair : tp_channels){
    uint32_t chan = tp_pair.first;
    size_t num_tps = tp_pair.second.size();
    if(num_tps > 5){
      std::cout << "Channel " << chan << ": " << num_tps << " TPs" << std::endl; 
    }
  }

  std::vector<uint32_t> high_tp_channels = {9330, 9406, 2306, 4172, 9527};

/*
  uint32_t t_channel_id = high_tp_channels.at(0);
  raw::RawDigit::ADCvector_t t_chan_v = rawdigit_channels[t_channel_id];
  std::vector<int32_t> t_chan_i(t_chan_v.size());
  for(size_t i = 0; i < t_chan_v.size(); i++){
    t_chan_i[i] = t_chan_v.at(i) - 900;
  }
  std::vector<int32_t> rd_timestamps(t_chan_v.size());
  for(size_t i = 0; i < t_chan_v.size(); i++){
    rd_timestamps[i] = i;
  }

  // make the tgraph, just testing for now
  auto tg_rds = new TGraph(rd_timestamps.size(), rd_timestamps.data(), t_chan_i.data());
  tg_rds->SetLineColor(2);
  tg_rds->SetLineWidth(3);
  //tg_rds->Draw("AL");

  auto tg_tps = new TGraph();
  for(tpEvent tp : tp_channels[t_channel_id]){
    double time = static_cast<double>(tp.tdc_time);
    double peak = static_cast<double>(tp.adc_peak);
    tg_tps->AddPoint(time, peak);
  }
  tg_tps->SetMarkerSize(1.25);
  tg_tps->SetMarkerStyle(21);
  tg_tps->SetMarkerColor(1);

  //TCanvas* c1 = new TCanvas();
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(tg_rds, "l");
  mg->Add(tg_tps, "p");
  mg->Draw("a");

  //tg_rds->Write("Test_TG");
  //mg->Write("Test_MG");
*/

  //std::unique_ptr<TFile> o_file(TFile::Open("test_plot_rd.root", "RECREATE"));
  art::ServiceHandle<art::TFileService> tfs;

  //std::vector<TMultiGraph*> graphs;
  // ok, now time to make multiple TGraphs for each channel
  for(uint32_t chan : high_tp_channels){
    size_t channel_size = rawdigit_channels[chan].size(); 
    std::vector<int32_t> rd_data(channel_size);
    for(size_t i = 0; i < channel_size; i++){
      rd_data[i] = rawdigit_channels[chan].at(i) - pedestal;
    }
    std::vector<int32_t> rd_timestamps(channel_size);
    for(size_t i = 0; i < channel_size; i++){
      rd_timestamps[i] = i;
    }
    TCanvas* c1 = new TCanvas();
    c1->SetCanvasSize(2088,1416);
    auto tg_rds = new TGraph(channel_size, rd_timestamps.data(), rd_data.data());
    tg_rds->SetLineColor(2);
    tg_rds->SetLineWidth(3);
    auto tg_tps = new TGraph();
    for(tpEvent tp : tp_channels[chan]){
      double time = static_cast<double>(tp.tdc_time);
      double peak = static_cast<double>(tp.adc_peak);
      tg_tps->AddPoint(time, peak);
    }
    tg_tps->SetMarkerSize(2.0);
    tg_tps->SetMarkerStyle(21);
    tg_tps->SetMarkerColor(1);

    //auto tg_el = new TGraph();

    TMultiGraph* mg;
    mg = tfs->make<TMultiGraph>();
    mg->Add(tg_rds, "l");
    mg->Add(tg_tps, "p");
    std::ostringstream tgraph_title;
    tgraph_title << "Channel " << chan << " TPs and RDs;Timestamp;ADC Counts";
    mg->SetTitle(tgraph_title.str().c_str());
    mg->Draw("a");
    //graphs.push_back(mg);
    std::ostringstream tgraph_img_name;
    tgraph_img_name << "TG_" << chan <<".png";
    c1->SaveAs(tgraph_img_name.str().c_str());
    std::ostringstream tgraph_name;
    tgraph_name << "TG_" << chan;
    mg->Write(tgraph_name.str().c_str());
  }

  //o_file.Close();

}

DEFINE_ART_MODULE(duneana::PlotTPsWaveform)
