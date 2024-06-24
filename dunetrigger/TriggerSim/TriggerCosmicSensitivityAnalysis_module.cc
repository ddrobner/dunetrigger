////////////////////////////////////////////////////////////////////////
// Class:       TriggerCosmicSensitivityAnalysis
// Plugin Type: analyzer (Unknown Unknown)
// File:        TriggerCosmicSensitivityAnalysis_module.cc
//
// Generated at Wed Jun 12 13:59:53 2024 by ddrobner using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

// base from cetskelgen
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"


#include "lardataobj/RawData/RawDigit.h"

// dune trigger structs
#include "detdataformats/trigger/TriggerPrimitive.hpp"

#include "lardataobj/Simulation/SimChannel.h"

#include "TGraph.h"
#include "TCanvas.h"

#include <fstream>
#include <stdint.h>
#include <map>
#include <cstdlib>
#include <ranges>

namespace duneana
{
  class TriggerCosmicSensitivityAnalysis;
};

class duneana::TriggerCosmicSensitivityAnalysis : public art::EDAnalyzer
{
public:
  explicit TriggerCosmicSensitivityAnalysis(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TriggerCosmicSensitivityAnalysis(TriggerCosmicSensitivityAnalysis const &) = delete;
  TriggerCosmicSensitivityAnalysis(TriggerCosmicSensitivityAnalysis &&) = delete;
  TriggerCosmicSensitivityAnalysis &operator=(TriggerCosmicSensitivityAnalysis const &) = delete;
  TriggerCosmicSensitivityAnalysis &operator=(TriggerCosmicSensitivityAnalysis &&) = delete;

  typedef std::tuple<unsigned int, uint32_t, uint16_t> tpEvent_t;
  typedef std::tuple<unsigned int, float> electronEvent_t;
  typedef unsigned int TDCTime_t;

  // Required functions.
  void analyze(art::Event const &e) override;

private:

  // Declare member data here.
  art::InputTag tp_tag_;
  art::InputTag tpc_tag_;

  std::vector<dunedaq::trgdataformats::TriggerPrimitive> fTriggerPrimitives;
  std::vector<uint32_t> fADCIntegrals;

  std::vector<sim::SimChannel> fSimChannels;

  std::map<raw::ChannelID_t, std::vector<electronEvent_t>> electron_channels;
  std::map<uint32_t, std::vector<tpEvent_t>> tp_channels;

  // absolute difference for any unsigned datatype
  // necessary to avoid integer underflow
  template <typename T> T u_absdiff(T a, T b){
    return (a > b) ? a - b : b - a;
  }

  // writing some methods we will use in the analysis below
  // TODO? templatize these to use on TPs/other data products?
  std::vector<electronEvent_t> get_electron_pulse(std::vector<electronEvent_t>& channel_data, unsigned int start_time, unsigned int max_time){
    std::vector<electronEvent_t> pulse_evts;
    for(electronEvent_t evt : channel_data){
      if(u_absdiff<unsigned int>(std::get<0>(evt), start_time) <= max_time){
        pulse_evts.push_back(evt);
      }
    }
    return pulse_evts;
  };

  float find_pulse_peak(std::vector<electronEvent_t>& channel_data, TDCTime_t start_time, TDCTime_t max_time){
    std::vector<electronEvent_t> e_pulse = get_electron_pulse(channel_data, start_time, max_time);
    // todo figure out if vectors have nice facilities for this
    unsigned int high = 0;
    for (electronEvent_t i : e_pulse){
      unsigned int val = std::get<1>(i);
      if(val > high){
        high = val;
      }
    } 
    return high;
  };

  std::vector<TDCTime_t> pulse_start_times(std::vector<electronEvent_t>& channel_data, TDCTime_t max_time){
    std::vector<TDCTime_t> start_times = {std::get<TDCTime_t>(channel_data.at(0)),};
    for(size_t i = 0; i < (channel_data.size()-1); i++){
      if( ( ( u_absdiff<unsigned int>(start_times.back(), std::get<0>(channel_data[i])) >= max_time ) || u_absdiff<unsigned int>(std::get<0>(channel_data[i]), std::get<0>(channel_data[i+1])) >= max_time ) & ( u_absdiff<unsigned int>(start_times.back(), std::get<0>(channel_data[i])) >= max_time)){
        start_times.push_back(std::get<0>(channel_data[i]));
      } 
    }
    return start_times;
  }

  uint32_t count_e_over_threshold(std::vector<electronEvent_t>& channel_data, TDCTime_t max_time, uint32_t threshold) {
    uint32_t e_over_threshold = 0U;
    std::vector<TDCTime_t> pulse_starts = pulse_start_times(channel_data, max_time);
    for (TDCTime_t t_start : pulse_starts){
      float e_peak = find_pulse_peak(channel_data, t_start, max_time);
      if (e_peak >= static_cast<float>(threshold)){
        e_over_threshold++;
      }
    }
    return e_over_threshold;
  }

};

duneana::TriggerCosmicSensitivityAnalysis::TriggerCosmicSensitivityAnalysis(fhicl::ParameterSet const &p)
    : EDAnalyzer{p},
      tp_tag_(p.get<art::InputTag>("tp_tag")),
      tpc_tag_(p.get<art::InputTag>("tpc_info_tag"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this
  // module.
  consumes<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>(tp_tag_);
  consumes <std::vector<sim::SimChannel>>(tpc_tag_);
}

void duneana::TriggerCosmicSensitivityAnalysis::analyze(art::Event const &e)
{
  // Implementation of required member function here.

  // get TPs
  auto fTPHandle = e.getValidHandle<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>(tp_tag_);
  fTriggerPrimitives = *fTPHandle;

  auto fChanHandle = e.getValidHandle<std::vector<sim::SimChannel>>(tpc_tag_);
  fSimChannels = *fChanHandle;

  // now let's do the same for the sim channels
  for(sim::SimChannel chan : fSimChannels){
    raw::ChannelID_t chanid = chan.Channel();
    if(electron_channels.count(chanid) == 0){
      std::vector<electronEvent_t> new_e_channel; 
      electron_channels[chanid] = new_e_channel;
    }
    for(sim::TDCIDE tdcide : chan.TDCIDEMap()){
      unsigned short tdc = tdcide.first;
      for(sim::IDE k : tdcide.second){
        float numElectrons = k.numElectrons;
        electron_channels[chanid].push_back({tdc, numElectrons});
      }
    }
  }

  size_t num_tps = fTriggerPrimitives.size();


/*
  uint32_t total_e_over_thresh = 0U;
  for(const auto& sm_pair : electron_channels){
    //auto cur_channel = sm_pair.first;
    auto cur_data = sm_pair.second;
    //std::vector<electronEvent_t> &current_e_channel = sm_pair.second;
    if(cur_data.size() != 0){
      total_e_over_thresh += count_e_over_threshold(cur_data, static_cast<TDCTime_t>(75), 1000U);
    }
  }
*/

  std::vector<double> electron_thresholds;
  for(uint32_t i = 0; i < 20000; i += 500){
    electron_thresholds.push_back(static_cast<double>(i));
  }

  std::vector<double> eff_per_t;
  for (uint32_t thresh : electron_thresholds){
    uint32_t c_el_over_t = 0;
    for(const auto& sm_pair : electron_channels){
      auto cur_data = sm_pair.second;
      if(cur_data.size() != 0) {
        c_el_over_t += count_e_over_threshold(cur_data, static_cast<TDCTime_t>(75), thresh);
      }
    }

    double ratio = static_cast<double>(num_tps)/static_cast<double>(c_el_over_t);
    if (ratio >= 1.0){
      eff_per_t.push_back(1.0);
    } else{
      eff_per_t.push_back(ratio);
    }
  }

  TCanvas* c1 = new TCanvas();
  auto tg = new TGraph(electron_thresholds.size(), electron_thresholds.data(), eff_per_t.data());
  tg->SetTitle("Trigger Sensitivity (200 ADC Threshold);Electron Threshold;Efficiency");
  tg->SetLineColor(2);
  tg->SetLineWidth(3);
  tg->SetMarkerSize(1.25);
  tg->SetMarkerStyle(21);
  tg->SetMarkerColor(1);
  tg->Draw("AC");
  c1->SaveAs("test_tgraph4.png");


  //std::cout << "Efficiency: " << static_cast<double>(num_tps)/static_cast<double>(total_e_over_thresh) << std::endl;

}

DEFINE_ART_MODULE(duneana::TriggerCosmicSensitivityAnalysis)
