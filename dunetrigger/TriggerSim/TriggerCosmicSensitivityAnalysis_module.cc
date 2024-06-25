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
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
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
#include "TTree.h"
#include "TBranch.h"

// I'm just saying screw it and ignoring what art/larsoft does usually here
#include "TFile.h"

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

  // Note to self: using a struct here would have been way easier than a tuple
  typedef std::tuple<unsigned int, uint32_t, uint16_t> tpEvent_t;
  typedef std::tuple<unsigned int, float, float> electronEvent_t;
  typedef unsigned int TDCTime_t;

  typedef struct tpEvent {
    TDCTime_t tdc_time;
    uint32_t adc_int;
    uint16_t adc_peak;
  };

  typedef struct electronEvent {
    TDCTime_t tdc_time;
    float num_electrons;
    float energy;
  };

  // Required functions.
  void analyze(art::Event const &e) override;

private:

  // Declare member data here.
  art::InputTag tp_tag_;
  art::InputTag tpc_tag_;

  TTree* fTGTree;
  TBranch* fPairBranch;

  unsigned int adc_threshold;
  unsigned int electron_pedestal;

  std::vector<dunedaq::trgdataformats::TriggerPrimitive> fTriggerPrimitives;
  std::vector<uint32_t> fADCIntegrals;

  std::vector<sim::SimChannel> fSimChannels;
  std::pair<std::vector<double>, std::vector<double>> paired_data;

  std::map<raw::ChannelID_t, std::vector<electronEvent>> electron_channels;
  std::map<uint32_t, std::vector<tpEvent>> tp_channels;

  // absolute difference for any unsigned datatype
  // necessary to avoid integer underflow
  template <typename T> T u_absdiff(T a, T b){
    return (a > b) ? a - b : b - a;
  }

  // writing some methods we will use in the analysis below
  // TODO? templatize these to use on TPs/other data products?
  std::vector<electronEvent> get_electron_pulse(std::vector<electronEvent>& channel_data, unsigned int start_time, unsigned int max_time){
    std::vector<electronEvent> pulse_evts;
    for(electronEvent evt : channel_data){
      if(u_absdiff<unsigned int>(evt.tdc_time, start_time) <= max_time){
        pulse_evts.push_back(evt);
      }
    }
    return pulse_evts;
  };

  float find_pulse_peak(std::vector<electronEvent>& channel_data, TDCTime_t start_time, TDCTime_t max_time){
    std::vector<electronEvent> e_pulse = get_electron_pulse(channel_data, start_time, max_time);
    // todo figure out if vectors have nice facilities for this
    unsigned int high = 0;
    for (electronEvent i : e_pulse){
      unsigned int val = i.num_electrons;
      if(val > high){
        high = val;
      }
    } 
    return high;
  };

  std::vector<TDCTime_t> pulse_start_times(std::vector<electronEvent>& channel_data, TDCTime_t max_time){
    std::vector<TDCTime_t> start_times = {channel_data.at(0).tdc_time,};
    for(size_t i = 0; i < (channel_data.size()-1); i++){
      if( ( ( u_absdiff<unsigned int>(start_times.back(), channel_data[i].tdc_time) >= max_time ) || u_absdiff<unsigned int>(channel_data[i].tdc_time, channel_data[i+1].tdc_time) >= max_time ) & ( u_absdiff<unsigned int>(start_times.back(), channel_data[i].tdc_time) >= max_time)){
        start_times.push_back(channel_data[i].tdc_time);
      } 
    }
    return start_times;
  }

  uint32_t count_e_over_threshold(std::vector<electronEvent>& channel_data, TDCTime_t max_time, uint32_t threshold) {
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
      tpc_tag_(p.get<art::InputTag>("tpc_info_tag")),
      adc_threshold(p.get<unsigned int>("adc_threshold")),
      electron_pedestal(p.get<unsigned int>("electron_pedestal"))
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

    //  process all of the TPs into a map
  for(dunedaq::trgdataformats::TriggerPrimitive i : fTriggerPrimitives){
    uint64_t adc_time;
    uint32_t adc_int;
    uint16_t adc_peak;
    int32_t channel;

    adc_time = i.time_start/32;
    adc_int = i.adc_integral;
    adc_peak = i.adc_peak;
    channel = i.channel;

    //tpEvent_t this_tp = {static_cast<TDCTime_t>(adc_time), adc_int, adc_peak};
    tpEvent this_tp;
    this_tp.tdc_time = static_cast<TDCTime_t>(adc_time); 
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

  // now let's do the same for the sim channels
  for(sim::SimChannel chan : fSimChannels){
    raw::ChannelID_t chanid = chan.Channel();
    if(electron_channels.count(chanid) == 0){
      std::vector<electronEvent> new_e_channel; 
      electron_channels[chanid] = new_e_channel;
    }
    for(sim::TDCIDE tdcide : chan.TDCIDEMap()){
      unsigned short tdc = tdcide.first;
      for(sim::IDE k : tdcide.second){
        float numElectrons = (k.numElectrons > electron_pedestal) ? k.numElectrons : 0;
        float energy = k.energy;

        electronEvent this_event;
        this_event.tdc_time = tdc;
        this_event.num_electrons = numElectrons;
        this_event.energy = energy;

        electron_channels[chanid].push_back(this_event);
      }
    }
  }

  size_t num_tps = fTriggerPrimitives.size();


  std::vector<double> electron_thresholds;
  for(uint32_t i = 0; i < 40000; i += 500){
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

  std::ostringstream tg_title;
  tg_title << adc_threshold << " ADC;Electron Threshold;Efficiency";
  //TCanvas* c1 = new TCanvas();
  auto tg = new TGraph(electron_thresholds.size(), electron_thresholds.data(), eff_per_t.data());
  tg->SetTitle(tg_title.str().c_str());
  tg->SetLineColor(2);
  tg->SetLineWidth(3);
  tg->SetMarkerSize(1.25);
  tg->SetMarkerStyle(21);
  tg->SetMarkerColor(1);
  tg->Draw("AC");
  //c1->SaveAs("test_tgraph4.png");

  std::ostringstream tg_outfname;
  tg_outfname << adc_threshold << "_sens.root";
  std::ostringstream tg_name;
  tg_name << "TG_" << adc_threshold;

  std::unique_ptr<TFile> o_file(TFile::Open(tg_outfname.str().c_str(), "RECREATE"));
  tg->Write(tg_name.str().c_str());
  //o_file->WriteObject(&tg, tg_name.str().c_str());
  o_file->Close();


  //std::cout << "Efficiency: " << static_cast<double>(num_tps)/static_cast<double>(total_e_over_thresh) << std::endl;

}

DEFINE_ART_MODULE(duneana::TriggerCosmicSensitivityAnalysis)
