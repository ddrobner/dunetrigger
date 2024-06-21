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

#include <fstream>
#include <stdint.h>
#include <map>
#include <cstdlib>

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

  const int threshold;

  const double ke_to_adc_factor = 41.649;

  std::vector<dunedaq::trgdataformats::TriggerPrimitive> fTriggerPrimitives;
  std::vector<uint32_t> fADCIntegrals;

  // object with simmed e hits
  size_t fNumTPs;
  size_t fNumTPCEvts;

  std::vector<sim::SimChannel> fSimChannels;

  std::vector<double> fEventADCs;

  std::map<raw::ChannelID_t, std::vector<electronEvent_t>> electron_channels;
  std::map<uint32_t, std::vector<tpEvent_t>> tp_channels;


  // writing some methods we will use in the analysis below
  // TODO? templatize these to use on TPs/other data products?
  std::vector<electronEvent_t> get_electron_pulse(std::vector<electronEvent_t>& channel_data, unsigned int start_time, unsigned int max_time){
    std::vector<electronEvent_t> pulse_evts;
    for(electronEvent_t evt : channel_data){
      if(std::abs(static_cast<int>(std::get<0>(evt) - start_time)) <= max_time){
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
  std::vector<TDCTime_t> pulse_start_times(std::vector<electronEvent_t>& channel_data, TDCTime_t start_time, TDCTime_t max_time){
    std::vector<TDCTime_t> start_times = {std::get<TDCTime_t>(channel_data.at(0)),};
    // do more stuff
    return start_times;
  }

};

duneana::TriggerCosmicSensitivityAnalysis::TriggerCosmicSensitivityAnalysis(fhicl::ParameterSet const &p)
    : EDAnalyzer{p},
      tp_tag_(p.get<art::InputTag>("tp_tag")),
      tpc_tag_(p.get<art::InputTag>("tpc_info_tag")),
      threshold(p.get<int>("threshold"))
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

    tpEvent_t this_tp = {static_cast<TDCTime_t>(adc_time), adc_int, adc_peak};

    if(tp_channels.count(channel) == 0){
      std::vector<tpEvent_t> new_ch_vect;
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

}

DEFINE_ART_MODULE(duneana::TriggerCosmicSensitivityAnalysis)
