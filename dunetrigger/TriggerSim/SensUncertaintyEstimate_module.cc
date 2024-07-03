////////////////////////////////////////////////////////////////////////
// Class:       SensUncertaintyEstimate
// Plugin Type: analyzer (Unknown Unknown)
// File:        SensUncertaintyEstimate_module.cc
//
// Generated at Thu Jun 27 09:43:08 2024 by ddrobner using cetskelgen
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

#include "EventStructs.hh"

#include "lardataobj/Simulation/SimChannel.h"

namespace duneana {
  class SensUncertaintyEstimate;
};


class duneana::SensUncertaintyEstimate : public art::EDAnalyzer {
public:
  explicit SensUncertaintyEstimate(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SensUncertaintyEstimate(SensUncertaintyEstimate const&) = delete;
  SensUncertaintyEstimate(SensUncertaintyEstimate&&) = delete;
  SensUncertaintyEstimate& operator=(SensUncertaintyEstimate const&) = delete;
  SensUncertaintyEstimate& operator=(SensUncertaintyEstimate&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  art::InputTag tp_tag;
  art::InputTag tpc_tag;


  std::vector<TriggerPrimitive_t> fTPs;
  std::vector<sim::SimChannel> fSimChannels;

  template <typename T> T u_absdiff(T a, T b){
    return (a > b) ? a - b : b - a;
  }

  // method definitions
  std::vector<electronEvent> get_electron_pulse(std::vector<electronEvent>& channel_data, TDCTime_t start_time, TDCTime_t max_time);
  float find_pulse_peak(std::vector<electronEvent>& channel_data, TDCTime_t start_time, TDCTime_t max_time);
  std::vector<TDCTime_t> pulse_start_times(std::vector<electronEvent>& channel_data, TDCTime_t max_time);
  uint32_t count_e_over_threshold(std::vector<electronEvent>& channel_data, TDCTime_t max_time, uint32_t threshold);

  std::map<raw::ChannelID_t, std::vector<electronEvent>> electron_channels;
  std::map<raw::ChannelID_t, std::vector<tpEvent>> tp_channels;
};


duneana::SensUncertaintyEstimate::SensUncertaintyEstimate(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    tp_tag(p.get<art::InputTag>("tp_tag")),
    tpc_tag(p.get<art::InputTag>("tpc_tag"))
    // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this
  // module.
  consumes<std::vector<TriggerPrimitive_t>>(tp_tag);
  consumes<std::vector<sim::SimChannel>>(tpc_tag);
}

std::vector<duneana::electronEvent> duneana::SensUncertaintyEstimate::get_electron_pulse(std::vector<electronEvent>& channel_data, TDCTime_t start_time, TDCTime_t max_time){
  std::vector<electronEvent> pulse_evts;
  for(electronEvent evt : channel_data){
    if(u_absdiff<unsigned int>(evt.tdc_time, start_time) <= max_time){
      pulse_evts.push_back(evt);
    }
  }
  return pulse_evts;
}


float duneana::SensUncertaintyEstimate::find_pulse_peak(std::vector<electronEvent>& channel_data, TDCTime_t start_time, TDCTime_t max_time){
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
}

std::vector<duneana::TDCTime_t> duneana::SensUncertaintyEstimate::pulse_start_times(std::vector<electronEvent>& channel_data, TDCTime_t max_time){
  std::vector<TDCTime_t> start_times = {channel_data.at(0).tdc_time,};
  for(size_t i = 0; i < (channel_data.size()-1); i++){
    if( ( ( u_absdiff<unsigned int>(start_times.back(), channel_data[i].tdc_time) >= max_time ) || u_absdiff<unsigned int>(channel_data[i].tdc_time, channel_data[i+1].tdc_time) >= max_time ) & ( u_absdiff<unsigned int>(start_times.back(), channel_data[i].tdc_time) >= max_time)){
      start_times.push_back(channel_data[i].tdc_time);
    } 
  }
  return start_times;
}

uint32_t duneana::SensUncertaintyEstimate::count_e_over_threshold(std::vector<electronEvent>& channel_data, TDCTime_t max_time, uint32_t threshold){
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

void duneana::SensUncertaintyEstimate::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  // now we have to read the TPs and SimChannels in

  auto fTPHandle = e.getValidHandle<std::vector<TriggerPrimitive_t>>(tp_tag);
  fTPs = *fTPHandle;

  for(TriggerPrimitive_t i : fTPs){
    uint64_t adc_time = i.time_start/32;
    uint32_t adc_int = i.adc_integral;
    uint16_t adc_peak = i.adc_peak;
    int32_t channel = i.channel;

    tpEvent this_tp;
    this_tp.tdc_time = static_cast<TDCTime_t>(adc_time);
    this_tp.adc_int = adc_int;
    this_tp.adc_peak = adc_peak;

    if(tp_channels.count(channel) == 0){
      std::vector<tpEvent> new_ch_vect;
      tp_channels[channel] = new_ch_vect;
    } 

    tp_channels[channel].push_back(this_tp);
  }

  auto fSChandle = e.getValidHandle<std::vector<sim::SimChannel>>(tpc_tag);
  fSimChannels = *fSChandle;

  for(sim::SimChannel chan : fSimChannels){
    raw::ChannelID_t chanid = chan.Channel();

    if(electron_channels.count(chanid) == 0){
      std::vector<electronEvent> new_e_channel;
      electron_channels[chanid] = new_e_channel;
    }

    for(sim::TDCIDE tdcide : chan.TDCIDEMap()){
      unsigned short tdc = tdcide.first;
      for(sim::IDE k : tdcide.second){
        float numElectrons = k.numElectrons;
        float energy = k.energy;

        electronEvent this_event;
        this_event.tdc_time = tdc;
        this_event.num_electrons = numElectrons;
        this_event.energy = energy;

        electron_channels[chanid].push_back(this_event);
      }
    }
  }

  // Ok that's done. Now time for the new stuff
  // the basic idea here is we go
  // over each channel and determine how pulses have multiple TPs
  size_t num_extra_tps = 0;
  for(auto& tp_pair : tp_channels){
    std::vector<tpEvent> tps = tp_pair.second;
    for(size_t i = 1; i<tps.size(); i++){
      if(u_absdiff(tps.at(i).tdc_time, tps.at(i-1).tdc_time) <= 75){
        num_extra_tps++;
      }
    }
  }

  size_t num_tps_total = fTPs.size();
  double extra_tps = static_cast<double>(num_extra_tps)/static_cast<double>(num_tps_total);
  std::cout << "PERCENT EXTRA TPs: " << extra_tps << std::endl;

}

DEFINE_ART_MODULE(duneana::SensUncertaintyEstimate)
