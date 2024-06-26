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

// datatypes for data products I need
// as well as misc other typedef'd types
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "detdataformats/DetID.hpp"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "EventStructs.hh"

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
  std::map<raw::ChannelID_t, raw::RawDigit::ADCvector_t> rawdigit_channels;

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

  // now we need to handle the raw digits :)
  auto rawdigit_handle = e.getValidHandle<std::vector<raw::RawDigit>>(rawdigit_tag);
  std::vector<raw::RawDigit> rd_vec = *rawdigit_handle;
  for(raw::RawDigit rd : rd_vec){
    raw::RawDigit::ADCvector_t adc_vect = rd.ADCs();
    raw::ChannelID_t chan = rd.Channel();
    if(rawdigit_channels.count(chan) != 0){
      std::cout << "Channels are not unique !!" << std::endl;
    }
    rawdigit_channels[chan] = adc_vect;
  }

}

DEFINE_ART_MODULE(duneana::PlotTPsWaveform)
