#ifndef DUNETRIGGER_EVENTSTRUCTS_hh
#define DUNETRIGGER_EVENTSTRUCTS_hh

#include <cstdint>
#include <vector>
#include "detdataformats/trigger/TriggerPrimitive.hpp"

namespace duneana {

  typedef unsigned int TDCTime_t;
  typedef dunedaq::trgdataformats::TriggerPrimitive TriggerPrimitive_t;
    
  struct tpEvent {
    TDCTime_t tdc_time;
    TDCTime_t peak_time;
    uint32_t adc_int;
    uint32_t adc_peak;
    uint16_t detid;

    tpEvent(TDCTime_t tdc_time, TDCTime_t peak_time, uint32_t adc_int, uint32_t adc_peak):
    tdc_time(tdc_time),
    peak_time(peak_time),
    adc_int(adc_int),
    adc_peak(adc_peak)
    {
      this->detid = 0;
    }
  };

  struct electronEvent {
    TDCTime_t tdc_time;
    float num_electrons;
    float energy;
    int trackID;

    electronEvent(TDCTime_t tdc_time, float num_electrons, float energy):
    tdc_time(tdc_time),
    num_electrons(num_electrons),
    energy(energy)
    {
      this->trackID = 0;
    }
  };

  struct pulse_t
  {
    TDCTime_t tdctime;
    uint32_t charge;
    float energy;
    int trackID;
    uint32_t channel;
    uint32_t peak_charge;
    float peak_energy;
    TDCTime_t peak_time;

    pulse_t(uint32_t channel):
    channel(channel)
    {
      this->tdctime = 0;
      this->charge = 0;
      this->energy = 0.0f;
      this->peak_charge = 0;
      this->peak_time = 0;
      this->peak_energy = 0.0f;
    }

    pulse_t(TDCTime_t tdctime, TDCTime_t peak_time, uint32_t charge, float energy, uint32_t channel, uint32_t peak_charge):
    tdctime(tdctime),
    charge(charge),
    energy(energy),
    channel(channel),
    peak_charge(peak_charge),
    peak_time(peak_time)
    {
      this->trackID = 0; 
    }
  };

}

#endif