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
    uint32_t adc_int;
    uint16_t adc_peak;
    uint16_t detid;
  };

  struct electronEvent {
    TDCTime_t tdc_time;
    float num_electrons;
    float energy;
    int trackID;
  };

}

#endif