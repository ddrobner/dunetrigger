#ifndef DUNETRIGGER_EVENTSTRUCTS_hh
#define DUNETRIGGER_EVENTSTRUCTS_hh

#include <cstdint>
#include <map>
#include <vector>

#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "lardataobj/Simulation/SimChannel.h"

namespace duneana {

typedef unsigned int TDCTime_t;
typedef dunedaq::trgdataformats::TriggerPrimitive TriggerPrimitive_t;

struct tpEvent {
  TDCTime_t tdc_time;
  TDCTime_t peak_time;
  uint32_t adc_int;
  uint32_t adc_peak;
  uint16_t detid;

  tpEvent() {
    this->tdc_time = 0;
    this->peak_time = 0;
    this->adc_int = 0;
    this->adc_peak = 0;
    this->detid = 0;
  }

  tpEvent(TDCTime_t tdc_time, TDCTime_t peak_time, uint32_t adc_int,
          uint32_t adc_peak)
      : tdc_time(tdc_time), peak_time(peak_time), adc_int(adc_int),
        adc_peak(adc_peak) {
    this->detid = 0;
  }
};

struct electronEvent {
  TDCTime_t tdc_time;
  float num_electrons;
  float energy;
  int trackID;

  electronEvent() {
    this->tdc_time = 0;
    this->num_electrons = 0.0f;
    this->energy = 0.0f;
    this->trackID = 0.0f;
  }

  electronEvent(TDCTime_t tdc_time, float num_electrons, float energy)
      : tdc_time(tdc_time), num_electrons(num_electrons), energy(energy) {
    this->trackID = 0;
  }
};

struct pulse_t {
  TDCTime_t tdctime;
  uint32_t charge;
  float energy;
  int trackID;
  uint32_t channel;
  uint32_t peak_charge;
  float peak_energy;
  TDCTime_t peak_time;

  pulse_t(uint32_t channel) : channel(channel) {
    this->tdctime = 0;
    this->charge = 0;
    this->energy = 0.0f;
    this->peak_charge = 0;
    this->peak_time = 0;
    this->peak_energy = 0.0f;
    this->trackID = 0;
  }

  pulse_t(TDCTime_t tdctime, TDCTime_t peak_time, uint32_t charge, float energy,
          uint32_t channel, uint32_t peak_charge)
      : tdctime(tdctime), charge(charge), energy(energy), channel(channel),
        peak_charge(peak_charge), peak_time(peak_time) {
    this->trackID = 0;
  }
};

// I probably shouldn't include an impl in the header but I've used
// this enough and don't want to make a whole module for it
class SensUtilityFunctions {
public:
  static std::map<raw::ChannelID_t, std::vector<tpEvent>> process_tps_into_map(
      std::vector<dunedaq::trgdataformats::TriggerPrimitive> &tps) {
    std::map<raw::ChannelID_t, std::vector<tpEvent>> processed_tps;
    for (dunedaq::trgdataformats::TriggerPrimitive i : tps) {
      raw::ChannelID_t channel = i.channel;

      tpEvent this_tp(static_cast<TDCTime_t>(i.time_start / 32),
                      static_cast<TDCTime_t>(i.time_peak / 32),
                      static_cast<uint32_t>(i.adc_integral),
                      static_cast<uint32_t>(i.adc_peak));

      if (processed_tps.count(channel) == 0) {
        std::vector<tpEvent> new_ch_vect;
        processed_tps[channel] = new_ch_vect;
      }

      processed_tps[channel].push_back(this_tp);
    }
    return processed_tps;
  }

  static std::map<raw::ChannelID_t, std::vector<electronEvent>>
  process_es_into_map(std::vector<sim::SimChannel> &electrons) {
    std::map<raw::ChannelID_t, std::vector<electronEvent>> processed_es;
    for (sim::SimChannel chan : electrons) {
      raw::ChannelID_t chanid = chan.Channel();
      if (processed_es.count(chanid == 0)) {
        std::vector<electronEvent> new_e_channel;
        processed_es[chanid] = new_e_channel;
      }
      for (sim::TDCIDE tdcide : chan.TDCIDEMap()) {
        TDCTime_t tdc = static_cast<TDCTime_t>(tdcide.first);
        for (sim::IDE k : tdcide.second) {
          electronEvent this_event(tdc, k.numElectrons, k.energy);
          processed_es[chanid].push_back(this_event);
        }
      }
    }
    return processed_es;
  }

  static std::vector<pulse_t>
  channel_pulse_charges(std::vector<electronEvent> &channel_data,
                        uint32_t threshold, uint32_t channel) {
    std::vector<pulse_t> pulse_info;
    // uint32_t hit_t = 1000;
    uint32_t hit_t = threshold;
    bool last_over = false;
    pulse_t this_pulse(channel);
    for (electronEvent e : channel_data) {
      uint32_t current_charge = static_cast<uint32_t>(e.num_electrons);
      float current_energy = e.energy;
      bool is_over = (current_charge > hit_t);
      if (is_over) {
        this_pulse.charge += current_charge;
        this_pulse.energy += current_energy;
        if (!last_over) {
          this_pulse.tdctime = e.tdc_time;
        }
        if (current_charge > this_pulse.peak_charge) {
          this_pulse.peak_charge = current_charge;
          this_pulse.peak_time = e.tdc_time;
        }
        if (current_energy > this_pulse.peak_energy) {
          this_pulse.peak_energy = current_energy;
        }
        last_over = true;
      } else if (last_over && !is_over) {
        // pulse_t this_pulse(start_time, peak_time, total_charge, total_energy,
        // channel, peak_charge);
        pulse_info.push_back(this_pulse);
        pulse_t this_pulse(channel);
      }
    }
    return pulse_info;
  }
};

} // namespace duneana

#endif