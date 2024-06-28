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


// ART root IO
//#include "art_root_io/TFileService.h"
// ^^ don't need this for now

// dune trigger structs
#include "detdataformats/trigger/TriggerPrimitive.hpp"

// Larsoft formats
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "EventStructs.hh"

// ROOT classes
#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMultiGraph.h"
#include "TH1F.h"

// I'm just saying screw it and ignoring what art/larsoft does usually here
#include "TFile.h"

// C++ stdlib stuff we need
#include <stdint.h>
#include <map>
#include <cstdlib>
#include <ranges>
#include <filesystem>

namespace duneana
{
  typedef std::pair<TDCTime_t, uint32_t> pulse_t;
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

  // Required functions.
  void analyze(art::Event const &e) override;

  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  art::InputTag tp_tag_;
  art::InputTag tpc_tag_;

  unsigned int adc_threshold;
  unsigned int electron_pedestal;

  std::vector<dunedaq::trgdataformats::TriggerPrimitive> fTriggerPrimitives;

  std::vector<sim::SimChannel> fSimChannels;

  std::map<raw::ChannelID_t, std::vector<electronEvent>> electron_channels;
  std::map<uint32_t, std::vector<tpEvent>> tp_channels;

  std::vector<std::pair<float, bool>> total_eff;

  TGraph* tg;
  TFile* o_file;


  // absolute difference for any unsigned datatype
  // necessary to avoid integer underflow
  template <typename T> T u_absdiff(T a, T b){
    return (a > b) ? a - b : b - a;
  }


  std::vector<pulse_t> channel_pulse_charges(std::vector<electronEvent>& channel_data, uint32_t threshold);
  

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

void duneana::TriggerCosmicSensitivityAnalysis::beginJob(){
  //TGraph* tg = new TGraph();
  tg = new TGraph();
  //std::unique_ptr<TFile> o_file(TFile::Open("sens_hist.root", "UPDATE"));
  //TFile* o_file(TFile::Open("sens_hit.root", "UPDATE"));
  o_file = TFile::Open("sens_hist.root", "UPDATE");
  std::ostringstream tg_title;
  tg_title << adc_threshold << " ADC;Num Electrons;Efficiency";
  tg->SetTitle(tg_title.str().c_str());
  tg->SetLineColor(2);
  tg->SetLineWidth(8);
}

void duneana::TriggerCosmicSensitivityAnalysis::endJob(){
  
  std::cout << "Running endJob()" << std::endl;

  std::vector<float> e_bin_edges;
  int max_electrons = 100000;
  int bin_size = 1250;
  for(int i = 0; i < max_electrons; i+= bin_size){
    e_bin_edges.push_back(static_cast<float>(i));
  }

  std::vector<std::pair<float, float>> bin_results;

  float last_eff = 0.0f;
  // now let's compute the ratios for each bin
  for(size_t i = 1; i<e_bin_edges.size(); i++){
    float bin_low = e_bin_edges.at(i-1);
    float bin_high = e_bin_edges.at(i);

    float hits_w_tp = 0.0f;
    float total_hits = 0.0f;

    for(std::pair<float, bool> p: total_eff){
      float num_es = p.first;
      bool has_tp = p.second;
      if( (num_es >= bin_low) && (num_es < bin_high) ){
        total_hits += 1.0f;
        if(has_tp){
          hits_w_tp += 1.0f;
        }
      }
    }
    float eff = hits_w_tp/total_hits;
    // check for NaN
    if(eff != eff){
      eff = last_eff;
    } else {
      eff = hits_w_tp/total_hits;
    }
   last_eff = eff;
    float bin_mid = static_cast<float>(e_bin_edges.at(i-1)) + static_cast<float>(bin_size)/2;
    std::pair<float, float> this_result(bin_mid, eff);
    //std::cout << "MID: " << bin_mid << " EFF: " << eff << std::endl;
    bin_results.push_back(this_result);
  }

  for(std::pair<float, float> r : bin_results){
    tg->AddPoint(r.first, r.second);
  }

  tg->Draw("AL");
  std::ostringstream tgname;
  tgname << "TG_" << adc_threshold;
  o_file->cd();
  tg->Write(tgname.str().c_str());
  o_file->Close();
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

  // ok - so we compute total pulse charge based on that method
  uint32_t electron_threshold = 1;
  std::vector<std::pair<float, bool>> hit_res;
  TDCTime_t tp_window = 100;
  for(auto& sm_pair : electron_channels){
    raw::ChannelID_t chan = sm_pair.first;
    std::vector<electronEvent> chan_data = sm_pair.second;
    //std::cout << "===================== CHANNEL " << chan << " ========================" << std::endl;
    std::vector<pulse_t> pulses = channel_pulse_charges(chan_data, electron_threshold);
    for(auto& p : pulses){
      TDCTime_t p_time = p.first;
      uint32_t num_electrons = p.second;
      bool has_tp = false;
      for(tpEvent t : tp_channels[chan]){
        if(u_absdiff(p_time, t.tdc_time) < tp_window){
          has_tp = true;
          break;
        } 
      }
      std::pair<float,bool> this_res(num_electrons, has_tp);
      hit_res.push_back(this_res);
      }
    }

    total_eff.reserve(total_eff.size() + std::distance(hit_res.begin(), hit_res.end()));
    total_eff.insert(total_eff.end(), hit_res.begin(), hit_res.end());


}

std::vector<duneana::pulse_t> duneana::TriggerCosmicSensitivityAnalysis::channel_pulse_charges(std::vector<electronEvent>& channel_data, uint32_t threshold){
  std::vector<pulse_t> pulse_info;
  uint32_t hit_t = threshold;
  bool last_over = false;
  uint32_t total_charge = 0;
  TDCTime_t start_time;
  for(electronEvent e : channel_data){
    uint32_t current_charge = e.num_electrons;
    //hit_t = static_cast<uint32_t>(current_charge/2);
    bool is_over = (current_charge > hit_t);
    if(is_over){
      if(last_over){
        total_charge += current_charge;
      }else if(!last_over){
        total_charge += current_charge;
        start_time = e.tdc_time;
      }
      last_over = true;
    } else if(last_over && !is_over){
      pulse_t this_pulse(start_time, total_charge);
      pulse_info.push_back(this_pulse);
      //std::cout << "START: " << start_time << " CHARGE: " << total_charge << std::endl;
      total_charge = 0;
      last_over = false;
      //hit_t = threshold;
    } else{
      last_over = false;
      //hit_t = threshold;
    }
  }
  return pulse_info;
}

DEFINE_ART_MODULE(duneana::TriggerCosmicSensitivityAnalysis)
