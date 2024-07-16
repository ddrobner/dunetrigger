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
#include "art_root_io/TFileService.h"
// ^^ don't need this for now

// dune trigger structs
#include "detdataformats/trigger/TriggerPrimitive.hpp"

// Larsoft formats
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "EventStructs.hh"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"

// ROOT classes
//#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMultiGraph.h"
#include "TH1F.h"
#include "TGraphErrors.h"

// I'm just saying screw it and ignoring what art/larsoft does usually here
#include "TFile.h"

// C++ stdlib stuff we need
#include <math.h>
#include <stdint.h>
#include <map>
#include <cstdlib>
#include <ranges>
#include <filesystem>

namespace duneana
{
  // typedef std::pair<TDCTime_t, uint32_t> pulse_t;

  class TriggerCosmicSensitivityAnalysis;
};

class duneana::TriggerCosmicSensitivityAnalysis : public art::EDAnalyzer
{
public:
  explicit TriggerCosmicSensitivityAnalysis(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  art::ServiceHandle<art::TFileService> tfs;


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
  int binsize;
  int max_electrons;
  unsigned int tp_window;
  size_t overall_tps = 0;
  size_t overall_pulses = 0;

  int min_hits = std::numeric_limits<int>::max();
  int max_hits = 0;

  std::vector<std::pair<float, bool>> total_eff;
  std::vector<double> conv_factors;

  TGraphErrors *tg;
  TFile *o_file;

  art::ServiceHandle<geo::Geometry> geom;
  geo::SigType_t coll_t = geo::_plane_sigtype::kCollection;

  // float adc_to_mev = 95.2f;
  // float adc_to_mev = 200.0f;
  unsigned short adc_tolerance;

  uint32_t electron_threshold;
  float adc_to_mev;

  // absolute difference for any unsigned datatype
  // necessary to avoid integer underflow
  template <typename T>
  T u_absdiff(T a, T b)
  {
    return (a > b) ? a - b : b - a;
  }

  struct effError_t {
    float electrons;
    float efficiency;
    float xErr;
    float yErr;
    int num_hits;

    effError_t(float electrons, float efficiency, float xErr, float yErr):
      electrons(electrons),
      efficiency(efficiency),
      xErr(xErr),
      yErr(yErr)
      {
        this->num_hits = 0;
      }
    effError_t(float electrons, float efficiency, int num_hits):
    electrons(electrons),
    efficiency(efficiency),
    num_hits(num_hits)
    {
      this->xErr = 0.0f;
      this->yErr = 0.0f;
    }
  };
  float compute_xerr(float binsize);
  float compute_yerr(float num_hits, float efficiency);
  float compute_yerr_oneside(float num_hits, float efficiency);

};

duneana::TriggerCosmicSensitivityAnalysis::TriggerCosmicSensitivityAnalysis(fhicl::ParameterSet const &p)
    : EDAnalyzer{p},
      tp_tag_(p.get<art::InputTag>("tp_tag")),
      tpc_tag_(p.get<art::InputTag>("tpc_info_tag")),
      adc_threshold(p.get<unsigned int>("adc_threshold")),
      electron_pedestal(p.get<unsigned int>("electron_pedestal")),
      binsize(p.get<int>("binsize")),
      max_electrons(p.get<int>("max_electrons")),
      tp_window(p.get<unsigned int>("tp_window")),
      adc_tolerance(p.get<unsigned short>("adc_tolerance")),
      electron_threshold(p.get<uint32_t>("electron_threshold")),
      adc_to_mev(p.get<float>("adc_to_mev"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this
  // module.
  consumes<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>(tp_tag_);
  consumes<std::vector<sim::SimChannel>>(tpc_tag_);
}

void duneana::TriggerCosmicSensitivityAnalysis::beginJob()
{
  std::stringstream out_name;
  out_name << "sens_hists/sens_" << adc_threshold << ".root";
  tg = new TGraphErrors();
  //o_file = TFile::Open("sens_hist.root", "UPDATE");
  o_file = TFile::Open(out_name.str().c_str(), "RECREATE");
  std::ostringstream tg_title;
  tg_title << adc_threshold << " ADC;Num Electrons;Efficiency";
  tg->SetTitle(tg_title.str().c_str());
  tg->SetLineColor(2);
  tg->SetLineWidth(8);
}

void duneana::TriggerCosmicSensitivityAnalysis::endJob()
{

  int round_binsize = binsize;
  int cur_point = 0;

  std::vector<float> e_bin_edges;
  for (int i = 0; i <= max_electrons; i += round_binsize)
  {
    e_bin_edges.push_back(static_cast<float>(i));
  }

  std::vector<effError_t> bin_results;
  float last_eff = 0.0f;
  // now let's compute the ratios for each bin
  for (size_t i = 1; i < e_bin_edges.size(); i++)
  {
    float bin_low = e_bin_edges.at(i - 1);
    float bin_high = e_bin_edges.at(i);

    float hits_w_tp = 0.0f;
    float total_hits = 0.0f;

    for (std::pair<float, bool> p : total_eff)
    {
      float num_es = p.first;
      bool has_tp = p.second;
      if ((num_es > bin_low) && (num_es <= bin_high))
      {
        total_hits += 1.0f;
        if (has_tp)
        {
          hits_w_tp += 1.0f;
        }
      }
    }
    float eff = hits_w_tp / total_hits;
    // check for NaN
    if (eff != eff)
    {
      eff = last_eff;
      std::cout << "eff is NaN" << std::endl;
    }
    else
    {
      eff = hits_w_tp / total_hits;
    }
    last_eff = eff;

    effError_t this_result(bin_low, eff, total_hits);
    bin_results.push_back(this_result);
  }
  //std::sort(bin_results.begin(), bin_results.end());

  for (effError_t r : bin_results)
  {

    if(r.num_hits < min_hits && r.num_hits > 0){
	min_hits = r.num_hits;	
    }
    if(r.num_hits > max_hits){
	max_hits = r.num_hits;
    }
    float xErr = 0.0f;
    float yErr = 0.0f;

    xErr = compute_xerr(r.electrons);
    if(r.efficiency <= 0.999f){
      // we can use the efficiency as the parameter p for the binomial, since
      // both the maximum likelihood estimator and method of moments give that
      // p_hat = x/n where x is the number of successes (ie. TPs) 
      yErr = compute_yerr(r.num_hits, r.efficiency);
      // ^^ this gives the uncertainty in the number of hits with TPs, so to get
      // an efficiency we need to normalize by the number of hits
    } else{
			// here we "take the limit" of a binomial approaching q=1 (or p=0) giving
			// us a poisson
      // eff of 1.0 is a degenerate case
			yErr = compute_yerr_oneside(r.num_hits, (r.efficiency == 1.0f) ? 0.99999f : r.efficiency);
    }
    tg->AddPoint(r.electrons, r.efficiency);
    tg->SetPointError(cur_point, xErr, yErr);
    cur_point += 1;
  }

  tg->SetFillColorAlpha(1, 0.15f);
  tg->SetFillStyle(3001);
  std::ostringstream tgname;
  tgname << "TG_" << adc_threshold;
  o_file->cd();
  tg->Write(tgname.str().c_str());
  o_file->Close();

  std::cout << "Min Hits: " << min_hits << " Max Hits: " << max_hits << std::endl;
}

void duneana::TriggerCosmicSensitivityAnalysis::analyze(art::Event const &e)
{

  // Implementation of required member function here.
  std::vector<dunedaq::trgdataformats::TriggerPrimitive> fTriggerPrimitives;
  std::vector<sim::SimChannel> fSimChannels;
  // get TPs
  auto fTPHandle = e.getValidHandle<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>(tp_tag_);
  fTriggerPrimitives = *fTPHandle;

  auto fChanHandle = e.getValidHandle<std::vector<sim::SimChannel>>(tpc_tag_);
  fSimChannels = *fChanHandle;

  overall_tps += fTriggerPrimitives.size();

  // process the TPs into a map
  std::map<raw::ChannelID_t, std::vector<tpEvent>> tp_channels = SensUtilityFunctions::process_tps_into_map(fTriggerPrimitives);

  // find which channels are collection channels
  // to avoid needing the geometry service in the utility classs
  std::vector<sim::SimChannel> collection_channels;
  for (sim::SimChannel chan : fSimChannels)
  {
    raw::ChannelID_t chanid = chan.Channel();
    auto chan_type = geom->SignalType(geom->ChannelToROP(chanid));
    if (chan_type == coll_t)
    {
      collection_channels.push_back(chan);
    }
  }

  // and now process the electron channels as we were before
  std::map<raw::ChannelID_t, std::vector<electronEvent>> electron_channels = SensUtilityFunctions::process_es_into_map(collection_channels);

  std::vector<pulse_t> pulses;

  for (auto &sm_pair : electron_channels)
  {
    uint32_t cur_chan = sm_pair.first;
    std::vector<pulse_t> cur_p = SensUtilityFunctions::channel_pulse_charges(sm_pair.second, electron_threshold, cur_chan);
    pulses.reserve(pulses.size() + std::distance(cur_p.begin(), cur_p.end()));
    pulses.insert(pulses.end(), cur_p.begin(), cur_p.end());
  }


  std::vector<std::pair<float, bool>> deposits_w_tps;

  for (pulse_t p : pulses)
  {
    uint32_t p_charge = p.charge;
    uint32_t p_chan = p.channel;

    bool has_tp = false;

    for (tpEvent tp : tp_channels[p_chan])
    {
      if (u_absdiff(p.peak_time, tp.peak_time) < tp_window)
      {
        has_tp = true;
        overall_pulses += 1;
        break;
      }
    }
    std::pair<float, bool> cur_res(static_cast<float>(p_charge), has_tp);
    deposits_w_tps.push_back(cur_res);
  }


  total_eff.reserve(total_eff.size() + std::distance(deposits_w_tps.begin(), deposits_w_tps.end()));
  total_eff.insert(total_eff.end(), deposits_w_tps.begin(), deposits_w_tps.end());
}

float duneana::TriggerCosmicSensitivityAnalysis::compute_xerr(float binsize){
  return sqrtf(binsize);
}

float duneana::TriggerCosmicSensitivityAnalysis::compute_yerr(float num_hits, float efficiency){
  //return sqrtf(num_hits*efficiency*(1.0f - efficiency));
	return sqrtf(efficiency*(1.0f-efficiency)/num_hits);
}

float duneana::TriggerCosmicSensitivityAnalysis::compute_yerr_oneside(float num_hits, float efficiency){
  //return sqrtf(x_point);
  // this is really (1 - efficiency)*num_hits/num_hits
  return sqrtf(((1.0f - efficiency)));

}

DEFINE_ART_MODULE(duneana::TriggerCosmicSensitivityAnalysis)
