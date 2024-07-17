////////////////////////////////////////////////////////////////////////
// Class:       HistDistHist
// Plugin Type: analyzer (Unknown Unknown)
// File:        HistDistHist_module.cc
//
// Generated at Wed Jul 17 11:11:43 2024 by ddrobner using cetskelgen
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
#include <algorithm>
#include <art_root_io/TFileService.h>

#include "EventStructs.hh"
#include "detdataformats/trigger/TriggerPrimitive.hpp"

#include <larcore/Geometry/Geometry.h>
#include <larcoreobj/SimpleTypesAndConstants/geo_types.h>
#include "canvas/Persistency/Common/FindOneP.h"

#include "TH1F.h"
#include "TCanvas.h"

namespace duneana {
  class HistDistHist;
}


class duneana::HistDistHist : public art::EDAnalyzer {
public:
  explicit HistDistHist(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HistDistHist(HistDistHist const&) = delete;
  HistDistHist(HistDistHist&&) = delete;
  HistDistHist& operator=(HistDistHist const&) = delete;
  HistDistHist& operator=(HistDistHist&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  //art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<geo::Geometry> geom;
  art::InputTag tp_tag;

  std::vector<uint32_t> peaks = {};

};



duneana::HistDistHist::HistDistHist(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  tp_tag(p.get<art::InputTag>("tp_tag"))  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this
  // module.
  consumes<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>(tp_tag);
}

void duneana::HistDistHist::beginJob(){
}

void duneana::HistDistHist::endJob(){
  TCanvas* c1 = new TCanvas();
  c1->SetCanvasSize(2088, 1416);
  //hits = tfs->make<TH1F>();
  //float vmin = *(std::min_element(peaks.begin(), peaks.end()));
  //float vmax = *(std::max_element(peaks.begin(), peaks.end()));
  TH1F* hits = new TH1F("ADC_Dist", "ADC Distribution;ADC;Count", 150, 0, 8000);
  c1->SetLogy(true);

  for(auto h : peaks){
    hits->Fill(static_cast<float>(h));
  }

  hits->SetTitle("Sim Int ADC Distribution;ADC;Count");
  hits->Draw();
  hits->SetLineWidth(4);
  
  c1->SaveAs("hit_distribution_int_sim_50t.png");

}

void duneana::HistDistHist::analyze(art::Event const& e)
{

  auto rawdigit_handle = e.getValidHandle<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>(tp_tag);
  std::vector<dunedaq::trgdataformats::TriggerPrimitive> tps = *rawdigit_handle;

  for(auto tp : tps){
    auto chan = tp.channel; 
    if(geom->SignalType(geom->ChannelToROP(chan)) == geo::SigType_t::kCollection){
      peaks.push_back(tp.adc_integral);
    }
  }
  // Implementation of required member function here.
}

DEFINE_ART_MODULE(duneana::HistDistHist)
