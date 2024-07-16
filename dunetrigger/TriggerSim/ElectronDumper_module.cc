////////////////////////////////////////////////////////////////////////
// Class:       ElectronDumper
// Plugin Type: analyzer (Unknown Unknown)
// File:        ElectronDumper_module.cc
//
// Generated at Wed Jul  3 14:56:29 2024 by ddrobner using cetskelgen
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
#include "detdataformats/trigger/TriggerPrimitive.hpp"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"

#include "EventStructs.hh"

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "EventStructs.hh"

#include <fstream>
#include <map>

namespace duneana {
  class ElectronDumper;
}


class duneana::ElectronDumper : public art::EDAnalyzer {
public:
  explicit ElectronDumper(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ElectronDumper(ElectronDumper const&) = delete;
  ElectronDumper(ElectronDumper&&) = delete;
  ElectronDumper& operator=(ElectronDumper const&) = delete;
  ElectronDumper& operator=(ElectronDumper&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginJob() override;
  void endJob() override;

private:

  geo::SigType_t coll_t = geo::_plane_sigtype::kCollection;

  art::InputTag electron_tag;
  art::InputTag tp_tag;

  art::ServiceHandle<geo::Geometry> geom;
  // Declare member data here.
  std::ofstream electron_outfile;
  std::ofstream tps_outfile;

};

void duneana::ElectronDumper::beginJob(){
  electron_outfile.open("electrons_mu_serialized_one.txt");
  tps_outfile.open("tps_serialized_one.txt");
}

void duneana::ElectronDumper::endJob(){
  electron_outfile.close();
  tps_outfile.close();
}

duneana::ElectronDumper::ElectronDumper(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  electron_tag(p.get<art::InputTag>("electron_tag")),
  tp_tag(p.get<art::InputTag>("tp_tag"))  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this
  // module.
  consumes<sim::SimChannel>(electron_tag);
  consumes<sim::SimChannel>(tp_tag);
}

void duneana::ElectronDumper::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  std::vector<sim::SimChannel> els = *(e.getValidHandle<std::vector<sim::SimChannel>>(electron_tag));
  std::vector<dunedaq::trgdataformats::TriggerPrimitive> tps = *(e.getValidHandle<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>(tp_tag));

  std::map<raw::ChannelID_t, std::vector<tpEvent>> tp_channels = SensUtilityFunctions::process_tps_into_map(tps);



  for(sim::SimChannel ec : els){
    raw::ChannelID_t chan = ec.Channel();
    //auto chan_rop = geom->ChannelToROP(chan);
    //auto chan_type = chan_rop.
    geo::SigType_t sigtype = geom->SignalType(geom->ChannelToROP(chan));
    if(sigtype == coll_t){
      for(sim::TDCIDE tdcide : ec.TDCIDEMap()){
        TDCTime_t ts = tdcide.first;
        for(sim::IDE ide : tdcide.second){
          electron_outfile << "(" << chan << "," << ts << "," << ide.numElectrons << ")" << std::endl;
        }
        // outputting tp channel, time, peak time, peak adc, int adc 
        for(tpEvent t : tp_channels[chan]){
          tps_outfile << "(" << chan << "," << t.tdc_time << "," << t.peak_time << "," << t.adc_peak << "," << t.adc_int << std::endl; 
        }
      }
    }
  }

}

DEFINE_ART_MODULE(duneana::ElectronDumper)
