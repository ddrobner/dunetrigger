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

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"

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

private:

  art::InputTag electron_tag;

  art::ServiceHandle<geo::Geometry> geom;
  // Declare member data here.

};


duneana::ElectronDumper::ElectronDumper(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  electron_tag(p.get<art::InputTag>("electron_tag"))  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this
  // module.
  consumes<sim::SimChannel>(electron_tag);
}

void duneana::ElectronDumper::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  if(e.event() == 13){
    std::ofstream outfile;
    outfile.open("electrons_mu_serialized.txt");
    auto e_handle = e.getValidHandle<std::vector<sim::SimChannel>>(electron_tag);
    std::vector<sim::SimChannel> els = *e_handle;
    for(sim::SimChannel ec : els){
      raw::ChannelID_t chan = ec.Channel();
      //auto chan_rop = geom->ChannelToROP(chan);
      //auto chan_type = chan_rop.
      geo::SigType_t sigtype = geom->SignalType(geom->ChannelToROP(chan));
      if(sigtype == geo::_plane_sigtype::kCollection){
        for(sim::TDCIDE tdcide : ec.TDCIDEMap()){
          TDCTime_t ts = tdcide.first;
          for(sim::IDE ide : tdcide.second){
            outfile << "(" << chan << "," << ts << "," << ide.numElectrons << ")" << std::endl;
          }
        }
      }
    }
    outfile.close();
  }
}

DEFINE_ART_MODULE(duneana::ElectronDumper)
