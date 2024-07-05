////////////////////////////////////////////////////////////////////////
// Class:       ConversionFactorCalc
// Plugin Type: analyzer (Unknown Unknown)
// File:        ConversionFactorCalc_module.cc
//
// Generated at Fri Jul  5 11:03:01 2024 by ddrobner using cetskelgen
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

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "EventStructs.hh"


#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"

namespace duneana {
  class ConversionFactorCalc;
}


class duneana::ConversionFactorCalc : public art::EDAnalyzer {
public:
  explicit ConversionFactorCalc(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ConversionFactorCalc(ConversionFactorCalc const&) = delete;
  ConversionFactorCalc(ConversionFactorCalc&&) = delete;
  ConversionFactorCalc& operator=(ConversionFactorCalc const&) = delete;
  ConversionFactorCalc& operator=(ConversionFactorCalc&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.

};


duneana::ConversionFactorCalc::ConversionFactorCalc(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void duneana::ConversionFactorCalc::analyze(art::Event const& e)
{
  // Implementation of required member function here.
}

DEFINE_ART_MODULE(duneana::ConversionFactorCalc)
