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

#include <map>

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

  typedef dunedaq::trgdataformats::TriggerPrimitive TriggerPrimitive_t;

  // input tags for the data products we need
  art::InputTag electron_tag;
  art::InputTag tp_tag;

  // getting a handle on the geometry service
  art::ServiceHandle<geo::Geometry> geom;

  // aliasing the enum for collection wire sigtype
  geo::SigType_t coll_t = geo::_plane_sigtype::kCollection;

  TDCTime_t tp_window;


  // abs diff for unsigned integral types template
  template <typename T>
  T u_absdiff(T a, T b)
  {
    return (a > b) ? a - b : b - a;
  }

  // storing the output of all of the conversion factors here
  std::vector<double> peak_conversion_factors;
  std::vector<double> integral_conversion_factors;

};


duneana::ConversionFactorCalc::ConversionFactorCalc(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  electron_tag(p.get<art::InputTag>("electron_tag")),
  tp_tag(p.get<art::InputTag>("tp_tag")),  // ,
  tp_window(p.get<TDCTime_t>("tp_window"))
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this
  // module.
  consumes<std::vector<sim::SimChannel>>("electron_tag");
  consumes<std::vector<TriggerPrimitive_t>>("tp_tag");
}

void duneana::ConversionFactorCalc::analyze(art::Event const& e)
{
  std::vector<sim::SimChannel> fElectrons = *(e.getValidHandle<std::vector<sim::SimChannel>>(electron_tag));
  std::vector<TriggerPrimitive_t> fTPs = *(e.getValidHandle<std::vector<TriggerPrimitive_t>>(tp_tag));

  std::map<raw::ChannelID_t, std::vector<electronEvent>> electron_channels;
  std::map<raw::ChannelID_t, std::vector<tpEvent>> tp_channels;
}

DEFINE_ART_MODULE(duneana::ConversionFactorCalc)
