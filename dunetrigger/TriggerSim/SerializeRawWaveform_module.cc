////////////////////////////////////////////////////////////////////////
// Class:       SerializeRawWaveform
// Plugin Type: analyzer (Unknown Unknown)
// File:        SerializeRawWaveform_module.cc
//
// Generated at Tue Jun 18 16:35:36 2024 by ddrobner using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "detdataformats/DetID.hpp"
#include "canvas/Persistency/Common/FindOneP.h"


#include <fstream>

class SerializeRawWaveform;


class SerializeRawWaveform : public art::EDAnalyzer {
public:
  explicit SerializeRawWaveform(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SerializeRawWaveform(SerializeRawWaveform const&) = delete;
  SerializeRawWaveform(SerializeRawWaveform&&) = delete;
  SerializeRawWaveform& operator=(SerializeRawWaveform const&) = delete;
  SerializeRawWaveform& operator=(SerializeRawWaveform&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  art::InputTag rawdigit_tag_;
  static const unsigned int ADC_SAMPLING_RATE_IN_DTS = 32; //32 DTS time ticks between adc samples
};


SerializeRawWaveform::SerializeRawWaveform(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    rawdigit_tag_(p.get<art::InputTag>("rawdigit_tag"))
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this
  // module.
  consumes<std::vector<raw::RawDigit>>(rawdigit_tag_);
}

void SerializeRawWaveform::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  auto rawdigit_handle = e.getValidHandle<std::vector<raw::RawDigit>>(rawdigit_tag_);

  std::vector<raw::RawDigit> rawdigit_vec = *rawdigit_handle;
  std::cout << "Found " << rawdigit_vec.size() << " raw::RawDigits" << std::endl;

  uint64_t default_timestamp = 0U;
  uint64_t start_time = default_timestamp;
  std::ofstream rawdigit_out;
  rawdigit_out.open("rawdigits_redux.txt");

  for(size_t i_digit = 0; i_digit < rawdigit_vec.size(); ++i_digit){
    auto const& digit = rawdigit_vec[i_digit];

    unsigned int channel_id = digit.Channel();

    std::vector<short> adcs = digit.ADCs();

    for(size_t i_t = 0; i_t < adcs.size(); ++i_t){
      int16_t sample = adcs[i_t];
      uint16_t ctimestamp = start_time + (i_t)*ADC_SAMPLING_RATE_IN_DTS;
      std::cout << static_cast<uint16_t>(sample) + ctimestamp << std::endl;
      std::cout << channel_id << std::endl;
      //rawdigit_out << "(" << i_digit << "," << adc_digits.at(i) << "," << channel_id << ")" << std::endl; 
    }
  }

  rawdigit_out.close();
}

DEFINE_ART_MODULE(SerializeRawWaveform)
