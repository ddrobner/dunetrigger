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

#include <numeric>
#include <map>

#include "art_root_io/TFileService.h"

#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"

namespace duneana
{
  class ConversionFactorCalc;
}

class duneana::ConversionFactorCalc : public art::EDAnalyzer
{
public:
  explicit ConversionFactorCalc(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ConversionFactorCalc(ConversionFactorCalc const &) = delete;
  ConversionFactorCalc(ConversionFactorCalc &&) = delete;
  ConversionFactorCalc &operator=(ConversionFactorCalc const &) = delete;
  ConversionFactorCalc &operator=(ConversionFactorCalc &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;
  void beginJob() override;
  void endJob() override;

private:
  typedef dunedaq::trgdataformats::TriggerPrimitive TriggerPrimitive_t;

  // input tags for the data products we need
  art::InputTag electron_tag;
  art::InputTag tp_tag;

  // getting a handle on the geometry service
  art::ServiceHandle<geo::Geometry> geom;

  // aliasing the enum for collection wire sigtype
  //geo::SigType_t coll_t = geo::_plane_sigtype::kCollection;


  art::ServiceHandle<art::TFileService> tfs;

  TDCTime_t tp_window;

  TH2F* peak_gr;
  TH2F* int_gr;
  //std::vector<double> w_ars;

  uint32_t electron_threshold;
  // abs diff for unsigned integral types template
  template <typename T>
  T u_absdiff(T a, T b)
  {
    return (a > b) ? a - b : b - a;
  }
};

duneana::ConversionFactorCalc::ConversionFactorCalc(fhicl::ParameterSet const &p)
    : EDAnalyzer{p},
      electron_tag(p.get<art::InputTag>("electron_tag")),
      tp_tag(p.get<art::InputTag>("tp_tag")), // ,
      tp_window(p.get<TDCTime_t>("tp_window")),
      electron_threshold(p.get<uint32_t>("electron_threshold"))
// More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this
  // module.
  consumes<std::vector<sim::SimChannel>>("electron_tag");
  consumes<std::vector<TriggerPrimitive_t>>("tp_tag");
}

void duneana::ConversionFactorCalc::beginJob()
{
  //peak_gr = tfs->make<TGraph>();
  peak_gr = tfs->make<TH2F>("Peak_Plot", "Peak ADC vs Electrons;ADC Count;Electrons Deposited",
   // x settings
   50, 0.0f, 600.0f,
   // y settings
   50, 0.0f, 25000);
  int_gr = tfs->make<TH2F>("Int_Plot", "Integral ADC vs Electrons;ADC Count;Electrons Deposited",
   50, 0.0f, 6000.0f,
   50, 0.0f, 150000.0f);
}

void duneana::ConversionFactorCalc::endJob()
{

  //auto cur_style = new TStyle("Default", "Default Style");
  //cur_style->SetOptFit(0001);
  // handling peak
  std::cout << "Peak Fit" << std::endl;
  TCanvas *c1 = new TCanvas();
  c1->SetGrid();
  c1->SetCanvasSize(4176, 1416);
  c1->Divide(2, 1);

  c1->cd(1);
  peak_gr->SetTitle("Peak ADC vs Electrons;ADC Count;Electrons Deposited");
  peak_gr->Fit("pol1");
  peak_gr->Draw("colz");
  peak_gr->SetStats(false);
  /*
  peak_gr->SetMarkerSize(14);
  peak_gr->SetMarkerColor(4);
  peak_gr->GetXaxis()->SetLimits(0,600);
  peak_gr->GetYaxis()->SetLimits(0,25000);
  peak_gr->SetMinimum(0.);
  peak_gr->SetMaximum(25000.);
  //peak_gr->GetHistogram()->SetStats(1);
  */


  TF1* peak_fit = peak_gr->GetFunction("pol1");
  peak_fit->SetLineWidth(8);
  //peak_fit->SetLineStyle(3);

  /*
  TLegend *peak_leg = new TLegend();
  peak_leg->AddEntry(peak_gr, "Peak Data");
  peak_leg->AddEntry(peak_fit, "Peak Linear Fit");
  peak_leg->Draw();
  */

  //peak_gr->Write("Peak_Conversion");

  std::cout << "Int Fit" << std::endl;

  c1->cd(2);
  int_gr->SetTitle("Integral ADC vs Electrons;ADC Count;Electrons Deposited");
  int_gr->Fit("pol1");
  int_gr->Draw("colz");
  int_gr->SetStats(false);
  /*
  int_gr->SetMarkerSize(8);
  int_gr->SetMarkerColor(4);
  int_gr->SetMinimum(0.);
  int_gr->SetMaximum(150000.);
  int_gr->GetXaxis()->SetLimits(0,6000);
  //int_gr->GetHistogram()->SetStats(1);
  */

  TF1* int_fit = int_gr->GetFunction("pol1");
  int_fit->SetLineWidth(8);
  //int_fit->SetLineStyle(3);

  /*
  TLegend* int_leg = new TLegend();
  int_leg->AddEntry(int_gr, "Integral Data");
  int_leg->AddEntry(int_fit, "Integral Linear Fit");
  */

  //int_gr->Write("Int_Conversion");

  c1->Draw();
  c1->SaveAs("conversion_plots_2.png");

  //std::cout << "============== " << std::reduce(w_ars.begin(), w_ars.end())/(static_cast<double>(w_ars.size())) << " =================" << std::endl;
}

void duneana::ConversionFactorCalc::analyze(art::Event const &e)
{
  std::vector<sim::SimChannel> fElectrons = *(e.getValidHandle<std::vector<sim::SimChannel>>(electron_tag));
  std::vector<TriggerPrimitive_t> fTPs = *(e.getValidHandle<std::vector<TriggerPrimitive_t>>(tp_tag));

  std::vector<sim::SimChannel> collection_channels;
  for (sim::SimChannel c : fElectrons)
  {
    raw::ChannelID_t chanid = c.Channel();
    auto chan_type = geom->SignalType(geom->ChannelToROP(chanid));
    if (chan_type == geo::_plane_sigtype::kCollection)
    {
      collection_channels.push_back(c);
    }
  }
  std::map<raw::ChannelID_t, std::vector<electronEvent>> electron_channels = SensUtilityFunctions::process_es_into_map(collection_channels);

  std::map<raw::ChannelID_t, std::vector<tpEvent>> tp_channels = SensUtilityFunctions::process_tps_into_map(fTPs);

  std::map<raw::ChannelID_t, std::vector<pulse_t>> pulses;
  for (auto &sm_pair : electron_channels)
  {
    raw::ChannelID_t chanid = sm_pair.first;
    std::vector<pulse_t> this_chan_pulses = SensUtilityFunctions::channel_pulse_charges(sm_pair.second, electron_threshold, chanid);
    pulses[chanid] = this_chan_pulses;
  }

 for (auto& pv : pulses)
  {
    raw::ChannelID_t chanid = pv.first;

  for (auto tp : tp_channels[chanid])
    {
      for (pulse_t p : pv.second)
      {
        if (u_absdiff(p.peak_time, tp.peak_time) < tp_window)
        {
          peak_gr->Fill(static_cast<double>(tp.adc_peak), static_cast<double>(p.peak_charge));
          int_gr->Fill(static_cast<double>(tp.adc_int), static_cast<double>(p.charge));
          //w_ars.push_back(static_cast<double>(p.energy)/static_cast<double>(p.charge));
        }
      }
    }
  }
}

DEFINE_ART_MODULE(duneana::ConversionFactorCalc)
