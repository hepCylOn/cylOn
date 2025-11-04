#include <TFile.h>
#include <TH1F.h>

#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>

TEST_CASE("ROOT write tests") {
  TFile f("test.root", "RECREATE");
  TH1F h("h", "Example histogram", 10, 0, 1);
  h.Fill(0.5);
  h.Write();
  REQUIRE(f.IsOpen());
  f.Close();
  REQUIRE(not f.IsOpen());
}

TEST_CASE("ROOT read tests") {
  TFile f("test.root", "READ");
  TH1F* h = (TH1F*)f.Get("h");
  
  REQUIRE(h->GetEntries() == 1);
  f.Close();
  REQUIRE(not f.IsOpen());
}
