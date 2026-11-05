#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>

#include <filesystem>
#include <string>

#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>


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

TEST_CASE("TTree basic write/read roundtrip", "[TTree]") {

  const std::string fileName = "test_ttree.root";

  // Clean up any previous file
  if (std::filesystem::exists(fileName))
    std::filesystem::remove(fileName);

  // --- create and fill the TTree
  {
    TFile outFile(fileName.c_str(), "RECREATE");
    REQUIRE(outFile.IsOpen());

    TTree tree("Validation", "PixelTrack Validation");

    int kind = 0;
    int coord = 0;
    float binCenter = 0.f;
    float value = 0.f;
    float num = 0.f;
    float den = 0.f;

    tree.Branch("kind", &kind, "kind/I");
    tree.Branch("coord", &coord, "coord/I");
    tree.Branch("binCenter", &binCenter, "binCenter/F");
    tree.Branch("value", &value, "value/F");
    tree.Branch("num", &num, "num/F");
    tree.Branch("den", &den, "den/F");

    // Fill 5 entries
    for (int i = 0; i < 5; ++i) {
      kind = i % 3;
      coord = i % 5;
      binCenter = 0.5f * i;
      value = 0.1f * i;
      num = static_cast<float>(i);
      den = static_cast<float>(i + 1);
      tree.Fill();
    }

    REQUIRE(tree.GetEntries() == 5);
    outFile.Write();
  }

  // --- reopen and read back
  {
    TFile inFile(fileName.c_str(), "READ");
    REQUIRE(inFile.IsOpen());

    TTree* tree = static_cast<TTree*>(inFile.Get("Validation"));
    REQUIRE(tree != nullptr);
    REQUIRE(tree->GetEntries() == 5);

    int kind = 0;
    int coord = 0;
    float binCenter = 0.f;
    float value = 0.f;
    float num = 0.f;
    float den = 0.f;

    tree->SetBranchAddress("kind", &kind);
    tree->SetBranchAddress("coord", &coord);
    tree->SetBranchAddress("binCenter", &binCenter);
    tree->SetBranchAddress("value", &value);
    tree->SetBranchAddress("num", &num);
    tree->SetBranchAddress("den", &den);

    // Check a few entries
    tree->GetEntry(0);
    CHECK(kind == 0);
    CHECK(coord == 0);
    CHECK(value == Catch::Approx(0.0f));

    tree->GetEntry(4);
    CHECK(kind == 1);  // since 4 % 3 == 1
    CHECK(coord == 4);
    CHECK(value == Catch::Approx(0.4f));
    CHECK(den == Catch::Approx(5.0f));

    // Compute a simple sum over entries
    double sum = 0.0;
    for (int i = 0; i < tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      sum += value;
    }
    CHECK(sum == Catch::Approx(1.0f));  // 0.0 + 0.1 + 0.2 + 0.3 + 0.4

    inFile.Close();
  }

  // Clean up
  std::filesystem::remove(fileName);
}

