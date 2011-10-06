
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Wed Jul 20 17:49:38 2011 by ROOT version5.27/06b)
//   from TTree Tree/EfficiencyMatrix
//   found on file: EfficiencyMatrix_DYToTauTau_test.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("EfficiencyMatrix_DYToTauTau_test.root");
   if (!f) {
      f = new TFile("EfficiencyMatrix_DYToTauTau_test.root");
      f->cd("EfficiencyMatrix_DYToTauTau_test.root:/effMatrix");
   }
   TTree *Tree = (TTree*)gDirectory->Get("Tree");

//Declaration of leaves types
   Double_t        recoTau_p4_fCoordinates_fX[87];
   Double_t        recoTau_p4_fCoordinates_fY[87];
   Double_t        recoTau_p4_fCoordinates_fZ[87];
   Double_t        recoTau_p4_fCoordinates_fT[87];
   Double_t        genTau_p4Vis_fCoordinates_fX[2];
   Double_t        genTau_p4Vis_fCoordinates_fY[2];
   Double_t        genTau_p4Vis_fCoordinates_fZ[2];
   Double_t        genTau_p4Vis_fCoordinates_fT[2];
   vector<int>     genTau_mode;
   vector<bool>    genTau_goodMother;
   vector<int>     recoTau_mode;
   vector<int>     recoTau_matchedGenMode;
   vector<bool>    recoTau_matchedGenGoodMother;

   // Set branch addresses.
   Tree->SetBranchAddress("recoTau_p4",&recoTau_p4_);
   Tree->SetBranchAddress("recoTau_p4.fCoordinates.fX",recoTau_p4_fCoordinates_fX);
   Tree->SetBranchAddress("recoTau_p4.fCoordinates.fY",recoTau_p4_fCoordinates_fY);
   Tree->SetBranchAddress("recoTau_p4.fCoordinates.fZ",recoTau_p4_fCoordinates_fZ);
   Tree->SetBranchAddress("recoTau_p4.fCoordinates.fT",recoTau_p4_fCoordinates_fT);
   Tree->SetBranchAddress("genTau_p4Vis",&genTau_p4Vis_);
   Tree->SetBranchAddress("genTau_p4Vis.fCoordinates.fX",genTau_p4Vis_fCoordinates_fX);
   Tree->SetBranchAddress("genTau_p4Vis.fCoordinates.fY",genTau_p4Vis_fCoordinates_fY);
   Tree->SetBranchAddress("genTau_p4Vis.fCoordinates.fZ",genTau_p4Vis_fCoordinates_fZ);
   Tree->SetBranchAddress("genTau_p4Vis.fCoordinates.fT",genTau_p4Vis_fCoordinates_fT);
   Tree->SetBranchAddress("genTau_mode",&genTau_mode);
   Tree->SetBranchAddress("genTau_goodMother",&genTau_goodMother);
   Tree->SetBranchAddress("recoTau_mode",&recoTau_mode);
   Tree->SetBranchAddress("recoTau_matchedGenMode",&recoTau_matchedGenMode);
   Tree->SetBranchAddress("recoTau_matchedGenGoodMother",&recoTau_matchedGenGoodMother);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// Tree->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   Long64_t nentries = Tree->GetEntries();

   Long64_t nbytes = 0;
//   for (Long64_t i=0; i<nentries;i++) {
//      nbytes += Tree->GetEntry(i);
//   }

