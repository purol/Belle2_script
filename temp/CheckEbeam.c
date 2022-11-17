
int CheckEbeam(){

    TH2D* Ebeamstar2D = new TH2D("Ebeamstar", ";E_{beam}^{*} from Mbc;E_{beam}^{*} from deltaE", 100, 10.3, 10.6, 100, 10.3, 10.6);

    double EBstar = -1;
    double PBstar = -1;
    double Mbc = -1;
    double deltaE = -1;
    double beamE = -1;

    TFile* input_file = new TFile("", "read");

    TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
    TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
    TTree* tree_Btag = (TTree*)input_file->Get("Btag");

    tree_Btag->SetBranchAddress("Btag_useCMSFrame_E", &EBstar);
    tree_Btag->SetBranchAddress("Btag_useCMSFrame_p", &PBstar);
    tree_Btag->SetBranchAddress("Btag_Mbc", &Mbc);
    tree_Btag->SetBranchAddress("Btag_deltaE", &deltaE);
    tree_Btag->SetBranchAddress("Upsilon", &beamE);

    for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
        tree_upsilon->GetEntry(j);
        tree_Bsig->GetEntry(j);
        tree_Btag->GetEntry(j);

        double EbeamstarfromMbc = std::sqrt(Mbc * Mbc - PBstar * PBstar);
        double EbeamstarfromdeltaE = EBstar - deltaE;

        Ebeamstar2D->Fill(EbeamstarfromMbc, EbeamstarfromdeltaE);
    }

    input_file->Close();

    TCanvas* c1 = new TCanvas("c1", 500, 500);
    Ebeamstar2D->Draw("COLZ");
    c1->SaveAs("cheack_Ebeam_star.png");

    return 0;
}
