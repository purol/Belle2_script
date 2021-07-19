void load_files(const char *dirname, std::vector<string>* names){
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (!file->IsDirectory() && fname.EndsWith(".root")) {
            names->push_back(fname.Data());
         }
      }
   }
}

typedef struct data{
    int __experiment__;
    int __run__;
    int __event__;
    int __candidate__;
    int __ncandidates__;

    int event_info[15];
    // 0: upsilon_experiment, 1: upsilon_run, 2: upsilon_event, 3: upsilon_candidate, 4: upsilon_ncandidates
    // 5: Bsig_experiment, 6: Bsig_run, 7: Bsig_event, 8: Bsig_candidate, 9: Bsig_ncandidates
    // 10: Btag_experiment, 11: Btag_run, 12: Btag_event, 13: Btag_candidate, 14: Btag_ncandidates

    double Upsilon_info[13];
    // 0: Upsilon_isSignal; 1: number of ECL clusters in ROE(cleanMask), 2: number of KLM clusters in ROE
    // 3: energy in ROE(cleanMask), 4: number of tracks in ROE(cleanMask), 5: roeEextra(cleanMask)
    // 6: nROE_NeutralECLClusters(cleanMask), 7: roeNeextra(cleanMask), 8: energy in ROE(cleanMask) at CMS
    // 9: roeExtra(cleanMask) at CMS, 10: roeNeextra(cleanMask) at CMS, 11: nROE_K_S0
    // 12: nROE_pi0, 

    double Bsig_info[14];
    // 0: Bsig_isSignal, 1: Bsig_E, 2: Bsig_E_CMS, 3: Bsig_E_Recoil, 4: Bsig_dmID
    // 5: Bsig_first_daughter's_actPID(3,2) 6: Bsig_first_daughter's_mcPDG
    // 7: Bsig_p, 8: Bsig_p_CMS, 9: Bsig_p_Recoil
    // 10: Kaon dr, 11: Kaon dz, 12: Bsig_first_daughter's eIDBelle
    // 13: Bsig_first_daughter's_muIDBelle

    double Btag_info[7];
    // 0: Btag_isSignal, 1:Btag_dmID, 2: Btag_Mbc, 3: Btag_deltaE
    // 4: Btag_E, 5: Btag_E_CMS, 6: Btag_signalprobability

} Data; 

bool event_info_is_valid(std::queue<Data> TotalData_){
    while(!TotalData_.empty()){
        Data temp_data = TotalData_.front();
        TotalData_.pop();
        if(temp_data.event_info[0] != temp_data.event_info[5] || temp_data.event_info[0] != temp_data.event_info[10] || temp_data.event_info[5] != temp_data.event_info[10]) return false;
        if(temp_data.event_info[1] != temp_data.event_info[6] || temp_data.event_info[1] != temp_data.event_info[11] || temp_data.event_info[6] != temp_data.event_info[11]) return false;
        if(temp_data.event_info[2] != temp_data.event_info[7] || temp_data.event_info[2] != temp_data.event_info[12] || temp_data.event_info[7] != temp_data.event_info[12]) return false;
        if(temp_data.event_info[3] != temp_data.event_info[8] || temp_data.event_info[3] != temp_data.event_info[13] || temp_data.event_info[8] != temp_data.event_info[13]) return false;
        if(temp_data.event_info[4] != temp_data.event_info[9] || temp_data.event_info[4] != temp_data.event_info[14] || temp_data.event_info[9] != temp_data.event_info[14]) return false;
    }
    return true;
}

void Draw_Upsilon_info(std::queue<Data> TotalData_, int index){
    double max = -std::numeric_limits<double>::max();
    double min = std::numeric_limits<double>::max();
    double interval;
    std::queue<double> temp;

    while(!TotalData_.empty()){
        Data temp_data = TotalData_.front();
        TotalData_.pop();
        if(temp_data.Upsilon_info[index] > max) max = temp_data.Upsilon_info[index];
        if(temp_data.Upsilon_info[index] < min) min = temp_data.Upsilon_info[index];
        temp.push(temp_data.Upsilon_info[index]);
    }

    interval = max - min;

    TCanvas* temp_c = new TCanvas("c", "", 1500, 1200);
    TH1F* temp_hist = new TH1F("temp",";;evt",100,min - interval*0.1 ,max + interval*0.1 );
    while(!temp.empty()){
        double temp_data = temp.front();
        temp.pop();
        temp_hist->Fill(temp_data);
    }
    temp_hist->Draw("Hist");
    temp_c->SaveAs( ("Upsilon_" + std::to_string(index) + ".png").c_str() );
    delete temp_hist;
    delete temp_c;

}

void Draw_Bsig_info(std::queue<Data> TotalData_, int index){
    double max = -std::numeric_limits<double>::max();
    double min = std::numeric_limits<double>::max();
    double interval;
    std::queue<double> temp;

    while(!TotalData_.empty()){
        Data temp_data = TotalData_.front();
        TotalData_.pop();
        if(temp_data.Bsig_info[index] > max) max = temp_data.Bsig_info[index];
        if(temp_data.Bsig_info[index] < min) min = temp_data.Bsig_info[index];
        temp.push(temp_data.Bsig_info[index]);
    }

    interval = max - min;

    TCanvas* temp_c = new TCanvas("c", "", 1500, 1200);
    TH1F* temp_hist = new TH1F("temp",";;evt",100,min - interval*0.1,max + interval*0.1);
    while(!temp.empty()){
        double temp_data = temp.front();
        temp.pop();
        temp_hist->Fill(temp_data);
    }
    temp_hist->Draw("Hist");
    temp_c->SaveAs( ("Bsig_" + std::to_string(index) + ".png").c_str() );
    delete temp_hist;
    delete temp_c;

}

void Draw_Bsig_info(std::queue<Data> TotalData_, int index, const char* arg1, const char* arg2){
    double max = -std::numeric_limits<double>::max();
    double min = std::numeric_limits<double>::max();
    double interval;
    std::queue<double> temp;

    while(!TotalData_.empty()){
        Data temp_data = TotalData_.front();
        TotalData_.pop();
        if(temp_data.Bsig_info[index] > max) max = temp_data.Bsig_info[index];
        if(temp_data.Bsig_info[index] < min) min = temp_data.Bsig_info[index];
        temp.push(temp_data.Bsig_info[index]);
    }

    interval = max - min;

    TCanvas* temp_c = new TCanvas("c", "", 1500, 1200);
    TH1F* temp_hist = new TH1F(arg1,arg2,100,min - interval*0.1,max + interval*0.1);
    while(!temp.empty()){
        double temp_data = temp.front();
        temp.pop();
        temp_hist->Fill(temp_data);
    }
    temp_hist->Draw("Hist");
    temp_c->SaveAs( ("Bsig_"+ std::string(arg1) + "_"  + std::to_string(index) + ".png").c_str() );
    delete temp_hist;
    delete temp_c;

}

void Draw_Btag_info(std::queue<Data> TotalData_, int index){
    double max = -std::numeric_limits<double>::max();
    double min = std::numeric_limits<double>::max();
    double interval;
    std::queue<double> temp;

    while(!TotalData_.empty()){
        Data temp_data = TotalData_.front();
        TotalData_.pop();
        if(temp_data.Btag_info[index] > max) max = temp_data.Btag_info[index];
        if(temp_data.Btag_info[index] < min) min = temp_data.Btag_info[index];
        temp.push(temp_data.Btag_info[index]);
    }

    interval = max - min;

    TCanvas* temp_c = new TCanvas("c", "", 1500, 1200);
    TH1F* temp_hist = new TH1F("temp",";;evt",100,min - interval*0.1 ,max + interval*0.1 );
    while(!temp.empty()){
        double temp_data = temp.front();
        temp.pop();
        temp_hist->Fill(temp_data);
    }
    temp_hist->Draw("Hist");
    temp_c->SaveAs( ("Btag_" + std::to_string(index) + ".png").c_str() );
    delete temp_hist;
    delete temp_c;

}

void Draw_Btag_info(std::queue<Data> TotalData_, int index, const char* arg1, const char* arg2){
    double max = -std::numeric_limits<double>::max();
    double min = std::numeric_limits<double>::max();
    double interval;
    std::queue<double> temp;

    while(!TotalData_.empty()){
        Data temp_data = TotalData_.front();
        TotalData_.pop();
        if(temp_data.Btag_info[index] > max) max = temp_data.Btag_info[index];
        if(temp_data.Btag_info[index] < min) min = temp_data.Btag_info[index];
        temp.push(temp_data.Btag_info[index]);
    }

    interval = max - min;

    TCanvas* temp_c = new TCanvas("c", "", 1500, 1200);
    TH1F* temp_hist = new TH1F(arg1,arg2,100,min - interval*0.1 ,max + interval*0.1 );
    while(!temp.empty()){
        double temp_data = temp.front();
        temp.pop();
        temp_hist->Fill(temp_data);
    }
    temp_hist->Draw("Hist");
    temp_c->SaveAs( ("Btag_" + std::string(arg1) + "_" + std::to_string(index) + ".png").c_str() );
    delete temp_hist;
    delete temp_c;

}

void Draw_Mbc_deltaE_Btag(std::queue<Data> TotalData_, int index, bool SetLogy = false){
    TH1F* temp_hist1 = new TH1F("deltaE_total","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-0.5,0.5);
    TH1F* temp_hist2 = new TH1F("Mbc_total","Mbc of B_{tag};M_{bc} [GeV];Num of candidate",100,5.2,5.3);
    TGraph* temp_MbcVSdeltaE_Btag = new TGraph();

    while(!TotalData_.empty()){
        Data temp_data = TotalData_.front();
        TotalData_.pop();
        temp_hist1->Fill(temp_data.Btag_info[3]);
        temp_hist2->Fill(temp_data.Btag_info[2]);
        temp_MbcVSdeltaE_Btag->SetPoint(temp_MbcVSdeltaE_Btag->GetN(), temp_data.Btag_info[2], temp_data.Btag_info[3]);
    }

    temp_MbcVSdeltaE_Btag->SetMarkerColor(1);
    temp_MbcVSdeltaE_Btag->SetMarkerStyle(20);
    temp_MbcVSdeltaE_Btag->SetMarkerSize(0.3);
    temp_MbcVSdeltaE_Btag->SetTitle("M_{bc} vs #DeltaE of B_{tag};M_{bc} [GeV];#DeltaE [GeV]");

    TCanvas* c_temp = new TCanvas("c","",900,900);
    if(SetLogy == true) c_temp->SetLogy();
    temp_hist1->Draw("Hist"); c_temp->SaveAs(("deltaE_Btag_" +std::to_string(index) + ".png").c_str());
    temp_hist2->Draw("Hist"); c_temp->SaveAs(("Mbc_Btag_" +std::to_string(index) + ".png").c_str());
    if(SetLogy == true) c_temp->SetLogy(0);
    temp_MbcVSdeltaE_Btag->Draw("AP"); temp_MbcVSdeltaE_Btag->GetXaxis()->SetLimits(5.2,5.3); temp_MbcVSdeltaE_Btag->GetHistogram()->SetMaximum(0.5); temp_MbcVSdeltaE_Btag->GetHistogram()->SetMinimum(-0.5);
    c_temp->SaveAs(("MbcVSdeltaE_Btag_" +std::to_string(index) + ".png").c_str());
    delete temp_hist1; delete temp_hist2; delete temp_MbcVSdeltaE_Btag; delete c_temp;
}

void Draw_Mbc_deltaE_Btag_projection(std::queue<Data> TotalData_, int index){
    TH2F* temp_hist1 = new TH2F("MbcVSdeltaE",";;", 500,5.2,5.3,500,-0.5,0.5);
    TH1D* projX;
    TH1D* projY;

    while(!TotalData_.empty()){
        Data temp_data = TotalData_.front();
        TotalData_.pop();
        temp_hist1->Fill(temp_data.Btag_info[2], temp_data.Btag_info[3]);
    }

    projX = temp_hist1->ProjectionX();
    projY = temp_hist1->ProjectionY();

    auto c_temp = new TCanvas("c1", "c1",1000,1000); c_temp->cd();
    gStyle->SetOptStat(0);

    TPad *center_pad = new TPad("center_pad", "center_pad",0.0,0.0,0.6,0.6);
    center_pad->Draw();

    TPad *right_pad = new TPad("right_pad", "right_pad",0.55,0.0,1.0,0.6);
    right_pad->Draw();

    TPad *top_pad = new TPad("top_pad", "top_pad",0.0,0.55,0.6,1.0);
    top_pad->Draw();

    center_pad->cd();
    gStyle->SetPalette(1);
    temp_hist1->Draw("Hist");

    top_pad->cd();
    projX->SetFillColor(33);
    projX->Draw("bar");

    right_pad->cd();
    projY->SetFillColor(33);
    projY->Draw("hbar");
    gPad->RedrawAxis();

    c_temp->SaveAs(("MbcVSdeltaE_Btag_proj_" +std::to_string(index) + ".png").c_str());

    delete temp_hist1; delete projX; delete projY; delete c_temp;
}

void BCS(std::queue<Data>* TotalData_, int index, bool select_highest = true){

    printf("===== BCS =====\n");

    std::queue<Data> new_container;

    while(!TotalData_->empty()){
        std::vector<Data> temp;
        Data initial_data = TotalData_->front();
        int experiment_ = initial_data.event_info[0];
        int run_ = initial_data.event_info[1];
        int event_ = initial_data.event_info[2];
        int ncandidates_ = initial_data.event_info[4];
        while(true){ // I suppose that the order of data exists
            Data temp_data = TotalData_->front();
            if(temp_data.event_info[0] == experiment_ && temp_data.event_info[1] == run_ && temp_data.event_info[2] == event_ && temp_data.event_info[4] == ncandidates_){
                TotalData_->pop();
                temp.push_back(temp_data);
            }
            else break;
        }

        if(select_highest == true){
            double max = -std::numeric_limits<double>::max();
            int best_candidate_index = -1;
            for(unsigned int i = 0; i<temp.size();i++){
                if(temp.at(i).Btag_info[index] > max) { max = temp.at(i).Btag_info[index]; best_candidate_index = i; }
            }
            if(best_candidate_index == -1) { printf("error!\n"); exit(1); }
            new_container.push(temp[best_candidate_index]);
        }
        else if(select_highest == false){
            double min = std::numeric_limits<double>::max();
            int best_candidate_index = -1;
            for(unsigned int i = 0; i<temp.size();i++){
                if(temp.at(i).Btag_info[index] < min) { min = temp.at(i).Btag_info[index]; best_candidate_index = i; }
            }
            if(best_candidate_index == -1) { printf("error!\n"); exit(1); }
            new_container.push(temp[best_candidate_index]);
        }

    }
    while(!TotalData_->empty()) TotalData_->pop();
    while(!new_container.empty()){
        Data temp_data = new_container.front();
        new_container.pop();
        TotalData_->push(temp_data);
    }
}

bool IsBCSValid(std::queue<Data> TotalData_){
    typedef struct labels{
        int __experiment__;
        int __run__;
        int __event__;
        int __ncandidates__;
    } Labels;

    std::vector<Labels> label_list;

    while(!TotalData_.empty()){
        Data temp = TotalData_.front();
        TotalData_.pop();
        for(unsigned int i=0; i<label_list.size(); i++){
            if(label_list.at(i).__experiment__ == temp.event_info[0] && label_list.at(i).__run__ == temp.event_info[1] && label_list.at(i).__event__ == temp.event_info[2] && label_list.at(i).__ncandidates__ == temp.event_info[4]) return false;
        }
        Labels temp_Labels;
        temp_Labels.__experiment__ = temp.event_info[0];
        temp_Labels.__run__ = temp.event_info[1];
        temp_Labels.__event__ = temp.event_info[2];
        temp_Labels.__ncandidates__ = temp.event_info[4];
        label_list.push_back(temp_Labels);
    }
    return true;
}

void PrintInformation(std::queue<Data> TotalData_){
    typedef struct labels{
        int __experiment__;
        int __run__;
        int __event__;
        int __ncandidates__;
    } Labels;
    std::vector<Labels> label_list;

    int N_Btag_isSignal_1 = 0;
    int N_Btag_isSignal_not1 = 0;

    int N_Bsig_isSignal_1 = 0;
    int N_Bsig_isSignal_not1 = 0;

    int N_Upsilon_isSignal_1 = 0;
    int N_Upsilon_isSignal_not1 = 0;

    int N_candidate = TotalData_.size();
    int N_event = 0;

    while(!TotalData_.empty()){
        Data temp = TotalData_.front();
        TotalData_.pop();
        if(temp.Upsilon_info[0] > 0.9 && temp.Upsilon_info[0] < 1.1) N_Upsilon_isSignal_1++;
        else {N_Upsilon_isSignal_not1++;}
        if(temp.Bsig_info[0] > 0.9 && temp.Bsig_info[0] < 1.1) N_Bsig_isSignal_1++;
        else {N_Bsig_isSignal_not1++;}
        if(temp.Btag_info[0] > 0.9 && temp.Btag_info[0] < 1.1) N_Btag_isSignal_1++;
        else {N_Btag_isSignal_not1++;}

        bool overlap = false;
        for(unsigned int i=0; i<label_list.size(); i++){
            if(label_list.at(i).__experiment__ == temp.event_info[0] && label_list.at(i).__run__ == temp.event_info[1] && label_list.at(i).__event__ == temp.event_info[2] && label_list.at(i).__ncandidates__ == temp.event_info[4]) { overlap = true; }
        }
        if(overlap == false){
            N_event++;
            Labels temp_Labels;
            temp_Labels.__experiment__ = temp.event_info[0];
            temp_Labels.__run__ = temp.event_info[1];
            temp_Labels.__event__ = temp.event_info[2];
            temp_Labels.__ncandidates__ = temp.event_info[4];
            label_list.push_back(temp_Labels);
        }

    }
    printf("Number of event: %d\n", N_event);
    printf("Number of candidate: %d\n", N_candidate);
    printf("Number of B_tag isSignal = 1: %d\n", N_Btag_isSignal_1);
    printf("Number of B_tag isSignal != 1: %d\n", N_Btag_isSignal_not1);
    printf("Number of B_sig isSignal = 1: %d\n", N_Bsig_isSignal_1);
    printf("Number of B_sig isSignal != 1: %d\n", N_Bsig_isSignal_not1);
    printf("Number of Upsilon isSignal = 1: %d\n", N_Upsilon_isSignal_1);
    printf("Number of Upsilon isSignal != 1: %d\n", N_Upsilon_isSignal_not1);
}

void ReadRootFiles_r(){

    std::vector<string> names;
    const char* dirname = "/home/jwpark/storage/Ntuple32";

    load_files(dirname, &names);

    std::queue<Data> TotalData;
    Data temp = {0};

    unsigned int totalnum_entry = 0;
    for(unsigned int i = 0; i<names.size(); i++){
        TFile *input_file = new TFile( (dirname+std::string("/")+names.at(i)).c_str(),"read");
        TTree *tree_upsilon = (TTree *) input_file -> Get("Upsilon");
        TTree *tree_Bsig = (TTree *) input_file -> Get("Bsig");
        TTree *tree_Btag = (TTree *) input_file -> Get("Btag");

        // get event_info
        tree_upsilon->SetBranchAddress("__experiment__", &temp.event_info[0]);
        tree_upsilon->SetBranchAddress("__run__", &temp.event_info[1]);
        tree_upsilon->SetBranchAddress("__event__", &temp.event_info[2]);
        tree_upsilon->SetBranchAddress("__candidate__", &temp.event_info[3]);
        tree_upsilon->SetBranchAddress("__ncandidates__", &temp.event_info[4]);
        tree_Bsig->SetBranchAddress("__experiment__", &temp.event_info[5]);
        tree_Bsig->SetBranchAddress("__run__", &temp.event_info[6]);
        tree_Bsig->SetBranchAddress("__event__", &temp.event_info[7]);
        tree_Bsig->SetBranchAddress("__candidate__", &temp.event_info[8]);
        tree_Bsig->SetBranchAddress("__ncandidates__", &temp.event_info[9]);
        tree_Btag->SetBranchAddress("__experiment__", &temp.event_info[10]);
        tree_Btag->SetBranchAddress("__run__", &temp.event_info[11]);
        tree_Btag->SetBranchAddress("__event__", &temp.event_info[12]);
        tree_Btag->SetBranchAddress("__candidate__", &temp.event_info[13]);
        tree_Btag->SetBranchAddress("__ncandidates__", &temp.event_info[14]);

        // get Upsilon_info
        tree_upsilon->SetBranchAddress("isSignal", &temp.Upsilon_info[0]);
        tree_upsilon->SetBranchAddress("nROE_ECLClusters__bocleanMask__bc",&temp.Upsilon_info[1]);
        tree_upsilon->SetBranchAddress("nROE_KLMClusters",&temp.Upsilon_info[2]);
        tree_upsilon->SetBranchAddress("roeE__bocleanMask__bc",&temp.Upsilon_info[3]);
        tree_upsilon->SetBranchAddress("nROE_Tracks__bocleanMask__bc",&temp.Upsilon_info[4]);
        tree_upsilon->SetBranchAddress("roeEextra__bocleanMask__bc",&temp.Upsilon_info[5]);
        tree_upsilon->SetBranchAddress("nROE_NeutralECLClusters__bocleanMask__bc",&temp.Upsilon_info[6]);
        tree_upsilon->SetBranchAddress("roeNeextra__bocleanMask__bc",&temp.Upsilon_info[7]);
        tree_upsilon->SetBranchAddress("useCMSFrame__boroeE__bocleanMask__bc__bc",&temp.Upsilon_info[8]);
        tree_upsilon->SetBranchAddress("useCMSFrame__boroeEextra__bocleanMask__bc__bc",&temp.Upsilon_info[9]);
        tree_upsilon->SetBranchAddress("useCMSFrame__boroeNeextra__bocleanMask__bc__bc",&temp.Upsilon_info[10]);
        tree_upsilon->SetBranchAddress("nROE_ParticlesInList__boK_S0__clgood__bc",&temp.Upsilon_info[11]);
        tree_upsilon->SetBranchAddress("nROE_ParticlesInList__bopi0__clgood__bc",&temp.Upsilon_info[12]);

        // get Bsig_info
        tree_Bsig->SetBranchAddress("Bsig_isSignal", &temp.Bsig_info[0]);
        tree_Bsig->SetBranchAddress("Bsig_E", &temp.Bsig_info[1]);
        tree_Bsig->SetBranchAddress("Bsig_useCMSFrame_E", &temp.Bsig_info[2]);
        tree_upsilon->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp.Bsig_info[3]);
        tree_Bsig->SetBranchAddress("Bsig_extraInfo_decayModeID", &temp.Bsig_info[4]);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_atcPIDBelle_3_2", &temp.Bsig_info[5]);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_mcPDG", &temp.Bsig_info[6]);
        tree_Bsig->SetBranchAddress("Bsig_p", &temp.Bsig_info[7]);
        tree_Bsig->SetBranchAddress("Bsig_useCMSFrame_p", &temp.Bsig_info[8]);
        tree_upsilon->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp.Bsig_info[9]);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_dr", &temp.Bsig_info[10]);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_dz", &temp.Bsig_info[11]);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_eIDBelle", &temp.Bsig_info[12]);
	tree_Bsig->SetBranchAddress("Bsig_daughter_0_muIDBelle", &temp.Bsig_info[13]);

	// get Btag_info
        tree_Btag->SetBranchAddress("Btag_isSignal", &temp.Btag_info[0]);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &temp.Btag_info[1]);
        tree_Btag->SetBranchAddress("Btag_Mbc", &temp.Btag_info[2]);
        tree_Btag->SetBranchAddress("Btag_deltaE", &temp.Btag_info[3]);
        tree_Btag->SetBranchAddress("Btag_E", &temp.Btag_info[4]);
        tree_Btag->SetBranchAddress("Btag_useCMSFrame_E", &temp.Btag_info[5]);
        tree_Btag->SetBranchAddress("Btag_extraInfo_SignalProbability", &temp.Btag_info[6]);

        printf("%s",("Read "+names.at(i) + "... ").c_str());
        printf("%lld entries...\n",tree_upsilon->GetEntries());
        totalnum_entry += tree_upsilon->GetEntries();
        for(unsigned int j = 0; j < tree_upsilon->GetEntries(); j++){ // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            TotalData.push(temp);
        }
    }
    if(event_info_is_valid(TotalData) == false) { printf("error!\n"); return; }
    printf("Total %d entries\n", totalnum_entry);

    PrintInformation(TotalData);

    // draw deltaE and Mbc of Btag
    Draw_Mbc_deltaE_Btag_projection(TotalData, 0);

    // cut Mbc > 5.2
    {
        printf("====== Mbc > 5.2 GeV =====\n");
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Btag_info[2] > 5.2) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }

    PrintInformation(TotalData);

    Draw_Mbc_deltaE_Btag_projection(TotalData, 1);

    // cut abs(deltaE) < 0.5
    {
        printf("===== abs(deltaE) < 0.5 =====\n");
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Btag_info[3] > -0.5 && temp_data.Btag_info[3] < 0.5) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }
    PrintInformation(TotalData);

    {
        TH1F* temp_hist1 = new TH1F("SignalProbability_Btag","SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate",100,-10,0);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist1->Fill(TMath::Log10(temp_data.Btag_info[6]));
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist1->Draw("Hist"); c_temp->SaveAs("SignalProbability_distribution_after_loose_deltaE_cut.png");
        delete temp_hist1; delete c_temp;
    }

    // cut SignalProbability > 0.01
    {
        printf("===== SignalProbability > 0.01 =====\n");
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Btag_info[6] > 0.01) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }
    PrintInformation(TotalData);

    Draw_Mbc_deltaE_Btag_projection(TotalData, 2);

    // draw atcPIDBelle_3_2 distribution with respect to mcPDG
    {
        TH1F* temp_hist1 = new TH1F("atcPID(3,2)_Kaon","atcPID(3,2) of daughter of B_{sig};atcPID(3,2);Num of candidate",100,-0.1,1.1);
        TH1F* temp_hist2 = new TH1F("atcPID(3,2)_non-Kaon","atcPID(3,2) of daughter of B_{sig};atcPID(3,2);Num of candidate",100,-0.1,1.1);
        TH1F* temp_hist3 = new TH1F("atcPID(3,2)","atcPID(3,2) of daughter of B_{sig};atcPID(3,2);Num of candidate",100,-0.1,1.1);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Bsig_info[6] > 320.5 && temp_data.Bsig_info[6] < 321.5) temp_hist1->Fill(temp_data.Bsig_info[5]);
            else if(temp_data.Bsig_info[6] > -321.5 && temp_data.Bsig_info[6] < -320.5) temp_hist1->Fill(temp_data.Bsig_info[5]);
            else { temp_hist2->Fill(temp_data.Bsig_info[5]); }
            temp_hist3->Fill(temp_data.Bsig_info[5]);
        }
        temp_hist1->SetStats(false);
        temp_hist1->SetLineColor(kBlue);
        temp_hist1->SetLineWidth(2);
        temp_hist1->SetFillColor(kBlue);
        temp_hist1->SetFillStyle(3013);
        temp_hist2->SetStats(false);
        temp_hist2->SetLineColor(kRed);
        temp_hist2->SetLineWidth(2);
        temp_hist2->SetFillColor(kRed);
        temp_hist2->SetFillStyle(3007);
        
        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist2->Draw("Hist"); temp_hist1->Draw("SAME"); c_temp->SetLogy(); c_temp->SaveAs("atcPID_distribution.png");
        temp_hist3->Draw("Hist"); c_temp->SaveAs("atcPID_distribution_total.png");
        delete temp_hist1; delete temp_hist2; delete temp_hist3; delete c_temp;
    }

    // cut: atcPID(3,2) > 0.6 for daughter of Bsig
    {
        printf("===== atcPID(3,2) of Bsig > 0.6 =====\n");
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Bsig_info[5] > 0.6) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }

    PrintInformation(TotalData);

    Draw_Mbc_deltaE_Btag_projection(TotalData, 3);

    { // print dr of Kaon
        TH1F* temp_hist = new TH1F("dr_Kaon_from_Bsig","dr of Kaon from B_{sig};dr [cm];candidates",100,-0.1,3.5);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Bsig_info[10]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("dr_Kaon_after_actPID_cut.png");
        delete temp_hist; delete c_temp;
    }

    // dr < 2 cm cut
    {
        printf("====== dr < 2 cm =====\n");
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Bsig_info[10] < 2) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }

    PrintInformation(TotalData);


    { // print dz of Kaon
        TH1F* temp_hist = new TH1F("dz_Kaon_from_Bsig","dz of Kaon from B_{sig};dz [cm];candidates",100,-6,6);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Bsig_info[11]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("dz_Kaon_after_dr_cut.png");
        delete temp_hist; delete c_temp;
    }

    // abs(dz) < 4 cm cut
    {
        printf("====== abs(dz) < 4 cm =====\n");
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Bsig_info[11] < 4 && temp_data.Bsig_info[11] > -4) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }

    PrintInformation(TotalData);

    { // print eIDBelle of Kaon
        TH1F* temp_hist = new TH1F("eIDBelle_Bsig_first_daughter","eIDBelle of Kaon from B_{sig};eIDBelle;candidates",100,0,1);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Bsig_info[12]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("eIDBelle_Kaon_after_dz_cut.png");
        delete temp_hist; delete c_temp;
    }

    // eIDBelle < 0.9 cm cut
    {
        printf("====== eIDBelle < 0.9 =====\n");
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Bsig_info[12] < 0.9) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }

    PrintInformation(TotalData);

    { // print muIDBelle of Kaon
        TH1F* temp_hist = new TH1F("muIDBelle_Bsig_first_daughter","muIDBelle of Kaon from B_{sig};muIDBelle;candidates",100,0,1);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Bsig_info[13]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("muIDBelle_Kaon_after_dz_cut.png");
        delete temp_hist; delete c_temp;
    }

    // muIDBelle < 0.9 cm cut
    {
        printf("====== muIDBelle < 0.9 =====\n");
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Bsig_info[13] < 0.9) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }

    PrintInformation(TotalData);

    { // print energy in ROE(cleanMask)
        TH1F* temp_hist = new TH1F("ROE_E_Upsilon","Energy in ROE of #Upsilon(4S) at LAB;energy [GeV];candidates",100,-0.1, 8);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Upsilon_info[3]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("ROE_E_after_muIDBelle_cut.png");
        delete temp_hist; delete c_temp;
    }

    { // print energy in ROE(cleanMask) at CMS
        TH1F* temp_hist = new TH1F("ROE_E_Upsilon_CMS","Energy in ROE of #Upsilon(4S) at CMS;energy [GeV];candidates",100,-0.1, 8);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Upsilon_info[8]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("ROE_E_CMS_after_muIDBelle_cut.png");
        delete temp_hist; delete c_temp;
    }

    { // print energy in ECLC clusters in ROE at CMS (cleanMask)
        TH1F* temp_hist = new TH1F("ROE_ECLC_Upsilon","Energy in ECLClusters in ROE of #Upsilon(4S) at CMS;energy [GeV];candidates",100,-0.1, 8);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Upsilon_info[9]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("ROE_ECLC_after_muIDBelle_cut.png");
        delete temp_hist; delete c_temp;
    }

    { // print energy in ECLC clusters in ROE at CMS (cleanMask)
        TH1F* temp_hist = new TH1F("ROE_NECLC_Upsilon","Energy in neutral ECLClusters in ROE of #Upsilon(4S) at CMS;energy [GeV];candidates",100,-0.1, 8);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Upsilon_info[10]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("ROE_NECLC_after_muIDBelle_cut.png");
        delete temp_hist; delete c_temp;
    }

    { // E_ROE < 1 cut
        printf("====== E_ROE < 1 GeV =====\n");
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Upsilon_info[3] < 1) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }

    PrintInformation(TotalData);

    { // print ntrack
        TH1F* temp_hist = new TH1F("nROE_track_Upsilon","number of tracks in ROE of #Upsilon(4S);number of tracks;evt",100,-0.5,13.5);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Upsilon_info[4]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("nROE_tracks_distribution_after_E_ROE_cut.png");
        delete temp_hist; delete c_temp;
    }

    { // ntrack = 0 cut
        printf("====== ntrack = 0 =====\n");
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Upsilon_info[4] < 0.5) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }

    PrintInformation(TotalData);

    { 
        TH1F* temp_hist1 = new TH1F("SignalProbability_Btag_correct","SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate",100,-10,0);
        TH1F* temp_hist2 = new TH1F("SignalProbability_Btag_non-correct","SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate",100,-10,0);
        TH1F* temp_hist3 = new TH1F("SignalProbability_Btag","SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate",100,-10,0);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Btag_info[0] > 0.9 && temp_data.Btag_info[0] < 1.1) temp_hist1->Fill(TMath::Log10(temp_data.Btag_info[6]));
            else { temp_hist2->Fill(TMath::Log10(temp_data.Btag_info[6])); }
            temp_hist3->Fill(TMath::Log10(temp_data.Btag_info[6]));
        }
        temp_hist1->SetStats(false);
        temp_hist1->SetLineColor(kBlue);
        temp_hist1->SetLineWidth(2);
        temp_hist1->SetFillColor(kBlue);
        temp_hist1->SetFillStyle(3013);
        temp_hist2->SetStats(false);
        temp_hist2->SetLineColor(kRed);
        temp_hist2->SetLineWidth(2);
        temp_hist2->SetFillColor(kRed);
        temp_hist2->SetFillStyle(3007);

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist2->Draw("Hist"); temp_hist1->Draw("SAME"); c_temp->SaveAs("SignalProbability_distribution_before_BCS.png");
        temp_hist3->Draw("Hist"); c_temp->SaveAs("SignalProbability_distribution_total_before_BCS.png");
        delete temp_hist1; delete temp_hist2; delete temp_hist3; delete c_temp;
    }

    // BCS
    BCS(&TotalData, 6, true);

    {
        TH1F* temp_hist1 = new TH1F("SignalProbability_Btag_correct","SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate",100,-10,0);
        TH1F* temp_hist2 = new TH1F("SignalProbability_Btag_non-correct","SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate",100,-10,0);
        TH1F* temp_hist3 = new TH1F("SignalProbability_Btag","SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate",100,-10,0);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Btag_info[0] > 0.9 && temp_data.Btag_info[0] < 1.1) temp_hist1->Fill(TMath::Log10(temp_data.Btag_info[6]));
            else { temp_hist2->Fill(TMath::Log10(temp_data.Btag_info[6])); }
            temp_hist3->Fill(TMath::Log10(temp_data.Btag_info[6]));
        }
        temp_hist1->SetStats(false);
        temp_hist1->SetLineColor(kBlue);
        temp_hist1->SetLineWidth(2);
        temp_hist1->SetFillColor(kBlue);
        temp_hist1->SetFillStyle(3013);
        temp_hist2->SetStats(false);
        temp_hist2->SetLineColor(kRed);
        temp_hist2->SetLineWidth(2);
        temp_hist2->SetFillColor(kRed);
        temp_hist2->SetFillStyle(3007);

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist2->Draw("Hist"); temp_hist1->Draw("SAME"); c_temp->SaveAs("SignalProbability_distribution_after_BCS.png");
        temp_hist3->Draw("Hist"); c_temp->SaveAs("SignalProbability_distribution_total_after_BCS.png");
        delete temp_hist1; delete temp_hist2; delete temp_hist3; delete c_temp;
    }

    Draw_Mbc_deltaE_Btag_projection(TotalData, 4);
    PrintInformation(TotalData);
    if(IsBCSValid(TotalData) == false) {printf("error!\n"); return;}

    // Mbc > 5.27 cut
    {
        printf("====== Mbc > 5.27 GeV =====\n");
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Btag_info[2] > 5.27) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }

    PrintInformation(TotalData);

    Draw_Mbc_deltaE_Btag_projection(TotalData, 5);

    // abs(deltaE) < 0.1 cut
    {
        printf("===== abs(deltaE) < 0.1 =====\n");
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Btag_info[3] > -0.1 && temp_data.Btag_info[3] < 0.1) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }
    PrintInformation(TotalData);

    Draw_Mbc_deltaE_Btag_projection(TotalData, 6);

    {
        TH1F* temp_hist = new TH1F("nROE_ECLcluster_Upsilon","number of ECL clusters in ROE of #Upsilon(4S);number of ECL clusters;evt",14,-0.5,13.5);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Upsilon_info[1]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("nROE_ECLcluster_distribution_after_BCS.png");
        delete temp_hist; delete c_temp;
    }

    {
        TH1F* temp_hist = new TH1F("nROE_KLMcluster_Upsilon","number of KLM clusters in ROE of #Upsilon(4S);number of KLM clusters;evt",14,-0.5,13.5);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Upsilon_info[2]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("nROE_KLMcluster_distribution_after_BCS.png");
        delete temp_hist; delete c_temp;
    }

    {
        TH1F* temp_hist = new TH1F("nROE_energy_Upsilon","energy of ROE of #Upsilon(4S);energy of ROE [GeV];evt",50,0,3);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Upsilon_info[3]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("ROE_energy_distribution_after_BCS.png");
        delete temp_hist; delete c_temp;
    }

    {
        TH1F* temp_hist = new TH1F("nROE_track_Upsilon","number of tracks in ROE of #Upsilon(4S);number of tracks;evt",14,-0.5,13.5);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Upsilon_info[4]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("nROE_tracks_distribution_after_BCS.png");
        delete temp_hist; delete c_temp;
    }

    {
        TH1F* temp_hist = new TH1F("ROE_Eextra","ROE Eextra of #Upsilon(4S);number of tracks;evt",100,-0.5,5.5);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Upsilon_info[5]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("ROE_Eextra_after_BCS.png");
        delete temp_hist; delete c_temp;
    }

    { // print energy in ECLC clusters in ROE at CMS (cleanMask)
        TH1F* temp_hist = new TH1F("ROE_NECLC_Upsilon","Energy in neutral ECLClusters in ROE of #Upsilon(4S) at CMS;energy [GeV];candidates",100,-0.1, 8);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Upsilon_info[10]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("ROE_CECLC_after_BCS.png");
        delete temp_hist; delete c_temp;
    }

    {
        TH1F* temp_hist = new TH1F("Bsig_p_LAB","momentum of B_{sig} at LAB frame;p [GeV];evt",50,-0.5,6);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Bsig_info[7]);
        }

        TCanvas* c_temp = new TCanvas("c","",900,900);
        temp_hist->Draw("Hist"); c_temp->SaveAs("Bsig_p_LAB_after_BCS.png");
        delete temp_hist; delete c_temp;
    }

    {
        TH1F* temp_hist = new TH1F("Bsig_p_CMS","momentum of B_{sig} at CMS frame;p [GeV];evt",50,-0.5,6);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Bsig_info[8]);
        }

        TCanvas* c_temp = new TCanvas("c","",900,900);
        temp_hist->Draw("Hist"); c_temp->SaveAs("Bsig_p_CMS_after_BCS.png");
        delete temp_hist; delete c_temp;
    }

    {
        TH1F* temp_hist = new TH1F("Bsig_p_RecoilRest","momentum of B_{sig} at rest frame of recoil system;p [GeV];evt",50,-0.5,6);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Bsig_info[9]);
        }

        TCanvas* c_temp = new TCanvas("c","",900,900);
        temp_hist->Draw("Hist"); c_temp->SaveAs("Bsig_p_recoil_after_BCS.png");
        delete temp_hist; delete c_temp;
    }

    {
        TH1F* temp_hist = new TH1F("Btag_dmID","decay ID of B_{tag};decay ID;evt",74,-0.5,36.5);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Btag_info[1]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("Btag_dmID_after_BCS.png");
        delete temp_hist; delete c_temp;
    }

    {
        TH1F* temp_hist = new TH1F("nROE_K_S0","number of K_S0:good candidates in ROE of #Upsilon(4S);number of good K_{S}^{0} candidates in ROE;evt",100,-0.5,5.5);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Upsilon_info[11]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("nROE_K_S0_after_BCS.png");
        delete temp_hist; delete c_temp;
    }

    {
        TH1F* temp_hist = new TH1F("nROE_pi0","number of #pi^{0} candidates in ROE of #Upsilon(4S);number of #pi^{0} candidates in ROE;evt",100,-0.5,8.5);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Upsilon_info[12]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("nROE_pi0_after_BCS.png");
        delete temp_hist; delete c_temp;
    }


}
