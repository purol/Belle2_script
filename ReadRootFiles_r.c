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

    double Upsilon_info[5];
    // 0: Upsilon_isSignal; 1: number of ECL clusters in ROE(cleanMask), 2: number of KLM clusters in ROE
    // 3: energy in ROE(cleanMask), 4: number of tracks in ROE(cleanMask)

    double Bsig_info[7];
    // 0: Bsig_isSignal, 1: Bsig_E, 2: Bsig_E_CMS, 3: Bsig_E_Recoil, 4: Bsig_dmID
    // 5: Bsig_first_daughter's_actPID(3,2) 6: Bsig_first_daughter's_mcPDG

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
    temp_c->SaveAs( ("Bsig_" + std::to_string(index) + ".png").c_str() );
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
    temp_c->SaveAs( ("Btag_" + std::to_string(index) + ".png").c_str() );
    delete temp_hist;
    delete temp_c;

}

void BCS(std::queue<Data>* TotalData_, int index, bool select_highest = true){

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

    int N_Btag_isSignal_1 = 0;
    int N_Btag_isSignal_not1 = 0;

    int N_Bsig_isSignal_1 = 0;
    int N_Bsig_isSignal_not1 = 0;

    int N_Upsilon_isSignal_1 = 0;
    int N_Upsilon_isSignal_not1 = 0;

    while(!TotalData_.empty()){
        Data temp = TotalData_.front();
        TotalData_.pop();
        if(temp.Upsilon_info[0] > 0.9 && temp.Upsilon_info[0] < 1.1) N_Upsilon_isSignal_1++;
        else {N_Upsilon_isSignal_not1++;}
        if(temp.Bsig_info[0] > 0.9 && temp.Bsig_info[0] < 1.1) N_Bsig_isSignal_1++;
        else {N_Bsig_isSignal_not1++;}
        if(temp.Btag_info[0] > 0.9 && temp.Btag_info[0] < 1.1) N_Btag_isSignal_1++;
        else {N_Btag_isSignal_not1++;}
    }
    printf("Number of B_tag isSignal = 1: %d\n", N_Btag_isSignal_1);
    printf("Number of B_tag isSignal != 1: %d\n", N_Btag_isSignal_not1);
    printf("Number of B_sig isSignal = 1: %d\n", N_Bsig_isSignal_1);
    printf("Number of B_sig isSignal != 1: %d\n", N_Bsig_isSignal_not1);
    printf("Number of Upsilon isSignal = 1: %d\n", N_Upsilon_isSignal_1);
    printf("Number of Upsilon isSignal != 1: %d\n", N_Upsilon_isSignal_not1);
}

void ReadRootFiles_r(){

    std::vector<string> names;
    const char* dirname = "/home/jwpark/Ntuple20";

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

        // get Bsig_info
        tree_Bsig->SetBranchAddress("Bsig_isSignal", &temp.Bsig_info[0]);
        tree_Bsig->SetBranchAddress("Bsig_E", &temp.Bsig_info[1]);
        tree_Bsig->SetBranchAddress("Bsig_useCMSFrame_E", &temp.Bsig_info[2]);
        tree_upsilon->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp.Bsig_info[3]);
        tree_Bsig->SetBranchAddress("Bsig_extraInfo_decayModeID", &temp.Bsig_info[4]);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_atcPIDBelle_3_2", &temp.Bsig_info[5]);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_mcPDG", &temp.Bsig_info[6]);

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

    for(int i=0; i<5; i++){
        Draw_Bsig_info(TotalData, i);
    }
    for(int i=0; i<6; i++){
        Draw_Btag_info(TotalData, i);
    }

    // draw deltaE and Mbc of Btag
    {
        TH1F* temp_hist1 = new TH1F("deltaE_total","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-1.0,1.0);
        TH1F* temp_hist2 = new TH1F("deltaE_Btag_isSig=1","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-1.0,1.0);
        TH1F* temp_hist3 = new TH1F("deltaE_Btag_isSig!=1","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-1.0,1.0);
        TH1F* temp_hist4 = new TH1F("Mbc_total","Mbc of B_{tag};M_{bc} [GeV];Num of candidate",100,5.1,5.3);
        TH1F* temp_hist5 = new TH1F("Mbc_Btag_isSig=1","M_{bc} of B_{tag};M_{bc} [GeV];Num of candidate",100,5.1,5.4);
        TH1F* temp_hist6 = new TH1F("Mbc_Btag_isSig!=1","M_{bc} of B_{tag};M_{bc} [GeV];Num of candidate",100,5.1,5.4);
        TGraph* temp_MbcVSdeltaE_Btag_isSig_one = new TGraph();
        TGraph* temp_MbcVSdeltaE_Btag_isSig_isnot_one = new TGraph();

        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Btag_info[0] > 0.9 && temp_data.Btag_info[0] < 1.1){
                temp_hist1->Fill(temp_data.Btag_info[3]);
                temp_hist2->Fill(temp_data.Btag_info[3]);
                temp_hist4->Fill(temp_data.Btag_info[2]);
                temp_hist5->Fill(temp_data.Btag_info[2]);
                temp_MbcVSdeltaE_Btag_isSig_one->SetPoint(temp_MbcVSdeltaE_Btag_isSig_one->GetN(),temp_data.Btag_info[2], temp_data.Btag_info[3]);
            }
            else {
                temp_hist1->Fill(temp_data.Btag_info[3]);
                temp_hist3->Fill(temp_data.Btag_info[3]);
                temp_hist4->Fill(temp_data.Btag_info[2]);
                temp_hist6->Fill(temp_data.Btag_info[2]);
                temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetPoint(temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetN(),temp_data.Btag_info[2], temp_data.Btag_info[3]);
            }
        }
        temp_hist2->SetStats(false);
        temp_hist2->SetLineColor(kBlue);
        temp_hist2->SetLineWidth(2);
        temp_hist2->SetFillColor(kBlue);
        temp_hist2->SetFillStyle(3013);
        temp_hist3->SetStats(false);
        temp_hist3->SetLineColor(kRed);
        temp_hist3->SetLineWidth(2);
        temp_hist3->SetFillColor(kRed);
        temp_hist3->SetFillStyle(3007);
        temp_hist5->SetStats(false);
        temp_hist5->SetLineColor(kBlue);
        temp_hist5->SetLineWidth(2);
        temp_hist5->SetFillColor(kBlue);
        temp_hist5->SetFillStyle(3013);
        temp_hist6->SetStats(false);
        temp_hist6->SetLineColor(kRed);
        temp_hist6->SetLineWidth(2);
        temp_hist6->SetFillColor(kRed);
        temp_hist6->SetFillStyle(3007);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerColor(1);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerStyle(20);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerSize(0.4);
        temp_MbcVSdeltaE_Btag_isSig_one->SetTitle("M_{bc} vs #DeltaE of B_{tag};M_{bc} [GeV];#DeltaE [GeV]");
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerColor(2);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerStyle(20);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerSize(0.4);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetTitle("M_{bc} vs #DeltaE of B_{tag};M_{bc} [GeV];#DeltaE [GeV]");

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist1->Draw("Hist"); c_temp->SaveAs("deltaE_total.png");
        temp_hist3->Draw("Hist"); temp_hist2->Draw("SAME"); c_temp->SetLogy(); c_temp->SaveAs("deltaE_total_classified.png"); c_temp->SetLogy(0);
        temp_hist4->Draw("Hist"); c_temp->SaveAs("Mbc_total.png");
        temp_hist6->Draw("Hist"); temp_hist5->Draw("SAME"); c_temp->SetLogy(); c_temp->SaveAs("Mbc_total_classified.png"); c_temp->SetLogy(0);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->Draw("AP"); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetXaxis()->SetLimits(5.1,5.4); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetHistogram()->SetMaximum(1.0); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetHistogram()->SetMinimum(-1.0);
        temp_MbcVSdeltaE_Btag_isSig_one->Draw("P");
        c_temp->SaveAs("MbcVSdeltaE.png");
        delete temp_hist1; delete temp_hist2; delete temp_hist3; delete temp_hist4; delete temp_hist5; delete temp_hist6; delete temp_MbcVSdeltaE_Btag_isSig_one; delete temp_MbcVSdeltaE_Btag_isSig_isnot_one; delete c_temp;
    }
    printf("1. number of candidate: %d\n",TotalData.size());

    // cut Mbc > 5.27
    {
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Btag_info[2] > 5.27) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }
    printf("2. number of candidate: %d\n",TotalData.size());

    {
        TH1F* temp_hist1 = new TH1F("deltaE_total","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-1.0,1.0);
        TH1F* temp_hist2 = new TH1F("deltaE_Btag_isSig=1","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-1.0,1.0);
        TH1F* temp_hist3 = new TH1F("deltaE_Btag_isSig!=1","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-1.0,1.0);
        TH1F* temp_hist4 = new TH1F("Mbc_total","Mbc of B_{tag};M_{bc} [GeV];Num of candidate",100,5.1,5.3);
        TH1F* temp_hist5 = new TH1F("Mbc_Btag_isSig=1","M_{bc} of B_{tag};M_{bc} [GeV];Num of candidate",100,5.1,5.4);
        TH1F* temp_hist6 = new TH1F("Mbc_Btag_isSig!=1","M_{bc} of B_{tag};M_{bc} [GeV];Num of candidate",100,5.1,5.4);
        TGraph* temp_MbcVSdeltaE_Btag_isSig_one = new TGraph();
        TGraph* temp_MbcVSdeltaE_Btag_isSig_isnot_one = new TGraph();

        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Btag_info[0] > 0.9 && temp_data.Btag_info[0] < 1.1){
                temp_hist1->Fill(temp_data.Btag_info[3]);
                temp_hist2->Fill(temp_data.Btag_info[3]);
                temp_hist4->Fill(temp_data.Btag_info[2]);
                temp_hist5->Fill(temp_data.Btag_info[2]);
                temp_MbcVSdeltaE_Btag_isSig_one->SetPoint(temp_MbcVSdeltaE_Btag_isSig_one->GetN(),temp_data.Btag_info[2], temp_data.Btag_info[3]);
            }
            else {
                temp_hist1->Fill(temp_data.Btag_info[3]);
                temp_hist3->Fill(temp_data.Btag_info[3]);
                temp_hist4->Fill(temp_data.Btag_info[2]);
                temp_hist6->Fill(temp_data.Btag_info[2]);
                temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetPoint(temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetN(),temp_data.Btag_info[2], temp_data.Btag_info[3]);
            }
        }
        temp_hist2->SetStats(false);
        temp_hist2->SetLineColor(kBlue);
        temp_hist2->SetLineWidth(2);
        temp_hist2->SetFillColor(kBlue);
        temp_hist2->SetFillStyle(3013);
        temp_hist3->SetStats(false);
        temp_hist3->SetLineColor(kRed);
        temp_hist3->SetLineWidth(2);
        temp_hist3->SetFillColor(kRed);
        temp_hist3->SetFillStyle(3007);
        temp_hist5->SetStats(false);
        temp_hist5->SetLineColor(kBlue);
        temp_hist5->SetLineWidth(2);
        temp_hist5->SetFillColor(kBlue);
        temp_hist5->SetFillStyle(3013);
        temp_hist6->SetStats(false);
        temp_hist6->SetLineColor(kRed);
        temp_hist6->SetLineWidth(2);
        temp_hist6->SetFillColor(kRed);
        temp_hist6->SetFillStyle(3007);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerColor(1);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerStyle(20);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerSize(0.4);
        temp_MbcVSdeltaE_Btag_isSig_one->SetTitle("M_{bc} vs #DeltaE of B_{tag};M_{bc} [GeV];#DeltaE [GeV]");
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerColor(2);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerStyle(20);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerSize(0.4);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetTitle("M_{bc} vs #DeltaE of B_{tag};M_{bc} [GeV];#DeltaE [GeV]");

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist1->Draw("Hist"); c_temp->SaveAs("deltaE_total_afterMbc_cut.png");
        temp_hist3->Draw("Hist"); temp_hist2->Draw("SAME"); c_temp->SetLogy(); c_temp->SaveAs("deltaE_total_classified_afterMbc_cut.png"); c_temp->SetLogy(0);
        temp_hist4->Draw("Hist"); c_temp->SaveAs("Mbc_total_afterMbc_cut.png");
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->Draw("AP"); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetXaxis()->SetLimits(5.1,5.4); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetHistogram()->SetMaximum(1.0); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetHistogram()->SetMinimum(-1.0);
        temp_MbcVSdeltaE_Btag_isSig_one->Draw("P");
        c_temp->SaveAs("MbcVSdeltaE_afterMbc_cut.png");
        delete temp_hist1; delete temp_hist2; delete temp_hist3; delete temp_hist4; delete temp_hist5; delete temp_hist6; delete temp_MbcVSdeltaE_Btag_isSig_one; delete temp_MbcVSdeltaE_Btag_isSig_isnot_one; delete c_temp;
    }
    printf("3. number of candidate: %d\n",TotalData.size());

    // cut abs(deltaE) < 0.1
    {
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Btag_info[3] > -0.1 && temp_data.Btag_info[3] < 0.1) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }
    printf("4. number of candidate: %d\n",TotalData.size());

    {
        TH1F* temp_hist1 = new TH1F("deltaE_total","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-1.0,1.0);
        TH1F* temp_hist2 = new TH1F("deltaE_Btag_isSig=1","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-1.0,1.0);
        TH1F* temp_hist3 = new TH1F("deltaE_Btag_isSig!=1","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-1.0,1.0);
        TH1F* temp_hist4 = new TH1F("Mbc_total","Mbc of B_{tag};M_{bc} [GeV];Num of candidate",100,5.1,5.3);
        TH1F* temp_hist5 = new TH1F("Mbc_Btag_isSig=1","M_{bc} of B_{tag};M_{bc} [GeV];Num of candidate",100,5.1,5.4);
        TH1F* temp_hist6 = new TH1F("Mbc_Btag_isSig!=1","M_{bc} of B_{tag};M_{bc} [GeV];Num of candidate",100,5.1,5.4);
        TGraph* temp_MbcVSdeltaE_Btag_isSig_one = new TGraph();
        TGraph* temp_MbcVSdeltaE_Btag_isSig_isnot_one = new TGraph();

        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Btag_info[0] > 0.9 && temp_data.Btag_info[0] < 1.1){
                temp_hist1->Fill(temp_data.Btag_info[3]);
                temp_hist2->Fill(temp_data.Btag_info[3]);
                temp_hist4->Fill(temp_data.Btag_info[2]);
                temp_hist5->Fill(temp_data.Btag_info[2]);
                temp_MbcVSdeltaE_Btag_isSig_one->SetPoint(temp_MbcVSdeltaE_Btag_isSig_one->GetN(),temp_data.Btag_info[2], temp_data.Btag_info[3]);
            }
            else {
                temp_hist1->Fill(temp_data.Btag_info[3]);
                temp_hist3->Fill(temp_data.Btag_info[3]);
                temp_hist4->Fill(temp_data.Btag_info[2]);
                temp_hist6->Fill(temp_data.Btag_info[2]);
                temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetPoint(temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetN(),temp_data.Btag_info[2], temp_data.Btag_info[3]);
            }
        }
        temp_hist2->SetStats(false);
        temp_hist2->SetLineColor(kBlue);
        temp_hist2->SetLineWidth(2);
        temp_hist2->SetFillColor(kBlue);
        temp_hist2->SetFillStyle(3013);
        temp_hist3->SetStats(false);
        temp_hist3->SetLineColor(kRed);
        temp_hist3->SetLineWidth(2);
        temp_hist3->SetFillColor(kRed);
        temp_hist3->SetFillStyle(3007);
        temp_hist5->SetStats(false);
        temp_hist5->SetLineColor(kBlue);
        temp_hist5->SetLineWidth(2);
        temp_hist5->SetFillColor(kBlue);
        temp_hist5->SetFillStyle(3013);
        temp_hist6->SetStats(false);
        temp_hist6->SetLineColor(kRed);
        temp_hist6->SetLineWidth(2);
        temp_hist6->SetFillColor(kRed);
        temp_hist6->SetFillStyle(3007);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerColor(1);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerStyle(20);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerSize(0.4);
        temp_MbcVSdeltaE_Btag_isSig_one->SetTitle("M_{bc} vs #DeltaE of B_{tag};M_{bc} [GeV];#DeltaE [GeV]");
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerColor(2);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerStyle(20);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerSize(0.4);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetTitle("M_{bc} vs #DeltaE of B_{tag};M_{bc} [GeV];#DeltaE [GeV]");

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist1->Draw("Hist"); c_temp->SaveAs("deltaE_total_afterMbc_cut_after_deltaE_cut.png");
        temp_hist3->Draw("Hist"); temp_hist2->Draw("SAME"); c_temp->SetLogy(); c_temp->SaveAs("deltaE_total_classified_afterMbc_cut_after_deltaE_cut.png"); c_temp->SetLogy(0);
        temp_hist4->Draw("Hist"); c_temp->SaveAs("Mbc_total_afterMbc_cut_after_deltaE_cut.png");
        temp_hist6->Draw("Hist"); temp_hist5->Draw("SAME"); c_temp->SetLogy(); c_temp->SaveAs("Mbc_total_classified_afterMbc_cut_after_deltaE_cut.png"); c_temp->SetLogy(0);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->Draw("AP"); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetXaxis()->SetLimits(5.1,5.4); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetHistogram()->SetMaximum(1.0); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetHistogram()->SetMinimum(-1.0);
        temp_MbcVSdeltaE_Btag_isSig_one->Draw("P");
        c_temp->SaveAs("MbcVSdeltaE_afterMbc_cut_after_deltaE_cut.png");
        delete temp_hist1; delete temp_hist2; delete temp_hist3; delete temp_hist4; delete temp_hist5; delete temp_hist6; delete temp_MbcVSdeltaE_Btag_isSig_one; delete temp_MbcVSdeltaE_Btag_isSig_isnot_one; delete c_temp;
    }
    printf("5. number of candidate: %d\n",TotalData.size());

    // draw atcPIDBelle_3_2 distribution with respect to mcPDG
    {
        TH1F* temp_hist1 = new TH1F("atcPID(3,2)_Kaon","atcPID(3,2) of daughter of B_{sig};atcPID(3,2);Num of candidate",100,-0.1,1.1);
        TH1F* temp_hist2 = new TH1F("atcPID(3,2)_non-Kaon","atcPID(3,2) of daughter of B_{sig};atcPID(3,2);Num of candidate",100,-0.1,1.1);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Bsig_info[6] > 320.5 && temp_data.Bsig_info[6] < 321.5) temp_hist1->Fill(temp_data.Bsig_info[5]);
            else if(temp_data.Bsig_info[6] > -321.5 && temp_data.Bsig_info[6] < -320.5) temp_hist1->Fill(temp_data.Bsig_info[5]);
            else { temp_hist2->Fill(temp_data.Bsig_info[5]); }
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
        delete temp_hist1; delete temp_hist2; delete c_temp;
    }

    // cut: atcPID(3,2) > 0.6 for daughter of Bsig
    {
        std::queue<Data> temp_queue;
        while(!TotalData.empty()){
            Data temp_data = TotalData.front();
            TotalData.pop();
            if(temp_data.Bsig_info[5] > 0.6) temp_queue.push(temp_data);
        }
        TotalData = temp_queue;
    }
    printf("6. number of candidate: %d\n",TotalData.size());

    {
        TH1F* temp_hist1 = new TH1F("deltaE_total","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-0.15,0.15);
        TH1F* temp_hist2 = new TH1F("deltaE_Btag_isSig=1","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-0.15, 0.15);
        TH1F* temp_hist3 = new TH1F("deltaE_Btag_isSig!=1","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-0.15, 0.15);
        TH1F* temp_hist4 = new TH1F("Mbc_total","Mbc of B_{tag};M_{bc} [GeV];Num of candidate",100,5.26,5.30);
        TH1F* temp_hist5 = new TH1F("Mbc_Btag_isSig=1","M_{bc} of B_{tag};M_{bc} [GeV];Num of candidate",100,5.26,5.30);
        TH1F* temp_hist6 = new TH1F("Mbc_Btag_isSig!=1","M_{bc} of B_{tag};M_{bc} [GeV];Num of candidate",100,5.26,5.30);
        TGraph* temp_MbcVSdeltaE_Btag_isSig_one = new TGraph();
        TGraph* temp_MbcVSdeltaE_Btag_isSig_isnot_one = new TGraph();

        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Btag_info[0] > 0.9 && temp_data.Btag_info[0] < 1.1){
                temp_hist1->Fill(temp_data.Btag_info[3]);
                temp_hist2->Fill(temp_data.Btag_info[3]);
                temp_hist4->Fill(temp_data.Btag_info[2]);
                temp_hist5->Fill(temp_data.Btag_info[2]);
                temp_MbcVSdeltaE_Btag_isSig_one->SetPoint(temp_MbcVSdeltaE_Btag_isSig_one->GetN(),temp_data.Btag_info[2], temp_data.Btag_info[3]);
            }
            else {
                temp_hist1->Fill(temp_data.Btag_info[3]);
                temp_hist3->Fill(temp_data.Btag_info[3]);
                temp_hist4->Fill(temp_data.Btag_info[2]);
                temp_hist6->Fill(temp_data.Btag_info[2]);
                temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetPoint(temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetN(),temp_data.Btag_info[2], temp_data.Btag_info[3]);
            }
        }
        temp_hist2->SetStats(false);
        temp_hist2->SetLineColor(kBlue);
        temp_hist2->SetLineWidth(2);
        temp_hist2->SetFillColor(kBlue);
        temp_hist2->SetFillStyle(3013);
        temp_hist3->SetStats(false);
        temp_hist3->SetLineColor(kRed);
        temp_hist3->SetLineWidth(2);
        temp_hist3->SetFillColor(kRed);
        temp_hist3->SetFillStyle(3007);
        temp_hist5->SetStats(false);
        temp_hist5->SetLineColor(kBlue);
        temp_hist5->SetLineWidth(2);
        temp_hist5->SetFillColor(kBlue);
        temp_hist5->SetFillStyle(3013);
        temp_hist6->SetStats(false);
        temp_hist6->SetLineColor(kRed);
        temp_hist6->SetLineWidth(2);
        temp_hist6->SetFillColor(kRed);
        temp_hist6->SetFillStyle(3007);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerColor(1);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerStyle(20);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerSize(0.4);
        temp_MbcVSdeltaE_Btag_isSig_one->SetTitle("M_{bc} vs #DeltaE of B_{tag};M_{bc} [GeV];#DeltaE [GeV]");
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerColor(2);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerStyle(20);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerSize(0.4);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetTitle("M_{bc} vs #DeltaE of B_{tag};M_{bc} [GeV];#DeltaE [GeV]");

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist1->Draw("Hist"); c_temp->SaveAs("deltaE_total_afterMbc_cut_after_deltaE_cut_after_atc_cut.png");
        temp_hist3->Draw("Hist"); temp_hist2->Draw("SAME"); c_temp->SetLogy(); c_temp->SaveAs("deltaE_total_classified_afterMbc_cut_after_deltaE_cut_after_atc_cut.png"); c_temp->SetLogy(0);
        temp_hist4->Draw("Hist"); c_temp->SaveAs("Mbc_total_afterMbc_cut_after_deltaE_cut_after_atc_cut.png");
        temp_hist6->Draw("Hist"); temp_hist5->Draw("SAME"); c_temp->SetLogy(); c_temp->SaveAs("Mbc_total_classified_afterMbc_cut_after_deltaE_cut_after_atc_cut.png"); c_temp->SetLogy(0);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->Draw("AP"); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetXaxis()->SetLimits(5.26,5.30); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetHistogram()->SetMaximum(0.15); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetHistogram()->SetMinimum(-0.15);
        temp_MbcVSdeltaE_Btag_isSig_one->Draw("P");
        c_temp->SaveAs("MbcVSdeltaE_afterMbc_cut_after_deltaE_cut_after_atc_cut.png");
        delete temp_hist1; delete temp_hist2; delete temp_hist3; delete temp_hist4; delete temp_hist5; delete temp_hist6; delete temp_MbcVSdeltaE_Btag_isSig_one; delete temp_MbcVSdeltaE_Btag_isSig_isnot_one; delete c_temp;
    }
    printf("7. number of candidate: %d\n",TotalData.size());

    { 
        TH1F* temp_hist1 = new TH1F("SignalProbability_Btag_correct","SignalProbability of B_{tag};SignalProbability;Num of candidate",100,-0.1,1.1);
        TH1F* temp_hist2 = new TH1F("SignalProbability_Btag_non-correct","SignalProbability of B_{tag};SignalProbability;Num of candidate",100,-0.1,1.1);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Btag_info[0] > 0.9 && temp_data.Btag_info[0] < 1.1) temp_hist1->Fill(temp_data.Btag_info[6]);
            else { temp_hist2->Fill(temp_data.Btag_info[6]); }
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
        temp_hist2->Draw("Hist"); temp_hist1->Draw("SAME"); c_temp->SetLogy(); c_temp->SaveAs("SignalProbability_distribution_before_BCS.png");
        delete temp_hist1; delete temp_hist2; delete c_temp;
    }

    // BCS
    BCS(&TotalData, 6, true);

    {
        TH1F* temp_hist1 = new TH1F("deltaE_total","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-0.15,0.15);
        TH1F* temp_hist2 = new TH1F("deltaE_Btag_isSig=1","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-0.15, 0.15);
        TH1F* temp_hist3 = new TH1F("deltaE_Btag_isSig!=1","#DeltaE of B_{tag};#DeltaE [GeV];Num of candidate",100,-0.15, 0.15);
        TH1F* temp_hist4 = new TH1F("Mbc_total","Mbc of B_{tag};M_{bc} [GeV];Num of candidate",100,5.26,5.30);
        TH1F* temp_hist5 = new TH1F("Mbc_Btag_isSig=1","M_{bc} of B_{tag};M_{bc} [GeV];Num of candidate",100,5.26,5.30);
        TH1F* temp_hist6 = new TH1F("Mbc_Btag_isSig!=1","M_{bc} of B_{tag};M_{bc} [GeV];Num of candidate",100,5.26,5.30);
        TGraph* temp_MbcVSdeltaE_Btag_isSig_one = new TGraph();
        TGraph* temp_MbcVSdeltaE_Btag_isSig_isnot_one = new TGraph();

        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Btag_info[0] > 0.9 && temp_data.Btag_info[0] < 1.1){
                temp_hist1->Fill(temp_data.Btag_info[3]);
                temp_hist2->Fill(temp_data.Btag_info[3]);
                temp_hist4->Fill(temp_data.Btag_info[2]);
                temp_hist5->Fill(temp_data.Btag_info[2]);
                temp_MbcVSdeltaE_Btag_isSig_one->SetPoint(temp_MbcVSdeltaE_Btag_isSig_one->GetN(),temp_data.Btag_info[2], temp_data.Btag_info[3]);
            }
            else {
                temp_hist1->Fill(temp_data.Btag_info[3]);
                temp_hist3->Fill(temp_data.Btag_info[3]);
                temp_hist4->Fill(temp_data.Btag_info[2]);
                temp_hist6->Fill(temp_data.Btag_info[2]);
                temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetPoint(temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetN(),temp_data.Btag_info[2], temp_data.Btag_info[3]);
            }
        }
        temp_hist2->SetStats(false);
        temp_hist2->SetLineColor(kBlue);
        temp_hist2->SetLineWidth(2);
        temp_hist2->SetFillColor(kBlue);
        temp_hist2->SetFillStyle(3013);
        temp_hist3->SetStats(false);
        temp_hist3->SetLineColor(kRed);
        temp_hist3->SetLineWidth(2);
        temp_hist3->SetFillColor(kRed);
        temp_hist3->SetFillStyle(3007);
        temp_hist5->SetStats(false);
        temp_hist5->SetLineColor(kBlue);
        temp_hist5->SetLineWidth(2);
        temp_hist5->SetFillColor(kBlue);
        temp_hist5->SetFillStyle(3013);
        temp_hist6->SetStats(false);
        temp_hist6->SetLineColor(kRed);
        temp_hist6->SetLineWidth(2);
        temp_hist6->SetFillColor(kRed);
        temp_hist6->SetFillStyle(3007);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerColor(1);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerStyle(20);
        temp_MbcVSdeltaE_Btag_isSig_one->SetMarkerSize(0.4);
        temp_MbcVSdeltaE_Btag_isSig_one->SetTitle("M_{bc} vs #DeltaE of B_{tag};M_{bc} [GeV];#DeltaE [GeV]");
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerColor(2);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerStyle(20);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetMarkerSize(0.4);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->SetTitle("M_{bc} vs #DeltaE of B_{tag};M_{bc} [GeV];#DeltaE [GeV]");

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist1->Draw("Hist"); c_temp->SaveAs("deltaE_total_afterMbc_cut_after_deltaE_cut_after_atc_cut_after_BCS.png");
        temp_hist3->Draw("Hist"); temp_hist2->Draw("SAME"); c_temp->SetLogy(); c_temp->SaveAs("deltaE_total_classified_afterMbc_cut_after_deltaE_cut_after_atc_cut_after_BCS.png"); c_temp->SetLogy(0);
        temp_hist4->Draw("Hist"); c_temp->SaveAs("Mbc_total_afterMbc_cut_after_deltaE_cut_after_atc_cut_after_BCS.png");
        temp_hist6->Draw("Hist"); temp_hist5->Draw("SAME"); c_temp->SetLogy(); c_temp->SaveAs("Mbc_total_classified_afterMbc_cut_after_deltaE_cut_after_atc_cut_after_BCS.png"); c_temp->SetLogy(0);
        temp_MbcVSdeltaE_Btag_isSig_isnot_one->Draw("AP"); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetXaxis()->SetLimits(5.26,5.30); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetHistogram()->SetMaximum(0.15); temp_MbcVSdeltaE_Btag_isSig_isnot_one->GetHistogram()->SetMinimum(-0.15);
        temp_MbcVSdeltaE_Btag_isSig_one->Draw("P");
        c_temp->SaveAs("MbcVSdeltaE_afterMbc_cut_after_deltaE_cut_after_atc_cut_after_BCS.png");
        delete temp_hist1; delete temp_hist2; delete temp_hist3; delete temp_hist4; delete temp_hist5; delete temp_hist6; delete temp_MbcVSdeltaE_Btag_isSig_one; delete temp_MbcVSdeltaE_Btag_isSig_isnot_one; delete c_temp;
    }
    printf("8. number of candidate: %d\n",TotalData.size());
    if(IsBCSValid(TotalData) == false) {printf("error!\n"); return;}

    { 
        TH1F* temp_hist1 = new TH1F("SignalProbability_Btag_correct","SignalProbability of B_{tag};SignalProbability;Num of candidate",100,-0.1,1.1);
        TH1F* temp_hist2 = new TH1F("SignalProbability_Btag_non-correct","SignalProbability of B_{tag};SignalProbability;Num of candidate",100,-0.1,1.1);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Btag_info[0] > 0.9 && temp_data.Btag_info[0] < 1.1) temp_hist1->Fill(temp_data.Btag_info[6]);
            else { temp_hist2->Fill(temp_data.Btag_info[6]); }
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
        temp_hist2->Draw("Hist"); temp_hist1->Draw("SAME"); c_temp->SetLogy(); c_temp->SaveAs("SignalProbability_distribution_after_BCS.png");
        delete temp_hist1; delete temp_hist2; delete c_temp;
    }

    {
        TH1F* temp_hist1 = new TH1F("nROE_ECLcluster_Upsilon_correct","number of ECL clusters in ROE of #Upsilon(4S);number of ECL clusters;evt",14,-0.5,13.5);
        TH1F* temp_hist2 = new TH1F("nROE_ECLcluster_Upsilon_non-correct","number of ECL clusters in ROE of #Upsilon(4S);number of ECL clusters;evt",14,-0.5,13.5);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Upsilon_info[0] > 0.9 && temp_data.Upsilon_info[0] < 1.1) temp_hist1->Fill(temp_data.Upsilon_info[1]);
            else { temp_hist2->Fill(temp_data.Upsilon_info[1]); }
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
        temp_hist2->Draw("Hist"); temp_hist1->Draw("SAME"); c_temp->SaveAs("nROE_ECLcluster_distribution_after_BCS.png");
        delete temp_hist1; delete temp_hist2; delete c_temp;
    }

    {
        TH1F* temp_hist1 = new TH1F("nROE_KLMcluster_Upsilon_correct","number of KLM clusters in ROE of #Upsilon(4S);number of KLM clusters;evt",14,-0.5,13.5);
        TH1F* temp_hist2 = new TH1F("nROE_KLMcluster_Upsilon_non-correct","number of KLM clusters in ROE of #Upsilon(4S);number of KLM clusters;evt",14,-0.5,13.5);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Upsilon_info[0] > 0.9 && temp_data.Upsilon_info[0] < 1.1) temp_hist1->Fill(temp_data.Upsilon_info[2]);
            else { temp_hist2->Fill(temp_data.Upsilon_info[2]); }
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
        temp_hist2->Draw("Hist"); temp_hist1->Draw("SAME"); c_temp->SaveAs("nROE_KLMcluster_distribution_after_BCS.png");
        delete temp_hist1; delete temp_hist2; delete c_temp;
    }

    {
        TH1F* temp_hist1 = new TH1F("nROE_energy_Upsilon_correct","energy of ROE of #Upsilon(4S);energy of ROE [GeV];evt",50,0,3);
        TH1F* temp_hist2 = new TH1F("nROE_energy_Upsilon_non-correct","energy of ROE of #Upsilon(4S);energy of ROE;evt",50,0,3);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Upsilon_info[0] > 0.9 && temp_data.Upsilon_info[0] < 1.1) temp_hist1->Fill(temp_data.Upsilon_info[3]);
            else { temp_hist2->Fill(temp_data.Upsilon_info[3]); }
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
        temp_hist2->Draw("Hist"); temp_hist1->Draw("SAME"); c_temp->SaveAs("ROE_energy_distribution_after_BCS.png");
        delete temp_hist1; delete temp_hist2; delete c_temp;
    }

    {
        TH1F* temp_hist1 = new TH1F("nROE_track_Upsilon_correct","number of tracks in ROE of #Upsilon(4S);number of tracks;evt",14,-0.5,13.5);
        TH1F* temp_hist2 = new TH1F("nROE_track_Upsilon_non-correct","number of tracks in ROE of #Upsilon(4S);number of tracks;evt",14,-0.5,13.5);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            if(temp_data.Upsilon_info[0] > 0.9 && temp_data.Upsilon_info[0] < 1.1) temp_hist1->Fill(temp_data.Upsilon_info[4]);
            else { temp_hist2->Fill(temp_data.Upsilon_info[4]); }
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
        temp_hist2->Draw("Hist"); temp_hist1->Draw("SAME"); c_temp->SaveAs("nROE_tracks_distribution_after_BCS.png");
        delete temp_hist1; delete temp_hist2; delete c_temp;
    }

    {
        TH1F* temp_hist = new TH1F("Bsig_E_LAB","energy of B_{sig} at LAB frame;energy [GeV];evt",50,-0.5,6);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Bsig_info[1]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("Bsig_E_LAB_after_BCS.png");
        delete temp_hist; delete c_temp;
    }

    {
        TH1F* temp_hist = new TH1F("Bsig_E_CMS","energy of B_{sig} at CMS frame;energy [GeV];evt",50,-0.5,6);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Bsig_info[2]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("Bsig_E_CMS_after_BCS.png");
        delete temp_hist; delete c_temp;
    }

    {
        TH1F* temp_hist = new TH1F("Bsig_E_RecoilRest","energy of B_{sig} at rest frame of recoil system;energy [GeV];evt",50,-0.5,6);
        std::queue<Data> temp_queue = TotalData;
        while(!temp_queue.empty()){
            Data temp_data = temp_queue.front();
            temp_queue.pop();
            temp_hist->Fill(temp_data.Bsig_info[3]);
        }

        TCanvas* c_temp = new TCanvas("c","",1500,1200);
        temp_hist->Draw("Hist"); c_temp->SaveAs("Bsig_E_recoil_after_BCS.png");
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
    PrintInformation(TotalData);
}
