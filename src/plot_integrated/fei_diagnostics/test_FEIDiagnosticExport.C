// ROOT smoke test independent of the Belle II calibration headers and data.
// root -l -b -q 'test_FEIDiagnosticExport.C("fei_diagnostics_smoke")'
#include "FEIDiagnosticExport.h"
#include <TMemFile.h>

void test_FEIDiagnosticExport(const char* directory="fei_diagnostics_smoke") {
    feidiag::start(directory,"synthetic-smoke-test");
    TH1D mcC("smoke_mcC","",3,0,3), dataC("smoke_dataC","",3,0,3);
    TH1D mcN("smoke_mcN","",1,0,1), dataN("smoke_dataN","",1,0,1);
    for (TH1D* h : {&mcC,&dataC,&mcN,&dataN}) { h->SetDirectory(nullptr); h->Sumw2(); }
    const std::vector<std::string> samples={"DATA","SIGNAL","CHG","MIX","UUBAR","DDBAR","SSBAR","CHARM"};
    for (int si=0;si<8;++si) {
        std::string sample=samples[si];
        TMemFile input((sample+".root").c_str(),"RECREATE");
        TTree u("Upsilon",""), s("Bsig",""), t("Btag",""), x("Xs","");
        double charge=0,mode=0,sig=0,mass=0.6,prob=.02,mva=.8,truth=0;
        double eu=0,es=0,et=0,ex=0;
        u.Branch("charge",&charge); u.Branch("MVA_BB",&mva); u.Branch("__event__",&eu);
        s.Branch("sig",&sig); s.Branch("mass",&mass); s.Branch("__event__",&es);
        t.Branch("mode",&mode); t.Branch("Btag_extraInfo_SignalProbability",&prob); t.Branch("__event__",&et);
        if (sample=="SIGNAL" || sample=="CHG" || sample=="MIX") t.Branch("Btag_isSignal",&truth);
        x.Branch("__event__",&ex);
        for (int c=0;c<4;++c) {
            charge=c==0 ? 1 : 0; mode=c==0 ? 0 : (c==1 ? 3 : (c==2 ? 18 : 31));
            sig=c%3; eu=es=et=ex=c; truth=sample=="SIGNAL" ? 1 : 0;
            u.Fill(); s.Fill(); t.Fill(); x.Fill();
        }
        bool missingRejected=false;
        double dummy=0;
        try { feidiag::bind(&u,"missing_required",&dummy); }
        catch (const std::runtime_error&) { missingRejected=true; }
        feidiag::require(missingRejected,"missing required branch was not rejected");
        feidiag::bind(&u,"missing_optional",&dummy,false);
        feidiag::require(std::isnan(feidiag::scalar(&u,"missing_optional")),"missing optional value was fabricated");
        const double norm=si==0 ? 1 : 2;
        feidiag::output().file(&input,sample+".root",sample,&u,&s,&t,sample=="SIGNAL"?&x:nullptr,norm);
        for (int c=0;c<4;++c) {
            feidiag::read(&u,c); feidiag::read(&s,c); feidiag::read(&t,c); feidiag::read(&x,c);
            feidiag::output().input(sample,charge,mode,sig,mass,norm,true,0);
            std::vector<double> factors(15,1.0); factors[0]=norm; factors[2]=si==0?1:0.5;
            const double official=(sample=="SIGNAL" || sample=="CHG" || sample=="MIX")?0.8:1.0;
            feidiag::output().selected(c,sample,charge,mode,sig,mass,1.0,factors,official);
            if (si==0) { if(c) dataC.Fill(c-.5); else dataN.Fill(.5); }
            else { if(c) mcC.Fill(c-.5); else mcN.Fill(.5); }
        }
        // A common event identifier must be checked, not just entry counts.
        et=-999;
        bool mismatchRejected=false;
        try { feidiag::output().sync(); }
        catch (const std::runtime_error&) { mismatchRejected=true; }
        feidiag::require(mismatchRejected,"event misalignment was not rejected");
    }
    feidiag::output().finish(&mcC,&dataC,&mcN,&dataN);
    std::cout << "ROOT exporter smoke test PASSED" << std::endl;
}
