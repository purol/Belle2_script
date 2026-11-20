#ifndef FEI_DIAGNOSTIC_EXPORT_H
#define FEI_DIAGNOSTIC_EXPORT_H

// CSV diagnostics around the original event loop; no replacement calibration model.
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TSystem.h>
#include <TH1D.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace feidiag {
inline double missing() { return std::numeric_limits<double>::quiet_NaN(); }
inline std::string quote(const std::string& s) {
    std::string r = "\"";
    for (char c : s) { if (c == '"') r += '"'; r += c; }
    return r + '"';
}
inline void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error("[FEI diagnostics] " + message);
}
inline bool integer(double x) { return std::isfinite(x) && std::abs(x-std::round(x)) < 1e-7; }
inline int category(double upsilon, double mode) {
    require(integer(upsilon) && (upsilon == 0 || upsilon == 1), "invalid Upsilon decayModeID");
    require(integer(mode) && mode >= 0, "invalid Btag decayModeID");
    if (upsilon == 1) return 0; // combined neutral tags
    int m = static_cast<int>(mode);
    if (m == 0 || m == 1 || m == 3 || m == 4) return 1;
    if (m == 15 || m == 16 || m == 18 || m == 19) return 2;
    return 3;
}
inline int sampleID(const std::string& s) {
    const std::vector<std::string> names = {"DATA", "SIGNAL", "CHG", "MIX", "UUBAR", "DDBAR", "SSBAR", "CHARM"};
    for (size_t i=0; i<names.size(); ++i) if (s == names[i]) return static_cast<int>(i);
    throw std::runtime_error("Unknown sample " + s);
}
inline double scalar(TTree* t, const char* name) {
    TLeaf* leaf = t ? t->GetLeaf(name) : nullptr;
    if (!leaf || leaf->GetLeafCount() || leaf->GetLenStatic() != 1) return missing();
    return leaf->GetValue();
}
inline bool scalarExists(TTree* t, const char* name) {
    TLeaf* leaf = t ? t->GetLeaf(name) : nullptr;
    return leaf && !leaf->GetLeafCount() && leaf->GetLenStatic() == 1;
}
template<class T> inline void bind(TTree* tree, const char* name, T* address, bool required=true) {
    require(tree != nullptr, "missing tree when binding " + std::string(name));
    if (!tree->GetBranch(name)) {
        require(!required, "missing required branch " + std::string(tree->GetName()) + "/" + name);
        return;
    }
    int result = tree->SetBranchAddress(name, address);
    require(result >= 0, "SetBranchAddress failed (status " + std::to_string(result) + ") for " + std::string(tree->GetName()) + "/" + name);
}
inline void read(TTree* tree, Long64_t entry) {
    require(tree && tree->GetEntry(entry) > 0, "cannot read tree entry " + std::to_string(entry));
}
struct OptionalVariable { const char* output; int tree; const char* branch; };
inline const std::vector<OptionalVariable>& variables() {
    static const std::vector<OptionalVariable> vars = {
        {"tag_is_signal",2,"Btag_isSignal"},
        {"sig_is_signal",1,"Bsig_isSignal"},
        {"tag_mc_errors",2,"Btag_mcErrors"},
        {"tag_mc_pdg",2,"Btag_mcPDG"},
        {"tag_probability",2,"Btag_extraInfo_SignalProbability"},
        {"tag_mbc",2,"Btag_Mbc"},
        {"tag_delta_e",2,"Btag_deltaE"},
        {"mva_bb",0,"MVA_BB"},
        {"xs_mass",1,"Bsig_daughter_0_M"},
        {"tag_wrong_daughters",2,"Btag_extraInfo_looseMCWrongDaughterN"},
        {"tag_missing_gamma_ok",2,"Btag_isSignalAcceptMissingGamma"},
        {"experiment",0,"__experiment__"},
        {"run",0,"__run__"},
        {"event",0,"__event__"},
        {"production",0,"__production__"},
        {"candidate",0,"__candidate__"}
    };
    return vars;
}
inline const std::vector<std::string>& factorNames() {
    static const std::vector<std::string> names = {"normalization", "input_fei", "pi0", "kaon_id", "pion_id", "fake_pid", "ks0", "knn", "xnn", "photon_multiplicity", "kpklkl", "ksklkl", "kstarklkl", "xklkl", "b_to_dkl"};
    return names;
}
struct Moments {
    Long64_t n = 0;
    double sumw = 0, sumw2 = 0;
    void add(double w) { ++n; sumw += w; sumw2 += w*w; }
};
class Exporter {
public:
    std::string directory, currentPath;
    std::ofstream events, files, availability;
    std::map<std::string, Moments> cutflow;
    TTree* trees[3] = {nullptr,nullptr,nullptr};
    TTree* xsTree = nullptr;
    int currentFile = -1;
    Exporter(const std::string& dir, const std::string& sourceHash) : directory(dir) {
        require(!dir.empty(), "empty output directory");
        gSystem->mkdir(dir.c_str(), true);
        require(gSystem->AccessPathName((dir + "/COMPLETE").c_str()), "output already complete; choose a new output directory");
        events.open(dir + "/selected_events.csv");
        files.open(dir + "/files.csv");
        availability.open(dir + "/branch_availability.csv");
        require(events && files && availability, "cannot create output files in " + dir);
        events << std::setprecision(17);
        events << "file_id,entry,sample,category,upsilon_id,tag_mode,sig_mode,selection_mass,weight,official_fei";
        for (const auto& v : variables()) events << ',' << v.output;
        for (const auto& name : factorNames()) events << ",f_" << name;
        events << '\n';
        files << "file_id,sample,path,entries,normalization,sync_keys\n" << std::setprecision(17);
        availability << "file_id,sample,variable,tree,branch,available\n";
        std::ofstream meta(dir + "/run_info.txt");
        meta << "schema_version=1\nsource_sha256=" << sourceHash << "\n"
             << "Selected entries use the original cuts and weight expression.\n"
             << "Cutflow uses normalization only (no detector/background corrections).\n"
             << "Preselection means the entries present in the supplied ntuple, not generated events.\n"
             << "The independently indexed Jpsi tree is deliberately not joined to event trees.\n";
    }
    void file(TFile* input, const std::string& path, const std::string& sample,
              TTree* u, TTree* s, TTree* t, TTree* xs, double norm) {
        require(input && !input->IsZombie(), "cannot open " + path);
        require(u && s && t, "missing Upsilon, Bsig or Btag tree in " + path);
        require(u->GetEntries() == s->GetEntries() && u->GetEntries() == t->GetEntries(), "tree entry counts differ in " + path);
        if (sample == "SIGNAL") require(xs && xs->GetEntries() == u->GetEntries(), "missing or misaligned Xs tree in " + path);
        require(std::isfinite(norm), "invalid normalization in " + path);
        currentPath = path; ++currentFile;
        trees[0]=u; trees[1]=s; trees[2]=t; xsTree=xs;
        std::string sync;
        for (const char* key : {"__experiment__", "__run__", "__event__", "__production__", "__candidate__"}) {
            if (scalarExists(u,key) && scalarExists(s,key) && scalarExists(t,key)) sync += std::string(key) + ';';
        }
        files << currentFile << ',' << sample << ',' << quote(path) << ',' << u->GetEntries() << ',' << norm << ',' << quote(sync) << '\n';
        for (const auto& v : variables()) {
            availability << currentFile << ',' << sample << ',' << v.output << ',' << trees[v.tree]->GetName()
                         << ',' << v.branch << ',' << scalarExists(trees[v.tree],v.branch) << '\n';
        }
    }
    void sync() {
        for (const char* key : {"__experiment__", "__run__", "__event__", "__production__", "__candidate__"}) {
            for (TTree* t : {trees[1],trees[2],xsTree}) {
                if (scalarExists(trees[0],key) && scalarExists(t,key)) {
                    require(scalar(trees[0],key) == scalar(t,key), "tree identifier mismatch for " + std::string(key) + " in " + currentPath);
                }
            }
        }
    }
    void input(const std::string& sample, double u, double tag, double sig, double mass,
               double norm, bool includeKstar, int option) {
        sync();
        int cat=category(u,tag);
        require(integer(sig) && std::isfinite(mass), "invalid signal mode or mass in " + currentPath);
        bool mode = (sig > -0.5 && sig < (includeKstar ? 2.5 : 0.5));
        bool massCut = mode && (!includeKstar || mass < 1.0);
        bool final = massCut && (option == 0 || (option == 1 && u == 0) || (option == 2 && u == 1));
        const bool stages[] = {true,mode,massCut,final};
        for (int stage=0; stage<4; ++stage) if (stages[stage]) {
            std::string key=sample + ',' + std::to_string(cat) + ',' + std::to_string(static_cast<int>(tag)) + ',' + std::to_string(stage);
            cutflow[key].add(norm);
        }
    }
    void selected(Long64_t entry, const std::string& sample, double u, double tag, double sig, double mass,
                  double weight, const std::vector<double>& factors, double official) {
        require(factors.size() == factorNames().size(), "wrong factor count");
        require(std::isfinite(weight) && std::isfinite(official) && official > 0, "non-finite weight or invalid official FEI factor in " + currentPath);
        double product=1;
        for (double f : factors) { require(std::isfinite(f), "non-finite correction in " + currentPath); product *= f; }
        require(std::abs(product-weight) <= 1e-10*std::max(1.0,std::abs(weight)), "factor product disagrees with original event weight");
        events << currentFile << ',' << entry << ',' << sample << ',' << category(u,tag) << ',' << u << ',' << tag << ',' << sig << ',' << mass << ',' << weight << ',' << official;
        for (const auto& v : variables()) events << ',' << scalar(trees[v.tree],v.branch);
        for (double f : factors) events << ',' << f;
        events << '\n';
        require(static_cast<bool>(events), "failed writing selected_events.csv");
    }
    void finish(TH1D* mcCharged, TH1D* dataCharged, TH1D* mcNeutral, TH1D* dataNeutral) {
        std::ofstream cuts(directory + "/input_cutflow.csv");
        cuts << "sample,category,tag_mode,stage,n,sumw,sumw2\n" << std::setprecision(17);
        for (const auto& kv : cutflow) cuts << kv.first << ',' << kv.second.n << ',' << kv.second.sumw << ',' << kv.second.sumw2 << '\n';
        std::ofstream reference(directory + "/reference.csv");
        reference << "category,data,mc,mc_sumw2\n" << std::setprecision(17);
        for (int c=0;c<4;++c) {
            TH1D* m=c ? mcCharged : mcNeutral;
            TH1D* d=c ? dataCharged : dataNeutral;
            int bin=c ? c : 1;
            reference << c << ',' << d->GetBinContent(bin) << ',' << m->GetBinContent(bin) << ',' << std::pow(m->GetBinError(bin),2) << '\n';
        }
        events.flush(); files.flush(); availability.flush(); cuts.flush(); reference.flush();
        require(events && files && availability && cuts && reference, "output write failed");
        events.close(); files.close(); availability.close(); cuts.close(); reference.close();
        std::ofstream complete(directory + "/COMPLETE");
        complete << "Export completed. Analyze with analyze_fei_diagnostics.py.\n";
        std::cout << "\nDiagnostics exported to " << directory << std::endl;
    }
};
inline Exporter*& instance() { static Exporter* p=nullptr; return p; }
inline Exporter& output() { require(instance()!=nullptr,"exporter not initialized"); return *instance(); }
inline void start(const std::string& dir, const std::string& hash) {
    delete instance(); instance()=nullptr;
    instance()=new Exporter(dir,hash);
}
inline void checkFiles(const std::vector<std::string>& names, const char* dirname) {
    require(!names.empty(), "no ROOT input files in " + std::string(dirname));
}
inline std::string inputDir(const char* env, const std::string& fallback, const std::string& suffix) {
    const char* value=gSystem->Getenv(env);
    return value && *value ? std::string(value) + "/" + suffix : fallback;
}
} // namespace feidiag
#endif
