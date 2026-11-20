"""Build the instrumented macro from the user's preserved source; no weight rewrite."""
from pathlib import Path
import hashlib
import re

HERE = Path(__file__).resolve().parent


def replace_exact(text, old, new, count=1):
    if text.count(old) != count:
        raise ValueError(f'Expected {count} instances of {old[:90]!r}; got {text.count(old)}')
    return text.replace(old, new)


def generate(source):
    digest = hashlib.sha256(source).hexdigest()
    s = source.decode('utf-8-sig').replace('\r\n', '\n')
    s = '#include "FEIDiagnosticExport.h"\n' + s
    s = s.replace('LetsFillKonlyMC', 'FEIDiagFillMC').replace('LetsFillKonlydata', 'FEIDiagFillData')
    # Check each original required branch. Historical pi0 arrays unused by MC15rd
    # remain optional; their absence cannot affect the MC15rd weight expression.
    s = re.sub(r'(tree_\w+)->SetBranchAddress\((.*)\);', r'feidiag::bind(\1, \2);', s)
    for fragment in ['N_pi0_syst; i_pi0++', 'N_pi0_syst_MC15ri; i_pi0++']:
        lines = s.splitlines()
        for i,line in enumerate(lines):
            if fragment in line and 'feidiag::bind' in line:
                lines[i] = line.replace(');', ', std::string(MCTYPE) == "MC15ri");')
        s = '\n'.join(lines) + '\n'
    s = s.replace('unsigned int j = 0; j < tree_upsilon->GetEntries()', 'Long64_t j = 0; j < tree_upsilon->GetEntries()')
    s = re.sub(r'(tree_\w+)->GetEntry\(j\);', r'feidiag::read(\1, j);', s)
    s = replace_exact(s, '        TTree* tree_Xs;', '        TTree* tree_Xs = nullptr;')
    for load in ['load_files(dirname, &names);', 'load_files(dirname, &names, included_string);']:
        s = replace_exact(s, load, load+'\n    feidiag::checkFiles(names, dirname);\n    std::sort(names.begin(), names.end());')
    anchor = '        feidiag::bind(tree_upsilon, "extraInfo__bodecayModeID__bc", &Upsilon_ID);'
    mc_register = '''        const double diag_norm = CAL * ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i)) * additional_weight;
        feidiag::output().file(input_file, std::string(dirname)+"/"+names.at(i), SampleName,
                              tree_upsilon, tree_Bsig, tree_Btag, tree_Xs, diag_norm);
'''
    data_register = '''        const double diag_norm = 1.0;
        feidiag::output().file(input_file, std::string(dirname)+"/"+names.at(i), "DATA",
                              tree_upsilon, tree_Bsig, tree_Btag, nullptr, diag_norm);
'''
    if s.count(anchor) != 2: raise ValueError('Missing binding anchors')
    parts=s.split(anchor)
    s=parts[0]+mc_register+anchor+parts[1]+data_register+anchor+parts[2]
    read_mc='            if (SampleName == "SIGNAL") feidiag::read(tree_Xs, j);'
    s=replace_exact(s,read_mc,read_mc+'''
            feidiag::output().input(SampleName, Upsilon_ID, Btag_ID, Bsig_ID, Bsig_M,
                                   diag_norm, IncludeKstar, option);''')
    data_start=s.index('void FEIDiagFillData')
    a,b=s[:data_start],s[data_start:]
    b=replace_exact(b,'            feidiag::read(tree_Btag, j);','''            feidiag::read(tree_Btag, j);
            feidiag::output().input("DATA", Upsilon_ID, Btag_ID, Bsig_ID, Bsig_M,
                                   diag_norm, IncludeKstar, option);''')
    s=a+b
    weight_line=next(line for line in s.splitlines() if 'double weights = ' in line)
    s=replace_exact(s,weight_line,weight_line+'''
            const double diag_official = (SampleName == "SIGNAL" || SampleName == "CHG" || SampleName == "MIX")
                ? corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE) : 1.0;
            feidiag::output().selected(j, SampleName, Upsilon_ID, Btag_ID, Bsig_ID, Bsig_M, weights,
                {diag_norm, FEI_calibration_factor, Correction_pi0, Correction_KID, Correction_PID,
                 Correction_fake, Correction_KS0, Correction_Knn, Correction_Xnn, Correction_multiplicity,
                 Correction_KpKLKL, Correction_KSKLKL, Correction_KstarKLKL, Correction_XKLKL, Correction_BtoDtoXKL}, diag_official);''')
    data_start=s.index('void FEIDiagFillData')
    a,b=s[:data_start],s[data_start:]
    anchor='            if ((Upsilon_ID > -0.5) && (Upsilon_ID < 0.5)) {'
    b=replace_exact(b,anchor,'''            feidiag::output().selected(j, "DATA", Upsilon_ID, Btag_ID, Bsig_ID, Bsig_M,
                                       1.0, std::vector<double>(15,1.0), 1.0);
'''+anchor)
    s=a+b
    histograms=re.findall(r'TH1D\* (\w+) = new TH1D',s)
    init='\n    feidiag::start(output_dir, "'+digest+'");\n'
    init+='    for (TH1D* h : {'+', '.join(histograms)+'}) { h->Reset(); h->Sumw2(); }\n'
    s=replace_exact(s,'void THStack_Jpsi_FEI_efficiency() {',
                    'void THStack_Jpsi_FEI_diagnostics(const char* output_dir="fei_diagnostics_output") {'+init)
    # Optional relocations preserve the source's paths by default.
    def relocate(m):
        variable,path=m.group(1),m.group(2)
        env='FEI_DIAG_DATA_BASE' if 'data_before' in variable else 'FEI_DIAG_MC_BASE'
        suffix=path.split('_Jpsi/',1)[1]
        return f'const std::string {variable}_value = feidiag::inputDir("{env}", "{path}", "{suffix}");\n    const char* {variable} = {variable}_value.c_str();'
    s=re.sub(r'const char\* (Jpsi_\w+_before_dirname) = "([^"]+)";',relocate,s)
    s=replace_exact(s,'FILE* fp = fopen("NEWFEIcal.txt", "w");',
                    'FILE* fp = fopen((std::string(output_dir)+"/NEWFEIcal_reproduced.txt").c_str(), "w");\n    feidiag::require(fp != nullptr, "cannot write reproduced factors");')
    s=replace_exact(s,'    fclose(fp);','''    fclose(fp);
    feidiag::output().finish(Nevt_MC_Bplus_My, Nevt_data_Bplus_My, Nevt_MC_Bzero_My, Nevt_data_Bzero_My);''')
    # Shared calibration headers and their working-directory data files are those
    # of the user's original ROOT environment, not copies from another analysis.
    return s


if __name__ == '__main__':
    source=(HERE/'source/THStack_Jpsi_FEI_efficiency.c').read_bytes()
    target=HERE/'THStack_Jpsi_FEI_diagnostics.C'
    target.write_text(generate(source),encoding='utf-8',newline='\n')
    print(target)
