#include <iostream>   
#include <algorithm>  
#include <string>
#include <vector>
#include <numeric>

// arXiv:1409.4557v2
# define TB0 1.5195 // (Table. 1)
# define TBp 1.6384 // (Table. 1)
# define BR_Kplus_nunubar 0.00000398 // (eq. 10)
# define BR_K0star_nunubar 0.00000919 // (eq. 11)
# define BR_K0_nunubar (BR_Kplus_nunubar*TB0/TBp) // under (eq. 15)
# define BR_Kplusstar_nunubar (BR_K0star_nunubar*TBp/TB0) // under (eq. 15)
# define BR_Xs_nunubar 0.000029 // (eq. 23)
# define BR_Xsu_nonresonant_nunubar (BR_Xs_nunubar - BR_Kplus_nunubar - BR_Kplusstar_nunubar)
# define BR_Xsd_nonresonant_nunubar (BR_Xs_nunubar - BR_K0_nunubar - BR_K0star_nunubar)

// https://confluence.desy.de/pages/viewpage.action?pageId=107054222
# define N_BpBp_1invab 565400000.0
# define N_B0B0_1invab 534600000.0

# define N_Kplus_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplus_nunubar)
# define N_Kplusstar_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplusstar_nunubar)
# define N_Xsu_nonresonant_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Xsu_nonresonant_nunubar)
# define N_K0_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0_nunubar)
# define N_K0star_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0star_nunubar)
# define N_Xsd_nunubar_1invab (2.0 * N_B0B0_1invab * BR_Xsd_nonresonant_nunubar)

// my MC sample number
# define N_Kplus_nunubar 10000000.0
# define N_K0_nunubar 10000000.0
# define N_Kplusstar_nunubar 10000000.0
# define N_K0star_nunubar 10000000.0
# define N_Xsu_nonresonant_nunubar 50000000.0
# define N_Xsd_nonresonant_nunubar 50000000.0

// scale factor for each MC sample
# define Scale_Kplus (N_Kplus_nunubar_1invab/N_Kplus_nunubar)
# define Scale_Kplusstar (N_Kplusstar_nunubar_1invab/N_Kplusstar_nunubar)
# define Scale_Xsu_nonresonant (N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_nunubar)
# define Scale_K0 (N_K0_nunubar_1invab/N_K0_nunubar)
# define Scale_K0star (N_K0star_nunubar_1invab/N_K0star_nunubar)
# define Scale_Xsd_nonresonant (N_Xsd_nunubar_1invab/N_Xsd_nonresonant_nunubar)

bool sorting(pair<double, double> a, pair<double, double> b) {
	
	return a.first < b.first;

}

void load_files(const char* dirname, std::vector<std::string>* names) {
	TSystemDirectory dir(dirname, dirname);
	TList* files = dir.GetListOfFiles();
	if (files) {
		TSystemFile* file;
		TString fname;
		TIter next(files);
		while ((file = (TSystemFile*)next())) {
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(".root")) {
				names->push_back(fname.Data());
			}
		}
	}
}

class CDF {
private:
	std::vector<double> m_values;
	std::vector<double> m_weights;
	double m_weights_sum;

	std::vector<double> m_signal_output;
	std::vector<double> m_signal_weight;
	std::vector<double> m_background_output;
	std::vector<double> m_background_weight;

	template <typename A, typename B>
	void zip(const std::vector<A>& a, const std::vector<B>& b, std::vector<std::pair<A, B>>& zipped)
	{
		for (size_t i = 0; i < a.size(); ++i)
		{
			zipped.push_back(std::make_pair(a[i], b[i]));
		}
	}

	template <typename A, typename B>
	void unzip(const std::vector<std::pair<A, B>>& zipped, std::vector<A>& a, std::vector<B>& b)
	{
		for (size_t i = 0; i < a.size(); i++)
		{
			a[i] = zipped[i].first;
			b[i] = zipped[i].second;
		}
	}


public:
	void init(std::vector<double> values, std::vector<double> weights) {
		// Zip the vectors together
		std::vector<std::pair<double, double>> zipped;
		zip(values, weights, zipped);

		sort(zipped.begin(), zipped.end(), sorting);

		// Write the sorted pairs back to the vectors
		unzip(zipped, m_values, m_weights);
		m_weights_sum = std::accumulate(m_weights.begin(), m_weights.end(), 0.0);
	}
	void initbypath(const char* dirname, double weight) {
		std::vector<string> names;

		load_files(dirname, &names);

		for (unsigned int i = 0; i < names.size(); i++) {
			float MVA_Continuum = 0;

			TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

			TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
			tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_Continuum);

			for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
				tree_upsilon->GetEntry(j);

				m_values.push_back(MVA_Continuum);
				m_weights.push_back(weight);

			}
			input_file->Close();
		}

		// Zip the vectors together
		std::vector<std::pair<double, double>> zipped;
		zip(m_values, m_weights, zipped);

		sort(zipped.begin(), zipped.end(), sorting);

		// Write the sorted pairs back to the vectors
		unzip(zipped, m_values, m_weights);
		m_weights_sum = std::accumulate(m_weights.begin(), m_weights.end(), 0.0);
	}
	double GetCDFvalue(double value) {
		auto larger_or_equal = std::lower_bound(m_values.begin(), m_values.end(), value);
		auto larger = std::upper_bound(m_values.begin(), m_values.end(), value);

		int until_index = -1;
		if (larger_or_equal - m_values.begin() == larger - m_values.begin()) { // input value is not same with one of value which is saved
			until_index = larger_or_equal - m_values.begin();
		}
		else { // input value is same with one of value which is saved
			until_index = larger_or_equal - m_values.begin();
			double weights_partial_sum = std::accumulate(m_weights.begin(), until_index + m_weights.begin(), 0.0);
			return weights_partial_sum / m_weights_sum;
		}

		if (larger == m_values.end()) { // it is last
			return 1.0;
		}

		double weights_partial_sum_back = std::accumulate(m_weights.begin(), until_index + m_weights.begin(), 0.0);
		double weights_partial_sum_front = std::accumulate(m_weights.begin(), until_index + m_weights.begin() + 1, 0.0);

		int index_back = until_index;
		int index_front = until_index + 1;

		double weights_partial_sum = weights_partial_sum_back + (value - m_values.at(index_back)) * (weights_partial_sum_front - weights_partial_sum_back) / (m_values.at(index_front) - m_values.at(index_back));

		return weights_partial_sum / m_weights_sum;

	}
	void print() {
		for (unsigned int i = 0; i < m_values.size(); i++) {
			printf("[%lf, %lf] ", m_values.at(i), m_weights.at(i));
		}
		printf("\n");
	}
	void Evaluate(const char* dirname, double weight, bool IsSignal) {
		std::vector<string> names;

		load_files(dirname, &names);

		for (unsigned int i = 0; i < names.size(); i++) {
			float MVA_Continuum = 0;

			TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

			TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
			tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_Continuum);

			for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
				tree_upsilon->GetEntry(j);

				if (IsSignal) {
					m_signal_output.push_back(GetCDFvalue(MVA_Continuum));
					m_signal_weight.push_back(weight);
				}
				else {
					m_background_output.push_back(GetCDFvalue(MVA_Continuum));
					m_background_weight.push_back(weight);
				}
			}
			input_file->Close();

		}
	}
	void PrintPNG() {
		TH1D* temp_signal_hist = new TH1D("CDF_signal", ";CDF;arbitrary unit", 100, 0.0, 1.0);
		TH1D* temp_background_hist = new TH1D("CDF_background", ";CDF;arbitrary unit", 100, 0.0, 1.0);

		for (unsigned int i = 0; i < m_signal_output.size(); i++) temp_signal_hist->Fill(m_signal_output.at(i), m_signal_weight.at(i));
		for (unsigned int i = 0; i < m_background_output.size(); i++) temp_background_hist->Fill(m_background_output.at(i), m_background_weight.at(i));

		temp_signal_hist->Scale(1.0 / temp_signal_hist->Integral(), "width");
		temp_background_hist->Scale(1.0 / temp_background_hist->Integral(), "width");

		TCanvas* c_temp = new TCanvas("c", "", 800, 800);

		temp_signal_hist->Draw("hist"); c_temp->SaveAs("signal_cdf.png");
		temp_background_hist->Draw("hist"); c_temp->SaveAs("background_cdf.png");

		delete c_temp;
	}
};

int GetCDF() {

	std::vector<double> values;
	std::vector<double> weights;

	CDF cdf;

	cdf.initbypath("./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu", Scale_Kplus);
	cdf.initbypath("./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu", Scale_Kplusstar);
	cdf.initbypath("./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu", Scale_Xsu_nonresonant);
	cdf.initbypath("./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu", Scale_K0);
	cdf.initbypath("./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu", Scale_K0star);
	cdf.initbypath("./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu", Scale_Xsd_nonresonant);

	cdf.Evaluate("./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu", Scale_Kplus, true);
	cdf.Evaluate("./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu", Scale_Kplusstar, true);
	cdf.Evaluate("./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu", Scale_Xsu_nonresonant, true);
	cdf.Evaluate("./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu", Scale_K0, true);
	cdf.Evaluate("./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu", Scale_K0star, true);
	cdf.Evaluate("./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu", Scale_Xsd_nonresonant, true);

	cdf.Evaluate("./CHG_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge", 1.0, false);
	cdf.Evaluate("./MIX_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge", 1.0, false);
	cdf.Evaluate("./UUBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge", 1.0, false);
	cdf.Evaluate("./DDBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge", 1.0, false);
	cdf.Evaluate("./SSBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge", 1.0, false);
	cdf.Evaluate("./CHARM_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge", 1.0, false);

	cdf.PrintPNG();

	return 0;
}