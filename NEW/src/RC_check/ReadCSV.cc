#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "TH1.h"
#include "TCanvas.h"

int main() {
    // Open the CSV file
    std::ifstream file("your_file.csv");
    std::string line;

    // Create a vector of vectors to store each column
    std::vector<std::vector<double>> columns;

    // Initialize columns with an appropriate number of empty vectors
    // You need to determine the number of columns first
    std::getline(file, line);
    std::stringstream ss(line);
    std::string value;
    while (std::getline(ss, value, ',')) {
        columns.push_back(std::vector<double>());
    }

    // Move back to the start of the file after reading the first line
    file.clear();
    file.seekg(0);

    // Read each line of the file
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        int colIndex = 0;

        // Split the line into values and add each one to the corresponding column
        while (std::getline(ss, value, ',')) {
            columns[colIndex].push_back(std::stod(value)); // Convert to double and store in the column
            colIndex++;
        }
    }

    // Close the file
    file.close();

    // Print the first column as an example
    for (const auto& val : columns[0]) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    TH1D* third_bin = new TH1D("third_bin", ";number of data in 3rd bin;count", 6, 1.5, 7.5);

    for (int i = 0; i < columns.at(2).size(); i++) {
        third_bin->Fill(columns.at(2).at(i));
    }

    TCanvas* c = new TCanvas("c", "", 800, 800);
    third_bin->Draf("Hist");
    c->SaveAs("third_bin.png");

    delete c;

    return 0;
}