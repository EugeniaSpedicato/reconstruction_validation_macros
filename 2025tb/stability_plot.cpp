#include <iostream>
#include <fstream>
#include <string>
#include <regex>

#include "TApplication.h"
#include "TGraph.h"
#include "TCanvas.h"

int main(int argc, char* argv[]) {



    if (argc < 2) {
        std::cerr << "Uso: " << argv[0] << " <file.txt>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];

    // Controlla se il nome del file contiene la stringa desiderata
    if (filename.find("single_muon_interaction_0") != std::string::npos) {
        std::cout << "Il file contiene 'single_muon_interaction_0' nel nome." << std::endl;

    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Errore nell'apertura del file: " << filename << std::endl;
        return 1;
    }

    std::string line;
    std::regex pattern(R"(Ratio good vrtx tr0:\s+([0-9]*\.?[0-9]+))");
    std::smatch match;

    while (std::getline(infile, line)) {
        if (std::regex_search(line, match, pattern)) {
            if (match.size() > 1) {
                std::cout << "Valore estratto: " << match[1] << std::endl;
                return 0;
            }
        }
    }
  }//if sing_int_1
  else  if (filename.find("single_muon_interaction_1") != std::string::npos) {
        std::cout << "Il file contiene 'single_muon_interaction_1' nel nome." << std::endl;

    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Errore nell'apertura del file: " << filename << std::endl;
        return 1;
    }

    std::string line;
    std::regex pattern(R"(Ratio good vrtx tr1:\s+([0-9]*\.?[0-9]+))");
    std::smatch match;

    while (std::getline(infile, line)) {
        if (std::regex_search(line, match, pattern)) {
            if (match.size() > 1) {
                std::cout << "Valore estratto: " << match[1] << std::endl;
                return 0;
            }
        }
    }

  }//if sing_int_1

    return 0;
}

