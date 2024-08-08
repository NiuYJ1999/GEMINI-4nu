#include <iostream>
#include <sstream>
#include "TCanvas.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TTreeReaderValue.h"
#include "TTreeReader.h"
#include <vector>
#include <algorithm>
#include <string>
#include <dirent.h>
bool compareNumbersInFileName(const std::string& a, const std::string& b) {
    auto extractNumber = [](const std::string& filename) {
        size_t start_index = filename.find("Ex_") + 3;
        size_t end_index = filename.find(".root");
        if (start_index != std::string::npos && end_index != std::string::npos) {
            return std::stod(filename.substr(start_index, end_index - start_index));
        }
        return std::numeric_limits<double>::max(); 
    };
    return extractNumber(a) < extractNumber(b);
}
void getRatioYosoi_15N(){
  int n_neutron = 0;
  int n_proton = 0;
  int n_deutron = 0;
  int n_triton = 0;
  int n_alpha = 0;
  int n2_neutron = 0;
  int n2_proton = 0;
  int n2_deutron = 0;
  int n2_triton = 0;
  int n2_alpha = 0;
  std::vector<std::string> file_names;
  std::string folderPath = "/home/niuyujie/GEMINI++4nu/ROOT/Yosoi_15N_1";
  // std::string folderPath = "/home/niuyujie/GEMINI++4nu/ROOT/Yosoi_15N_0.5";
    DIR* dir = opendir(folderPath.c_str());
    struct dirent *ent;
    if (dir != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            std::string filename = ent->d_name;
            if (filename.find(".root") != std::string::npos) {
                file_names.push_back(filename);
            }
        }
        closedir(dir);
    } else {
        std::cerr << "Error opening directory" << std::endl;
        // return EXIT_FAILURE;
    }

    std::sort(file_names.begin(), file_names.end(), compareNumbersInFileName);

    for (const auto& name : file_names) {
      std::cout << name << std::endl;
      std::string numberStr;
      for (char c : name) {
          if (std::isdigit(c) || c == '.') {
              numberStr += c;
          }
      }
      float number;
      std::istringstream(numberStr) >> number;
      std::cout<<number<<std::endl;
      if (number < 20 || number > 40){
        continue;
      }
      TString filename = folderPath + "/" + name;
      TFile* file = new TFile(filename, "READ");
      if (file->IsZombie()){
        std::cout << "Error: File is a Zombie!" << std::endl;
      }else{
        TTree* tree = (TTree*)file->Get("deexTree");
        if (tree == nullptr){
          std::cout << "Error: Tree not found!" << std::endl;
          continue;
        }
        int nEv = tree->GetEntries();
        std::cout<<nEv<<std::endl;
        int nParticles;
        float finalEx;
        int treePDGid[10] = {0};
        float treePx[10] = {0};
        float treePy[10] = {0};
        float treePz[10] = {0};
        float treeE[10] = {0};
        float treeKE[10] = {0};
        tree->SetBranchAddress("nParticles", &nParticles);
        tree->SetBranchAddress("finalEx", &finalEx);
        tree->SetBranchAddress("PDGid", treePDGid);
        tree->SetBranchAddress("Px", treePx);
        tree->SetBranchAddress("Py", treePy);
        tree->SetBranchAddress("Pz", treePz);
        tree->SetBranchAddress("E", treeE);
        tree->SetBranchAddress("KE", treeKE);
        for (int i = 0; i < nEv; ++i){
          tree->GetEntry(i); 
          // std::cout<<nParticles<<std::endl;
          for (int j = 0; j < nParticles; j++){
            // std::cout<<treePDGid[j]<<std::endl;
            if (treePDGid[j] == 2212 && treeKE[j] > 3.2){
              n_proton++;
              double R_SE = number - treeKE[j] * 11./10.;
              if( R_SE - 11.2  < 8.2) {n2_proton++;}
            }else if (treePDGid[j] == 1000010020 && treeKE[j] > 4.2){
              n_deutron++;
              double R_SE = number - treeKE[j] * 11./9.;
              if( R_SE - 15.8  < 8) {n2_deutron++;}
            }else if (treePDGid[j] == 1000010030 && treeKE[j] > 4.8){
              n_triton++;
              double R_SE = number - treeKE[j] * 11./8.;
              if( R_SE - 11.2  < 8) {n2_triton++;}
            }else if ((treePDGid[j] == 1000020040 && treeKE[j] > 5.5) || (treePDGid[j] == 1000020030 && treeKE[j] > 5.5)){
              n_alpha++;
              double R_SE = number - treeKE[j] * 11./7.;
              if( R_SE - 8.7  < 8.7) {n2_alpha++;}
            }else if (treePDGid[j] == 2112 && treeKE[j] > 3.2){
              n_neutron++;
              double R_SE = number - treeKE[j] * 11./10.;
              if( R_SE - 11.5  < 8.2) {n2_neutron++;}
            }
          }
        }
      }
    }
    double total = 126965.0;
    std::cout<< n_neutron / total<< " "<<std::endl;
    std::cout<< n_proton / total<< " "<<std::endl;
    std::cout<< n_deutron / total<< " "<<std::endl;
    std::cout<< n_triton / total<< " "<<std::endl;
    std::cout<< n_alpha / total<< " "<<std::endl;

    std::cout<< n2_neutron / total<< " "<<std::endl;
    std::cout<< n2_proton / total<< " "<<std::endl;
    std::cout<< n2_deutron / total<< " "<<std::endl;
    std::cout<< n2_triton / total<< " "<<std::endl;
    std::cout<< n2_alpha / total<< " "<<std::endl;
}
