// g++ getResult.cpp $(root-config --cflags --libs) -o getResult
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

int main(int argc, char* argv[]){
    std::string checkType = argv[1];
    std::string nucleus = argv[2];
    std::string folderPath = argv[3];
    std::vector<std::string> file_names;
    std::sort(file_names.begin(), file_names.end(), compareNumbersInFileName);
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
        //return EXIT_FAILURE;
    }
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
    if (checkType == "--Panin"){
        std::cout<<"Comparison to Panin et al."<<std::endl;
        if (nucleus != "11B"){
            std::cerr<<"Invalid nucleus"<<std::endl;
            return EXIT_FAILURE;
        }
        int total_n = 0;
        int total_dOrA = 0;
        for (const auto& name : file_names) {
            std::cout << name << std::endl;
            std::string numberStr;
            for (char c : name) {
                if (std::isdigit(c) || c == '.') { numberStr += c; }
            }
        float number;
        std::istringstream(numberStr) >> number;
        std::cout<<number<<std::endl;
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
            // std::cout<<nEv<<std::endl;
            int nParticles;
            int PDGid[10] = {0};
            float Px[10] = {0};
            float Py[10] = {0};
            float Pz[10] = {0};
            float E[10] = {0};
            float KE[10] = {0};
            tree->SetBranchAddress("nParticles", &nParticles);
            tree->SetBranchAddress("PDGid", PDGid);
            tree->SetBranchAddress("Px", Px);
            tree->SetBranchAddress("Py", Py);
            tree->SetBranchAddress("Pz", Pz);
            tree->SetBranchAddress("E", E);
            tree->SetBranchAddress("KE", KE);
            int needN = 0;
            int needA = 0;
            int needD = 0;
            int needT = 0;
            int needP = 0;
            for (int i = 0; i < nEv; ++i){
                tree->GetEntry(i); 
                //cout<<productNames->at(j)<<endl;
                if( PDGid[0] == 2112 && PDGid[1] == 1000050100){
                    needN++;
                }
                if((PDGid[0] == 1000020040 && PDGid[1] == 1000030070 ) ||( PDGid[0] == 1000020030 && PDGid[1] == 1000030080) ){
                    needA++;
                }
                if(PDGid[0] == 1000010020 && PDGid[1] == 1000040090){
                    needD++;
                }
                if(PDGid[0] == 1000010030 && (number - KE[0] * 11./8. - 11.2 < 5)){//similar to Yosoi et al.
                    needT++;
                }
                if(PDGid[0] == 2212 && PDGid[1] == 1000040100){
                    needP++;
                }
            }
                std::cout<<"For this Ex, its n+10B counter is:  "  <<needN <<std::endl;
                std::cout<<"For this Ex, its d+9Be counter is:  "  <<needD <<std::endl;
                std::cout<<"For this Ex, its 4He+7Li counter is:  "<<needA <<std::endl;
                std::cout<<"For this Ex, its t+8Be counter is:  "  <<needT <<std::endl;
                std::cout<<"For this Ex, its p+10Be counter is:  " <<needP <<std::endl;
                if (number >= 16 && number <= 35){
                    total_n += needN;
                    total_dOrA += needD + needA;
                }
            }
        }
        std::cout<<"Total n+10B counter is:  "  <<total_n <<std::endl;
        std::cout<<"Total d+9Be+4He+7Li counter is:  "  <<total_dOrA <<std::endl;
        std::cout<< "Ratio of n+10B is: " << total_n*1.0 / (total_dOrA +  total_n) << std::endl;
    } else if (checkType == "--Yosoi"){
        std::cout<<"Comparison to Yosoi et al."<<std::endl;
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
        if (nucleus == "15N"){
            std::cout<<"15N result: "<<std::endl;
            for (const auto& name : file_names) {
              std::cout << name << std::endl;
              int temp_n_neutron = 0;
              int temp_n_proton = 0;
              int temp_n_deutron = 0;
              int temp_n_triton = 0;
              int temp_n_alpha = 0;
              int temp_n2_neutron = 0;
              int temp_n2_proton = 0;
              int temp_n2_deutron = 0;
              int temp_n2_triton = 0;
              int temp_n2_alpha = 0;
              std::string numberStr;
              for (char c : name) {
                  if (std::isdigit(c) || c == '.') { numberStr += c;}
              }
              float number;
              std::istringstream(numberStr) >> number;
              std::cout<<number<<std::endl;
            //   if (number < 20 || number > 40){ continue; }
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
                int nParticles;
                int treePDGid[10] = {0};
                float treePx[10] = {0};
                float treePy[10] = {0};
                float treePz[10] = {0};
                float treeE[10] = {0};
                float treeKE[10] = {0};
                tree->SetBranchAddress("nParticles", &nParticles);
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
                    if (treePDGid[j] == 2212 && treeKE[j] > 3.1){
                      n_proton++;
                      temp_n_proton++;
                      double R_SE = number - treeKE[j] * 11./10.;
                      if( R_SE - 11.2  < 8.2) {
                        n2_proton++;
                        temp_n2_proton++;
                        }
                    }else if (treePDGid[j] == 1000010020 && treeKE[j] > 4.){
                      n_deutron++;
                      temp_n_deutron++;
                      double R_SE = number - treeKE[j] * 11./9.;
                      if( R_SE - 15.8  < 8) {
                        n2_deutron++;
                        temp_n2_deutron++;
                        }
                    }else if (treePDGid[j] == 1000010030 && treeKE[j] > 4.6){
                      n_triton++;
                      temp_n_triton++;
                      double R_SE = number - treeKE[j] * 11./8.;
                      if( R_SE - 11.2  < 8) {
                        n2_triton++;
                        temp_n2_triton++;
                        }
                    }else if ((treePDGid[j] == 1000020040 && treeKE[j] > 4.5) || (treePDGid[j] == 1000020030 && treeKE[j] > 4.5)){
                      n_alpha++;
                      temp_n_alpha++;
                      double R_SE = number - treeKE[j] * 11./7.;
                      if( R_SE - 8.7  < 8.7) {
                        n2_alpha++;
                        temp_n2_alpha++;
                        }
                    }else if (treePDGid[j] == 2112 && treeKE[j] > 3.2){
                      n_neutron++;
                      temp_n_neutron++;
                      double R_SE = number - treeKE[j] * 11./10.;
                      if( R_SE - 11.5  < 8.2) {
                        n2_neutron++;
                        temp_n2_neutron++;
                        }
                    }
                  }
                }
              }
              std::cout << "The Ex is: " << number << std::endl;
              std::cout << "For this Ex, its n counter is: " << temp_n_neutron << std::endl;
              std::cout << "For this Ex, its p counter is: " << temp_n_proton << std::endl;
              std::cout << "For this Ex, its d counter is: " << temp_n_deutron << std::endl;
              std::cout << "For this Ex, its t counter is: " << temp_n_triton << std::endl;
              std::cout << "For this Ex, its 4He counter is: " << temp_n_alpha << std::endl;
              if (number < 20 || number > 40){ 
                n_neutron -= temp_n_neutron;
                n_proton -= temp_n_proton;
                n_deutron -= temp_n_deutron;
                n_triton -= temp_n_triton;
                n_alpha -= temp_n_alpha;
                n2_neutron -= temp_n2_neutron;
                n2_proton -= temp_n2_proton;
                n2_deutron -= temp_n2_deutron;
                n2_triton -= temp_n2_triton;
                n2_alpha -= temp_n2_alpha;
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
        } else if (nucleus == "11B"){
            std::cout<<"11B result: "<<std::endl;
            for (const auto& name : file_names) {
              std::cout << name << std::endl;
              int temp_n_neutron = 0;
              int temp_n_proton = 0;
              int temp_n_deutron = 0;
              int temp_n_triton = 0;
              int temp_n_alpha = 0;
                int temp_n2_neutron = 0;
                int temp_n2_proton = 0;
                int temp_n2_deutron = 0;
                int temp_n2_triton = 0;
                int temp_n2_alpha = 0;
              std::string numberStr;
              for (char c : name) {
                  if (std::isdigit(c) || c == '.') { numberStr += c; }
              }
              float number;
              std::istringstream(numberStr) >> number;
              // std::cout<<number<<std::endl;
            //   if (number < 16 || number > 35){ continue; }
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
                int treePDGid[10] = {0};
                float treePx[10] = {0};
                float treePy[10] = {0};
                float treePz[10] = {0};
                float treeE[10] = {0};
                float treeKE[10] = {0};
                tree->SetBranchAddress("nParticles", &nParticles);
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
                    if (treePDGid[j] == 2212 && treeKE[j] > 3.1){
                      n_proton++;
                        temp_n_proton++;
                      double R_SE = number - treeKE[j] * 11./10.;
                      if( R_SE - 11.2  < 6.8) {
                        n2_proton++;
                        temp_n2_proton++;
                        }
                    }else if (treePDGid[j] == 1000010020 && treeKE[j] > 4){
                      n_deutron++;
                        temp_n_deutron++;
                      double R_SE = number - treeKE[j] * 11./9.;
                      if( R_SE - 15.8  < 5) {
                        n2_deutron++;
                        temp_n2_deutron++;
                        }
                    }else if (treePDGid[j] == 1000010030 && treeKE[j] > 4.6){
                      n_triton++;
                        temp_n_triton++;
                      double R_SE = number - treeKE[j] * 11./8.;
                      if( R_SE - 11.2  < 5) {
                        n2_triton++;
                        temp_n2_triton++;
                        }
                    }else if ((treePDGid[j] == 1000020040 && treeKE[j] > 4.5) || (treePDGid[j] == 1000020030 && treeKE[j] > 4.5)){
                      n_alpha++;
                        temp_n_alpha++;
                      double R_SE = number - treeKE[j] * 11./7.;
                      if( R_SE - 8.7  < 5) {
                        n2_alpha++;
                        temp_n2_alpha++;
                        }
                    }else if (treePDGid[j] == 2112 && treeKE[j] > 3.1){
                      n_neutron++;
                        temp_n_neutron++;
                      double R_SE = number - treeKE[j] * 11./10.;
                      if( R_SE - 11.5  < 5) {
                        n2_neutron++;
                        temp_n2_neutron++;
                        }
                    }
                  }
                }
              }
                std::cout << "The Ex is: " << number << std::endl;
                std::cout << "For this Ex, its n counter is: " << temp_n_neutron << std::endl;
                std::cout << "For this Ex, its p counter is: " << temp_n_proton << std::endl;
                std::cout << "For this Ex, its d counter is: " << temp_n_deutron << std::endl;
                std::cout << "For this Ex, its t counter is: " << temp_n_triton << std::endl;
                std::cout << "For this Ex, its 4He counter is: " << temp_n_alpha << std::endl;
              if (number < 16 || number > 35){ 
                n_neutron -= temp_n_neutron;
                n_proton -= temp_n_proton;
                n_deutron -= temp_n_deutron;
                n_triton -= temp_n_triton;
                n_alpha -= temp_n_alpha;
                n2_neutron -= temp_n2_neutron;
                n2_proton -= temp_n2_proton;
                n2_deutron -= temp_n2_deutron;
                n2_triton -= temp_n2_triton;
                n2_alpha -= temp_n2_alpha;
              }

            }
            double total = 277925.0;
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
        } else {
            std::cerr<<"Invalid nucleus"<<std::endl;
            return EXIT_FAILURE;
        }
    } else {
        std::cerr<<"Invalid type"<<std::endl;
        return EXIT_FAILURE;
    }
}