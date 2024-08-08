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
void get2Rof11B(){
  std::vector<std::string> file_names;
    //std::string folderPath = "/home/niuyujie/GEMINI++4nu/ROOT/Yosoi_11B_1";
    std::string folderPath = "/home/niuyujie/GEMINI++4nu/ROOT/Yosoi_11B_0.5";
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
       if (number < 16 && number > 35){
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
        int PDGid[10] = {0};
        float Px[10] = {0};
        float Py[10] = {0};
        float Pz[10] = {0};
        float E[10] = {0};
        float KE[10] = {0};
        tree->SetBranchAddress("nParticles", &nParticles);
        tree->SetBranchAddress("finalEx", &finalEx);
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
            if(PDGid[0] == 1000010030 && finalEx < 5){
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
      }
    }
}
