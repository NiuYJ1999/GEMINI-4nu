#include "CNucleus.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <array>
#include <random>
#include <unordered_map>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH1F.h"

// this is an example of using GEMINI CNucleus class to give the
// statistical decay of a compound nucleus

// level information
struct levelMap{
    std::string nucleusName;
    std::vector<float> levelEnergy;
    std::vector<std::string> emitting;
    std::vector<float> branchingRatios;
};

// nucleus information designed for only nucleus A < 15 and Z < 7
struct nucleusTable{
  std::string nucelusName;
  int Z;
  int A;
  double mass;
  double massExcess;
};

// generate random number
double randomDouble() {
    std::random_device rd;  
    std::mt19937 gen(rd());  
    std::uniform_real_distribution<double> dis(0.0, 1.0); 
    return dis(gen);  
}

// get the information from the mass table and level table, load the information
void loadMap(std::vector<levelMap> &levels,  std::vector<nucleusTable> &massInfo, std::unordered_map<std::string, nucleusTable> &Mass){
    std::ifstream levelFile;
    std::ifstream massFile;
    const char* Gpath = std::getenv("GINPUT");
    // std::string Gpath = GEMINIpath;
    if (Gpath == nullptr) {
        std::cerr << "Environment variable GINPUT is not set" << std::endl;
    } else {
        // std::cout << "Custom path set to: " << Gpath << std::endl;
        levelFile.open(std::string(Gpath) + "/Level.csv");
        massFile.open(std::string(Gpath) + "/Nucleus.csv");
    }
    if (! levelFile.is_open()) {
        std::cout << "Failed to open file." << std::endl;
        return;
    }else{
        // std::cout << "File opened successfully." << std::endl;
        levelMap mapL;
        std::string parentName0 = "";
        std::string parentName;
        double levelE;
        double levelE0 = -1.;
        std::string daughterName;
        double branchingRatio;
        // Skip the first line, to read the real information
        std::string line;
        std::getline(levelFile, line);
        while (std::getline(levelFile, line)) {
            std::stringstream ss(line);
            std::getline(ss, parentName, ',');
            ss >> levelE;
            ss.ignore(1);// ignore the comma
            std::getline(ss, daughterName, ',');
            ss >> branchingRatio;
            // important to ensure the correct order of the levels
            // the variable levelE0 is used to store the previous level energy, at first it is
            // not the reason to be the same as levelE, so that the last time level can be stored
            // the parentName changes means the last level of the previous parent nucleus
            if(parentName0 != parentName){
                parentName0 = parentName;
                levels.push_back(mapL);
                mapL.levelEnergy.clear();
                mapL.emitting.clear();
                mapL.branchingRatios.clear();
            }
            mapL.nucleusName = parentName;
            mapL.levelEnergy.push_back(levelE);
            mapL.emitting.push_back(daughterName);
            mapL.branchingRatios.push_back(branchingRatio);
            // std::cout << parentName << " "
            //           << levelE     << " "
            //           << daughterName   << " "
            //           << branchingRatio <<std::endl;
        }
        levels.push_back(mapL); // last mode
        levelFile.close();
    }
    // for (int i = 0; i < levels.size(); i++){
    //     std::cout<<levels[i].nucleusName<<std::endl;
    //     for (int j = 0; j < levels[i].levelEnergy.size(); j++){
    //         std::cout<<levels[i].levelEnergy[j]<<" "<<levels[i].emitting[j]<<" "<<levels[i].branchingRatios[j]<<std::endl;
    //     }
    // }

    nucleusTable nucleusInfo;
    std::string name;
    double massExcess;
    double mass;
    int Z;
    int A;
    if (!massFile.is_open()) {
      // std::cout << "Failed to open file." << std::endl;
      return;
    }else{
      // std::cout << "File opened successfully." << std::endl;
      // Skip the first line, to read the real information
        std::string line;   
        std::getline(massFile, line);
        while (std::getline(massFile, line)) {
            std::stringstream ss(line);
            std::getline(ss, name, ',');
            ss >> Z ;
            ss.ignore();
            ss >> A;
            ss.ignore();
            ss >> massExcess;
            ss.ignore();
            ss >> mass;
            nucleusInfo.nucelusName = name;
            nucleusInfo.Z = Z;
            nucleusInfo.A = A;
            nucleusInfo.massExcess = massExcess;
            nucleusInfo.mass = mass * 931.4936 - Z * 0.511; // convert the mass to MeV
            Mass[name] = nucleusInfo;
            massInfo.push_back(nucleusInfo);
        }
        massFile.close();
    }
    
}

std::string getResidueName(std::string parentName, std::string productName){
    std::vector<levelMap>                           levels;
    std::vector<nucleusTable>                       massInfo;
    std::unordered_map<std::string, nucleusTable>   Mass;
    loadMap(levels, massInfo, Mass);
    std::string residueName;
    int residueZ, residueA;
    residueZ = Mass[parentName].Z - Mass[productName].Z;
    residueA = Mass[parentName].A - Mass[productName].A;
    for (const auto mass : massInfo){
        if (mass.Z == residueZ && mass.A == residueA){
            residueName = mass.nucelusName;
        }
    }
    return residueName;
}
//modify the kinematcis
// parentMass is the mass of the parent, parentEx is the excited E of parent, parentKE is the kinematic of parent
// product is the same
// residueMass is the residue mass plus its excited energy
std::array<std::array<float, 4>, 2> modifyKinematic(float parentMass, float parentEx, float parentKE, float productMass, float residueM, float residueEx){
    std::array<std::array<float, 4>, 2> result;
    double parentE0 = parentMass + parentEx;
    double residueMass = residueM + residueEx;
    if (parentE0 - productMass - residueMass < 0){
        // std::cout<<"The decay is not allowed!"<<std::endl;
        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 4; j++){
                result[i][j] = 0;
            }
        }
        return result;
    }
    // std::cout<<"parentE0: "<<parentE0<<std::endl;
    // std::cout<<"parentMass: "<<parentMass<<std::endl;
    // std::cout<<"parentEx: "<<parentEx<<std::endl;
    double parentP = sqrt(parentKE * (parentKE+ 2 * parentE0));
    double massDelta = parentE0 * parentE0 + productMass * productMass - residueMass * residueMass;
    double P2 = sqrt ( massDelta * massDelta / (4 * parentE0 * parentE0) - productMass * productMass );
    // double massDelta = parentE0 * parentE0 - productMass * productMass - residueMass * residueMass;
    // double P2 = sqrt (massDelta * massDelta - 4 * productMass * productMass * residueMass * residueMass) / (2 * parentE0);
    // std::cout<<"parentE0: "<<parentE0 <<std::endl;
    // std::cout<<"productMass: "<<productMass <<std::endl;
    // std::cout<<"residueMass: "<<residueMass <<std::endl;
    // std::cout<<"massDelta: "<<massDelta<<std::endl;
    // std::cout<<massDelta * massDelta / (4 * parentE0 * parentE0)<<std::endl;
    // std::cout<<productMass * productMass<<std::endl;
    double productE = sqrt(productMass * productMass + P2 * P2);
    double residueE = sqrt(residueMass * residueMass + P2 * P2);
    float emitTheta = acos(2 * randomDouble() - 1);
    float emitPhi = 2 * CNucleus::pi * randomDouble();
    float ExTheta = acos(2 * randomDouble() - 1);
    float ExPhi = 2 * CNucleus::pi * randomDouble();

    TLorentzVector productVector(P2 * sin(emitTheta) * cos(emitPhi), P2 * sin(emitTheta) * sin(emitPhi), P2 * cos(emitTheta), productE);
    TLorentzVector residueVector(-P2 * sin(emitTheta) * cos(emitPhi), -P2 * sin(emitTheta) * sin(emitPhi), -P2 * cos(emitTheta), residueE);
    float parentE = sqrt(parentP * parentP + parentMass * parentMass);
    float beta = parentP / parentE;
    float betaX = beta * sin(ExTheta) * cos(ExPhi);
    float betaY = beta * sin(ExTheta) * sin(ExPhi);
    float betaZ = beta * cos(ExTheta);

    productVector.Boost(betaX, betaY, betaZ);
    residueVector.Boost(betaX, betaY, betaZ);

    result[0][0] = productVector.Px();
    result[0][1] = productVector.Py();
    result[0][2] = productVector.Pz();
    result[0][3] = productVector.E() - productMass;
    result[1][0] = residueVector.Px();
    result[1][1] = residueVector.Py();
    result[1][2] = residueVector.Pz();
    result[1][3] = residueVector.E() - residueMass;

    // std::cout<< "After modified, KE of product: " <<  result[0][3] << ", residue: "<< result[1][3] << std::endl;
    return result;
}

void setSupressionF(std::string fs_n, std::string fs_p, std::string fs_d, std::string fs_tr, std::string fs_a, std::string fs_h){
    const char* Gpath = std::getenv("GINPUT");
    if (Gpath == nullptr) {
        std::cerr << "Environment variable GINPUT is not set" << std::endl;
    } else {
        // std::cout << "Custom path set to: " << Gpath << std::endl;
        std::ifstream file;
        file.open(std::string(Gpath) + "/tbl/evap.inp");
        std::string name;
        std::string Z;
        std::string A;
        std::string spin;
        std::string Ex;
        std::string Ek;
        std::stringstream modifiedContent;
        float suppress;
        if (file.is_open()) {
            std::vector<std::string> lines;
            std::string line;
            while (std::getline(file, line)) {
                lines.push_back(line);
            }
            file.close();
            lines[2] = "0 1 0.5  0.      neutron   "+ fs_n  +"\t\t0.";
            lines[3] = "1 1 0.5  0.      proton    "+ fs_p  +"\t\t0.";
            lines[4] = "2 4 0.0  0.      alpha     "+ fs_a  +"\t\t0.";
            lines[5] = "1 2 1.0  0.      deuteron  "+ fs_d  +"\t\t0.";
            lines[6] = "1 3 0.5  0.      triton    "+ fs_tr +"\t\t0.";
            lines[7] = "2 3 0.5  0.      he3       "+ fs_h  +"\t\t0.";
            lines[8] = "3 6 1.0  0.      li6       "+ fs_a  +"\t\t0.";
            lines[9] = "3 6 3.0  2.186   li6       "+ fs_a  +"\t\t0.71170   1  2  1.  0.  2    3.  .024";
            lines[10]= "3 6 0.0  3.563   li6       "+ fs_a  +"\t\t0.";
            lines[11]= "3 6 2.0  4.31    li6       "+ fs_a  +"\t\t2.835     1  2  1.  0.  2    2.  1.3";
            lines[12]= "3 7 1.5  0.      li7       "+ fs_a  +"\t\t0.";
            lines[13]= "3 7 0.5  0.47761 li7       "+ fs_a  +"\t\t0.";
            lines[14]= "3 7 3.5  4.630   li7       "+ fs_a  +"\t\t2.185     1  3  0.5 0.0 3    3.5  .093";
            std::ofstream outFile(std::string(Gpath) + "/tbl/evap.inp");
            if (outFile.is_open()) {
                for (const std::string& modifiedLine : lines) {
                    outFile << modifiedLine << std::endl;
                }
                outFile.close();
                // std::cout << "Complete!" << std::endl;
            } else {
                std::cout << "Can not open file to modify" << std::endl;
            }
        } else {
            std::cout << "Can not open file to read" << std::endl;
        }
    }
}

int main(int argc, char const *argv[]){
    TH1F* h1 = new TH1F("h1", "h1", 1000, 0, 15000);
    // load the fact information
    std::vector<levelMap>                           levels;
    std::vector<nucleusTable>                       massInfo;
    std::unordered_map<std::string, nucleusTable>   Mass;
    loadMap(levels, massInfo, Mass);
    
    //input the information to specify one nucleus
    std::string modeSelect      = argv[1];
    // specify the compound nucleus
    int Z_input                 = std::stoi(argv[2]);
    int A_input                 = std::stoi(argv[3]);
    float Ex_input              = std::stof(argv[4]);
    float J_input               = std::stof(argv[5]);
    int simulationTimes         = 0;
    TString pathName = " ";
    if (modeSelect == "--event"){
        if (argc != 8){
            std::cout << "Wrong number of arguments!" << std::endl;
        }
        std::cout<<"Event-by-event mode"<<std::endl;
        std::string supNotice = argv[6];
        if (supNotice == "--suppress"){
            std::string suppression = argv[7];
            if (suppression == "1"){
                std::cout<<"Suppression factor has been set to type I."<<std::endl;
                setSupressionF("1.", "1.", "1.", "1.", "1.", "1.");
            }else if (suppression == "0.5"){
                std::cout<<"Suppression factor has been set to type II."<<std::endl;
                setSupressionF("1.", "0.5", "0.5", "0.5", "0.5", "0.5");
            }else if (suppression == "free"){
                std::string fs_n;
                std::string fs_p;
                std::string fs_d;
                std::string fs_tr;
                std::string fs_a;
                std::string fs_h;
                std::cout<<"Please choose your favor number:"<<std::endl;
                std::cout<<"Give suppression factor of neutron:"<<std::endl;
                std::cin>>fs_n;
                std::cout<<"Give suppression factor of proton:"<<std::endl;
                std::cin>>fs_p;
                std::cout<<"Give suppression factor of deuteron:"<<std::endl;
                std::cin>>fs_d;
                std::cout<<"Give suppression factor of triton:"<<std::endl;
                std::cin>>fs_tr;
                std::cout<<"Give suppression factor of alpha:"<<std::endl;
                std::cin>>fs_a;
                std::cout<<"Give suppression factor of 3he:"<<std::endl;
                std::cin>>fs_h;
                setSupressionF(fs_n, fs_p, fs_d, fs_tr, fs_a, fs_h);
            }
        }else {
            std::cout << "Please choose your suppression factor" << std::endl;
        }
        simulationTimes = 1;
    }else if (modeSelect == "--batch"){
        if (argc != 10){
            std::cout << "Wrong number of arguments!" << std::endl;
        }
        std::cout<<"Batch mode"<<std::endl;
        simulationTimes = std::stoi(argv[6]);
        pathName = argv[7];
        std::string supNotice = argv[8];
        if (supNotice == "--suppress"){
            std::string suppression = argv[9];
            if (suppression == "1"){
                std::cout<<"Suppression factor has been set to type I."<<std::endl;
                setSupressionF("1.", "1.", "1.", "1.", "1.", "1.");
            }else if (suppression == "0.5"){
                std::cout<<"Suppression factor has been set to type II."<<std::endl;
                setSupressionF("1.", "0.5", "0.5", "0.5", "0.5", "0.5");
            }else if (suppression == "default"){
                std::cout<<"Suppression factor has been set to GEMINI++"<<std::endl;
                setSupressionF("1.", "1.", "0.5", "0.5", "1.", "0.5");
            }else if (suppression == "free"){
                std::string fs_n;
                std::string fs_p;
                std::string fs_d;
                std::string fs_tr;
                std::string fs_a;
                std::string fs_h;
                std::cout<<"Please choose your favor number:"<<std::endl;
                std::cout<<"Give suppression factor of neutron:"<<std::endl;
                std::cin>>fs_n;
                std::cout<<"Give suppression factor of proton:"<<std::endl;
                std::cin>>fs_p;
                std::cout<<"Give suppression factor of deuteron:"<<std::endl;
                std::cin>>fs_d;
                std::cout<<"Give suppression factor of triton:"<<std::endl;
                std::cin>>fs_tr;
                std::cout<<"Give suppression factor of alpha:"<<std::endl;
                std::cin>>fs_a;
                std::cout<<"Give suppression factor of 3he:"<<std::endl;
                std::cin>>fs_h;
                setSupressionF(fs_n, fs_p, fs_d, fs_tr, fs_a, fs_h);
            }
        }else {
            std::cout << "Please choose your suppression factor" << std::endl;
        }
    }else if (modeSelect == "--h"){
        std::cout<<"Two modes provided:"<<std::endl;
        std::cout<<"--event: event-by-event mode"<<std::endl;
        std::cout<<"--batch: batch mode"<<std::endl;
        std::cout<<"Four basic parameter for both modes:"<<std::endl;
        std::cout<<"Z, A of specific nucleus, now available and valid for 11B and 15N; Excited erngy of nucleus and its spin"<<std::endl;
        std::cout<< "Other parameter is for suppression factor, you can choose 1, 0.5 or free to set your own suppression factor"<<std::endl;
    }else{
        std::cout<<"Please specify the mode of the simulation: --event or --batch, or use --h to know the parameter setting"<<std::endl;
        return 0;
    }

    int iZCN = Z_input;            // proton number of compound nucleus
    int iACN = A_input;           // mass number of compound nucleus
    CNucleus CN(iZCN, iACN); // constructor
    float fEx = Ex_input;         // excitation energy of compound nucleus
    float fJ = J_input;          // spin of compound nucleus

    CN.setCompoundNucleus(fEx, fJ); // specify the excitation energy and spin

    // CN.setVelocityCartesian(0.32206, 0, 0); // set initial CN velocity to non-zero
    CN.setVelocityCartesian(0); // set initial CN velocity to zero
    CAngle spin(CNucleus::pi / 2, (float)0.);
    CN.setSpinAxis(spin); // set the direction of the CN spin vector

    //batch mode need
    TString     geminiName;
    TFile       *rootFile;
    TTree       *tree;
        
    int treePDGid[10]       = {0};
    float treePx[10]        = {0};
    float treePy[10]        = {0};
    float treePz[10]        = {0};
    float treeE[10]         = {0};
    float treeKE[10]        = {0};
    int                     nParticles;

    if (modeSelect == "--batch") {
        geminiName = pathName + Form("/Ex_%f.root", fEx);
        rootFile = new TFile(geminiName, "RECREATE");
        std::cout<<geminiName<<" had been created! "<<std::endl;
        tree = new TTree("deexTree", "Data from GEMINI++4nu simulation");
        tree->Branch("nParticles",  &nParticles, "nParticles/I");
        tree->Branch("PDGid",       treePDGid,   "PDGid[nParticles]/I");
        tree->Branch("Px",          treePx,      "Px[nParticles]/F");
        tree->Branch("Py",          treePy,      "Py[nParticles]/F");
        tree->Branch("Pz",          treePz,      "Pz[nParticles]/F");
        tree->Branch("E",           treeE,       "E[nParticles]/F");
        tree->Branch("KE",          treeKE,      "KE[nParticles]/F");
    }
    for (int i = 0; i < simulationTimes; ++i){
        // std::cout << "Event: " << i << " of " << simulationTimes << std::endl;
        // CN.setWeightIMF(); // turn on enhanced IMF emission
        CN.decay();        // decay the compound nucleus
        
        // give the parent and product info. to use
        std::vector<std::string> productNamesL;
        std::vector<std::string> productNamesH;
        std::vector<std::string> parentNames;
        std::vector<double> productKEsL;
        std::vector<double> productKEsH;
        std::vector<double> productExsL;
        std::vector<double> productExsH;
        std::vector<double> parentKEs;
        std::vector<double> parentExs;
        std::vector<double> parentPx;
        std::vector<double> parentPy;
        std::vector<double> parentPz;
        std::vector<double> productPLx;
        std::vector<double> productPLy;
        std::vector<double> productPLz;
        std::vector<double> productPHx;
        std::vector<double> productPHy;
        std::vector<double> productPHz;

        // std::vector<std::string>            productName;
        std::vector<int>                    PDGid;
        std::vector<float>                  productPx;
        std::vector<float>                  productPy;
        std::vector<float>                  productPz;
        std::vector<float>                  productKE;
        std::vector<float>                  productE;

        std::vector<double> GammaRayEnergy;
        std::vector<double> GammaRayPx;
        std::vector<double> GammaRayPy;
        std::vector<double> GammaRayPz;

        if (CN.abortEvent){
            CN.reset();
            // --i;
            // continue;
        }

        CNucleus *products = CN.getProducts(0); // set pointer to firt
                                                // stable product
        int j = 0;
        //std::cout<<"This is the decay event: "<<i<<std::endl;
        std::string tempParentName = " ";
        float ExofFinaldecay = 0;
        //std::cout<<"tempParentName "<<std::endl;
        // COUT the de-excitation products -- stable
        // and its parent and excited energy of parent to decide the discrete levels
        for (;;)
        {
            CNucleus *parent;
            parent = products->getParent();
            if(parent->getName() != tempParentName){// avoid the same parent
                j++;
                if(products->getName() == parent->getName()){
                    productNamesL.push_back(" ");
                    productKEsL.push_back(0);
                    productExsL.push_back(0);
                    productPLx.push_back(0);
                    productPLy.push_back(0);
                    productPLz.push_back(0);
                }else{
                    productNamesL.push_back(parent->getLightDaughter()->getName());
                    productKEsL.push_back(parent->getLightDaughter()->getKE());
                    productExsL.push_back(parent->getLightDaughter()->fEx);
                    float* tempL = parent->getLightDaughter()->getMomentumVector();
                    productPLx.push_back(tempL[0]);
                    productPLy.push_back(tempL[1]);
                    productPLz.push_back(tempL[2]);
                }
                parentNames.push_back(parent->getName());
                productNamesH.push_back(parent->getHeavyDaughter()->getName()); 
                parentKEs.push_back(parent->getKE());                                                  
                parentExs.push_back(parent->fEx); 
                float* tempP = parent->getMomentumVector();
                parentPx.push_back(tempP[0]);
                parentPy.push_back(tempP[1]);
                parentPz.push_back(tempP[2]);                                         
                productKEsH.push_back(parent->getHeavyDaughter()->getKE());                                                                                            
                productExsH.push_back(parent->getHeavyDaughter()->fEx);
                float* tempH = parent->getHeavyDaughter()->getMomentumVector();
                productPHx.push_back(tempH[0]);
                productPHy.push_back(tempH[1]);
                productPHz.push_back(tempH[2]);

                tempParentName = parent->getName();
            }
            //tree->Fill();
            products = CN.getProducts(); // go to the next product
            if (products == NULL) {break;}
        }
        //laod the information before modify the kinematics
        int loadtimes = productNamesL.size();
        std::string laodparentname[loadtimes];
        std::string laodproductnameL[loadtimes];
        std::string laodproductnameH[loadtimes];
        for (int i = 0; i < loadtimes; i++){
            laodparentname[i]      = parentNames[i];
            laodproductnameL[i]    = productNamesL[i];
            laodproductnameH[i]    = productNamesH[i];
        }
        // std::cout<<"The total energy of the products: "<<sumE<<std::endl;
        // // add the discrete levels
        double finalLevel;
        std::string finalName = parentNames.back();
        std::string firstName = parentNames.front();
        int firstZ = Mass[firstName].Z;
        int firstA = Mass[firstName].A;
        float modifiedEx4Finaldecay = 0;
        // std::cout << "This has " << j << " times particle emission decay" << std::endl;
        // one thing needs to be attention, for 7Li decay in 15N, we do not have 7Li level in the table, so
        // we just pay  attention to the last decay, and modify the last decay of the residue of 7Li, 8Be
        // // not modify the 2-body decay
        if (j > 1 && firstZ == iZCN && firstA == iACN){
            finalLevel = productExsH.at(j - 2); // j for the decay times, j - 2 is the index of the last but second decay
            // std::cout<< "For the last decay, the excited energy of parent: "<< finalLevel<<std::endl;
            for (int i = 0; i < levels.size(); i++){
                if (levels[i].nucleusName == finalName){
                    //case 1 : small than all the level energies
                    if (finalLevel <= levels[i].levelEnergy[0]){
                        // std::cout<< "case 1"<<std::endl;
                        modifiedEx4Finaldecay = 0;
                        // no "last deacy" to modify 
                        productNamesL.pop_back();
                        productKEsL.pop_back();
                        productExsL.pop_back();
                        productPLx.pop_back();
                        productPLy.pop_back();
                        productPLz.pop_back();
                        parentNames.pop_back();
                        parentKEs.pop_back();
                        parentExs.pop_back();
                        parentPx.pop_back();
                        parentPy.pop_back();
                        parentPz.pop_back();
                        productNamesH.pop_back();
                        productKEsH.pop_back();
                        productExsH.pop_back(); 
                        productPHx.pop_back();
                        productPHy.pop_back();
                        productPHz.pop_back();
                        // modify the second to last decay
                        //(float parentMass, float parentKE, 
                        // float productMass, float residueMass, 
                        // float residueEx)
                        // std::cout<< parentM[rightIndex] << "" << parentKE[rightIndex] << " "
                        // << productM[rightIndex] << " " << Mass[finalName].mass << std::endl;
                        std::array<std::array<float, 4>, 2> result = modifyKinematic(
                            Mass[parentNames.back()].mass , parentExs.back(), parentKEs.back(),
                            Mass[productNamesL.back()].mass , Mass[productNamesH.back()].mass, 0
                        );
                        productPLx.back() = result[0][0];
                        productPLy.back() = result[0][1];
                        productPLz.back() = result[0][2];
                        productKEsL.back() = result[0][3];
                        productPHx.back() = result[1][0];
                        productPHy.back() = result[1][1];
                        productPHz.back() = result[1][2];
                        productKEsH.back() = result[1][3];
                        break;
                    } 
                    //case 2 : larger than all the level energies
                    // else if (finalLevel > levels[i].levelEnergy.back() + 1.5){
                    else if (finalLevel > levels[i].levelEnergy.back() + 2.){
                        // std::cout<< "Beyond all levels! Do not modify!" <<std::endl;
                        // gamma rays produced in decay, this is on-going
                        int numGammaRays = CN.getnGammaRays();
                        // sumGammaRayEnergy = CN.getSumGammaEnergy();
                        for (int kk = 0; kk < numGammaRays; kk++){
                            GammaRayEnergy.push_back(CN.getGammaRayEnergy(kk));
                            float emitTheta = acos(2 * randomDouble() - 1);
                            float emitPhi = 2 * CNucleus::pi * randomDouble();
                            GammaRayPx.push_back(GammaRayEnergy.back() * sin(emitTheta) * cos(emitPhi));
                            GammaRayPy.push_back(GammaRayEnergy.back() * sin(emitTheta) * sin(emitPhi));
                            GammaRayPz.push_back(GammaRayEnergy.back() * cos(emitTheta));
                        }
                        break;
                    }
                    //case 3 : between two level energies or others
                    else{
                        // std::cout<< "case 3"<<std::endl;
                        // preparing modifying the last decay
                        // cut all the final profucts
                        // decide which final level to decay
                        std::string tempProductName;
                        for (int jj = 0; jj < levels[i].levelEnergy.size(); ++jj){
                            float averageLevel = (levels[i].levelEnergy[jj] + levels[i].levelEnergy[jj + 1]) / 2;
                            if (finalLevel >= levels[i].levelEnergy[jj] && finalLevel < averageLevel){
                                modifiedEx4Finaldecay = levels[i].levelEnergy[jj];
                                tempProductName = levels[i].emitting[jj];
                            }else if(finalLevel > averageLevel && finalLevel <= levels[i].levelEnergy[jj + 1]){
                                // modifiedEx4Finaldecay = levels[i].levelEnergy[jj+1];
                                modifiedEx4Finaldecay = levels[i].levelEnergy[jj];
                                // tempProductName = levels[i].emitting[jj + 1];
                                tempProductName = levels[i].emitting[jj];
                            }
                            else if (finalLevel > levels[i].levelEnergy.back() && finalLevel <= levels[i].levelEnergy.back() + 2.){
                                modifiedEx4Finaldecay = levels[i].levelEnergy.back();
                                tempProductName = levels[i].emitting.back();
                            }
                        }
                        // modify second but last decay
                        productNamesL.pop_back();
                        productKEsL.pop_back();
                        productExsL.pop_back();
                        productPLx.pop_back();
                        productPLy.pop_back();
                        productPLz.pop_back();
                        parentNames.pop_back();
                        parentKEs.pop_back();
                        parentExs.pop_back();
                        parentPx.pop_back();
                        parentPy.pop_back();
                        parentPz.pop_back();
                        productNamesH.pop_back();
                        productKEsH.pop_back();
                        productExsH.pop_back(); 
                        productPHx.pop_back();
                        productPHy.pop_back();
                        productPHz.pop_back();
                        std::array<std::array<float, 4>, 2> result = modifyKinematic(
                            Mass[parentNames.back()].mass , parentExs.back(), parentKEs.back(),
                            Mass[productNamesL.back()].mass , Mass[productNamesH.back()].mass, modifiedEx4Finaldecay
                        );
                        productPLx.back() = result[0][0];
                        productPLy.back() = result[0][1];
                        productPLz.back() = result[0][2];
                        productKEsL.back() = result[0][3];
                        productPHx.back() = result[1][0];
                        productPHy.back() = result[1][1];
                        productPHz.back() = result[1][2];
                        productKEsH.back() = result[1][3];
                        //modify the last decay
                        if (tempProductName != "gamma"){
                            std::string residueName = getResidueName(productNamesH.back(), tempProductName);
                            result = modifyKinematic(
                                Mass[productNamesH.back()].mass, modifiedEx4Finaldecay, productKEsH.back(),
                                Mass[tempProductName].mass, Mass[residueName].mass, 0
                            );
                            productNamesL.push_back(tempProductName);
                            productKEsL.push_back(result[0][3]);
                            productExsL.push_back(0);
                            productPLx.push_back(result[0][0]);
                            productPLy.push_back(result[0][1]);
                            productPLz.push_back(result[0][2]);
                            parentNames.push_back(productNamesH.back());
                            parentKEs.push_back(productKEsH.back());
                            parentExs.push_back(modifiedEx4Finaldecay);
                            parentPx.push_back(productPHx.back());
                            parentPy.push_back(productPHy.back());
                            parentPz.push_back(productPHz.back());
                            productNamesH.push_back(residueName);
                            productKEsH.push_back(result[1][3]);
                            productExsH.push_back(0);
                            productPHx.push_back(result[1][0]);
                            productPHy.push_back(result[1][1]);
                            productPHz.push_back(result[1][2]);
                            if (residueName == "8Be"){
                                // std::cout<<"8Be decay"<<std::endl;
                                std::array<std::array<float, 4>, 2>resultMore = modifyKinematic(
                                    Mass["8Be"].mass, 0, result[1][3],
                                    Mass["4He"].mass, Mass["4He"].mass, 0
                                );
                                productNamesL.push_back("4He");
                                productKEsL.push_back(resultMore[0][3]);
                                productExsL.push_back(0);
                                productPLx.push_back(resultMore[0][0]);
                                productPLy.push_back(resultMore[0][1]);
                                productPLz.push_back(resultMore[0][2]);
                                parentNames.push_back(residueName);
                                parentKEs.push_back(result[1][3]);
                                parentExs.push_back(0);
                                parentPx.push_back(result[1][0]);
                                parentPy.push_back(result[1][1]);
                                parentPz.push_back(result[1][2]);
                                productNamesH.push_back("4He");
                                productKEsH.push_back(resultMore[1][3]);
                                productExsH.push_back(0);
                                productPHx.push_back(resultMore[1][0]);
                                productPHy.push_back(resultMore[1][1]);
                                productPHz.push_back(resultMore[1][2]);
                            } else if (residueName == "9B"){
                                // std::cout<<"9B decay + 8Be decay"<<std::endl;
                                std::array<std::array<float, 4>, 2>resultMore = modifyKinematic(
                                    Mass["9B"].mass, 0, result[1][3],
                                    Mass["p"].mass, Mass["8Be"].mass, 0
                                );
                                productNamesL.push_back("p");
                                productKEsL.push_back(resultMore[0][3]);
                                productExsL.push_back(0);
                                productPLx.push_back(resultMore[0][0]);
                                productPLy.push_back(resultMore[0][1]);
                                productPLz.push_back(resultMore[0][2]);
                                parentNames.push_back(residueName);
                                parentKEs.push_back(result[1][3]);
                                parentExs.push_back(0);
                                parentPx.push_back(result[1][0]);
                                parentPy.push_back(result[1][1]);
                                parentPz.push_back(result[1][2]);
                                productNamesH.push_back("8Be");
                                productKEsH.push_back(resultMore[1][3]);
                                productExsH.push_back(0);
                                productPHx.push_back(resultMore[1][0]);
                                productPHy.push_back(resultMore[1][1]);
                                productPHz.push_back(resultMore[1][2]);
                                std::array<std::array<float, 4>, 2>resultMore2 = modifyKinematic(
                                    Mass["8Be"].mass, 0, resultMore[1][3],
                                    Mass["4He"].mass, Mass["4He"].mass, 0
                                );
                                productNamesL.push_back("4He");
                                productKEsL.push_back(resultMore2[0][3]);
                                productExsL.push_back(0);
                                productPLx.push_back(resultMore2[0][0]);
                                productPLy.push_back(resultMore2[0][1]);
                                productPLz.push_back(resultMore2[0][2]);
                                parentNames.push_back("8Be");
                                parentKEs.push_back(resultMore[1][3]);
                                parentExs.push_back(0);
                                parentPx.push_back(resultMore[1][0]);
                                parentPy.push_back(resultMore[1][1]);
                                parentPz.push_back(resultMore[1][2]);
                                productNamesH.push_back("4He");
                                productKEsH.push_back(resultMore2[1][3]);
                                productExsH.push_back(0);
                                productPHx.push_back(resultMore2[1][0]);
                                productPHy.push_back(resultMore2[1][1]);
                                productPHz.push_back(resultMore2[1][2]);
                            } 
                        }
                        if (tempProductName == "gamma"){
                            GammaRayEnergy.push_back(modifiedEx4Finaldecay);
                            float emitTheta = acos(2 * randomDouble() - 1);
                            float emitPhi = 2 * CNucleus::pi * randomDouble();
                            GammaRayPx.push_back(GammaRayEnergy.back() * sin(emitTheta) * cos(emitPhi));
                            GammaRayPy.push_back(GammaRayEnergy.back() * sin(emitTheta) * sin(emitPhi));
                            GammaRayPz.push_back(GammaRayEnergy.back() * cos(emitTheta));
                        }
                        break;
                    }
                }
                else if (i == levels.size() - 1){
                    // std::cout << "No satisfied nucleus level!!!!"<<std::endl;
                    // gamma rays produced in decay, this is on-going
                    int numGammaRays = CN.getnGammaRays();
                    for (int kk = 0; kk < numGammaRays; kk++){
                        GammaRayEnergy.push_back(CN.getGammaRayEnergy(kk));
                        float emitTheta = acos(2 * randomDouble() - 1);
                        float emitPhi = 2 * CNucleus::pi * randomDouble();
                        GammaRayPx.push_back(GammaRayEnergy.back() * sin(emitTheta) * cos(emitPhi));
                        GammaRayPy.push_back(GammaRayEnergy.back() * sin(emitTheta) * sin(emitPhi));
                        GammaRayPz.push_back(GammaRayEnergy.back() * cos(emitTheta));
                    }
                    break;
                }
            }
        } else {
            // std::cout<< "Just one emission decay, no need to modify!!!"<<std::endl;
            int numGammaRays = CN.getnGammaRays();
            for (int kk = 0; kk < numGammaRays; kk++){
                GammaRayEnergy.push_back(CN.getGammaRayEnergy(kk));
                float emitTheta = acos(2 * randomDouble() - 1);
                float emitPhi = 2 * CNucleus::pi * randomDouble();
                GammaRayPx.push_back(GammaRayEnergy.back() * sin(emitTheta) * cos(emitPhi));
                GammaRayPy.push_back(GammaRayEnergy.back() * sin(emitTheta) * sin(emitPhi));
                GammaRayPz.push_back(GammaRayEnergy.back() * cos(emitTheta));
            }
            // break;
        }
        for (int p = 0; p < parentNames.size() - 1; p++){
            if (productNamesL[p] != " "){
                PDGid.push_back(1000000000 + Mass[productNamesL[p]].Z * 10000 + Mass[productNamesL[p]].A * 10);
                productPx.push_back(productPLx[p]);
                productPy.push_back(productPLy[p]);
                productPz.push_back(productPLz[p]);
                productKE.push_back(productKEsL[p]);
                productE.push_back(sqrt( Mass[productNamesL[p]].mass * Mass[productNamesL[p]].mass +
					productPLx[p] * productPLx[p] + productPLy[p]* productPLy[p] +
					productPLz[p] * productPLz[p]) );
            }
            if (productNamesH[p] != parentNames[p + 1]){
                PDGid.push_back(1000000000 + Mass[productNamesH[p]].Z * 10000 + Mass[productNamesH[p]].A * 10);
                productPx.push_back(productPHx[p]);
                productPy.push_back(productPHy[p]);
                productPz.push_back(productPHz[p]);
                productKE.push_back(productKEsH[p]);
                productE.push_back( sqrt(Mass[productNamesH[p]].mass * Mass[productNamesH[p]].mass + 
					productPHx[p] * productPHx[p] + productPHy[p] * productPHy[p] +
					productPHz[p] * productPHz[p]) );
            }
        }
        if (productNamesL[parentNames.size() - 1] != " "){
            PDGid.push_back(1000000000 + Mass[productNamesL.back()].Z * 10000 + Mass[productNamesL.back()].A * 10);
            productPx.push_back(productPLx.back());
            productPy.push_back(productPLy.back());
            productPz.push_back(productPLz.back());
            productKE.push_back(productKEsL.back());
	    productE.push_back(sqrt( Mass[productNamesL.back()].mass * Mass[productNamesL.back()].mass +                                                                                
                            productPLx.back() * productPLx.back() + productPLy.back()* productPLy.back() +
                            productPLz.back() * productPLz.back()) );
        }
        PDGid.push_back(1000000000 + Mass[productNamesH.back()].Z * 10000 + Mass[productNamesH.back()].A * 10);
        productPx.push_back(productPHx.back());
        productPy.push_back(productPHy.back());
        productPz.push_back(productPHz.back());
        productKE.push_back(productKEsH.back());
        productE.push_back( sqrt(Mass[productNamesH.back()].mass * Mass[productNamesH.back()].mass + 
                            productPHx.back() * productPHx.back() + productPHy.back() * productPHy.back() +          
                            productPHz.back() * productPHz.back()) );
	    nParticles = PDGid.size();
        // std::cout<< "For this decay products:"<<std::endl;
        for (int i = 0 ; i < nParticles; i++){
            // std::cout<< "PDGid: "<<PDGid[i]<<" Px: "<<productPx[i]<<" Py: "<<productPy[i]<<" Pz: "<<productPz[i]<<" KE: "<<productKE[i]<<" E: "<<productE[i]<<std::endl;
            if (PDGid[i] == 1000010010){
                PDGid[i] = 2212;
            }else if (PDGid[i] == 1000000010){
                PDGid[i] = 2112;
            }
            treePDGid[i] = PDGid[i];
            treePx[i] = productPx[i];
            treePy[i] = productPy[i];
            treePz[i] = productPz[i];
            treeE[i] = productE[i];
            treeKE[i] = productKE[i];
            // std::cout<< treePDGid[i] << std::endl;
        }
        for (int i = 0; i < GammaRayEnergy.size(); i++){
            PDGid.push_back(22);
            nParticles++;
            treePDGid[nParticles - 1] = 22;
            treePx[nParticles - 1] = GammaRayPx[i];
            treePy[nParticles - 1] = GammaRayPy[i];
            treePz[nParticles - 1] = GammaRayPz[i];
            treeE[nParticles - 1] = GammaRayEnergy[i];
            treeKE[nParticles - 1] = GammaRayEnergy[i];
        }
        // std::cout << "Before Fill:" << std::endl;
        int Esum = 0;
        for (int j = 0; j < nParticles; j++) {
            if (modeSelect == "--event"){
                std::cout << "PDGid: " << treePDGid[j] << " Px: " << treePx[j] << " Py: " << treePy[j] << " Pz: " << treePz[j] << " KE: " << treeKE[j] << " E: " << treeE[j] << std::endl;
            }
            // std::cout << "PDGid: " << treePDGid[j] << " Px: " << treePx[j] << " Py: " << treePy[j] << " Pz: " << treePz[j] << " KE: " << treeKE[j] << " E: " << treeE[j] << std::endl;
            Esum += treeE[j];
        }
        h1->Fill(Esum);
        // if (Esum < 10000){
        //     printf("The sum of the final products energy is: %d\n", Esum);
        //     for (int j = 0; j < nParticles; j++) {
        //         std::cout << "PDGid: " << treePDGid[j] << " Px: " << treePx[j] << " Py: " << treePy[j] << " Pz: " << treePz[j] << " KE: " << treeKE[j] << " E: " << treeE[j] << std::endl;
        //     }
        //     std::cout<<"Before modification: "<<std::endl;
        //     for (int j = 0 ; j < loadtimes; j++){
        //         std::cout<< "loadparentname: "<<laodparentname[j]<<" loadproductnameL: "<<laodproductnameL[j]<<" loadproductnameH: "<<laodproductnameH[j]<<std::endl;
        //     }
        // }
        if (modeSelect == "--batch"){
            tree->Fill();
        }
        // tree->Fill();
        CN.reset();
    }
        // rootFile->Write();
        // rootFile->Close();
        if (modeSelect == "--batch"){
            rootFile->Write();
            rootFile->Close();
            std::cout << "The mean energy of the final products is: " << h1->GetMean() << std::endl;
        }
}
