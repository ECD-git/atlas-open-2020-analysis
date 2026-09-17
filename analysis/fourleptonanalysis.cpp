// An analysis program designed to search for the higgs boson through the decay of two Z bosons into 4 leptons

#include<filesystem>
#include<map>
#include<string>
#include <fstream>
#include <iostream>

#include "json.hpp"
#include<TFile.h>
#include<TTree.h>
#include<TH1F.h>
#include<Math/Vector4D.h>
#include<THStack.h>

// enter data directory path

std::filesystem::path FILEPATH = __FILE__;
std::filesystem::path DATADIRECTORY = FILEPATH.parent_path()/"4lep/Data/";
std::filesystem::path MCDIRECTORY = FILEPATH.parent_path()/"4lep/MC/";
std::filesystem::path testLepA = "data_A.4lep.root";
std::filesystem::path testMC = "Zee";
std::filesystem::path MCINFOPATH = FILEPATH.parent_path()/"mcinfofile.json";

// any data specific data

std::map<std::string, float> Luminosities; // in fb-1, int lumi for each data set
Luminosities["data_A"] = 0.5;
Luminosities["data_B"] = 1.9;
Luminosities["data_C"] = 2.9;
Luminosities["data_D"] = 4.7; // sum these for all 4 data sets total

// TODO: MAP OF DIFFERENT SIGNALS WE WANT TO ANALYSE

// MY OWN FOUR VECTOR STRUCT
// Deprecated in favour of ROOT 4vector library, but useful to see exactly how the calculation goes
struct FourVector
{
    double px,py,pz,E;
};
FourVector Make_Four_Vector(double pt, double eta, double phi, double E)
{
    FourVector v;
    v.px = pt * cos(phi);
    v.py = pt * sin(phi);
    v.pz = pt * sinh(eta);
    v.E = E;
    return v;
}
FourVector Sum_Four_Vector(const std::vector<FourVector> &vectors)
{
    FourVector sum = {0,0,0,0};
    for (auto &v : vectors)
    {
        sum.px += v.px;
        sum.py += v.py;
        sum.pz += v.pz;
        sum.E += v.E;
    }
    return sum;
}
double My_Calc_Invariant_Mass(const FourVector &vector)
{
    double pSq = vector.px*vector.px + vector.py*vector.py + vector.pz*vector.pz;
    // M^2 = E^2 - P^2
    double mSq = vector.E*vector.E - pSq;
    return sqrt(mSq);
}

// MC DATA FILE STRUCT
struct MCInfo 
{
    long DSID;
    double events;
    double red_eff;
    double sumw;
    double xsec;
};
std::map<std::string, MCInfo> Load_MC_Info(const std::filesystem::path &jsonPath)
{
    std::map<std::string, MCInfo> result;
    std::ifstream inFile(jsonPath);
    if (!inFile.is_open())
    {
        std::cout<<"Could not read json file at "<<jsonPath.string()<<std::endl;
        return result;
    }
    
    nlohmann::json j;
    inFile >> j;

    for(auto &[sampleName, fields]: j.items())
    {
        MCInfo info;
        info.DSID = fields.value("DSID", 0.0);
        info.events = fields.value("events", 0.0);
        info.red_eff = fields.value("red_eff", 1.0);
        info.sumw = fields.value("sumw", 0.0);
        info.xsec = fields.value("xsec", 0.0);
        result[sampleName] = info;
    }
    std::cout<<"SUCESSFULLY LOADED "<<jsonPath.string()<<std::endl;
    return result;
}

// ANALYSIS FUNCTIONS
bool Cut_Lep_Type(std::vector<unsigned int> *lep_type, bool print=false)
{
    // sum the lepton types for each entry, electron is 11, muon is 13,
    // cut off if the sum != 44,52 or 48 ie eeee,uuuu,eeuu
    int sumLep_Type = std::accumulate(lep_type->begin(),lep_type->end(),0);
    bool typeCutOff = (sumLep_Type != 44) & (sumLep_Type != 48) & (sumLep_Type != 52);
    if(print)
    {
        // print stuff out for debugging
        std::cout<<std::left<<std::setfill(' ')<<std::setw(7)<<"Types" << " = ";
        for(int i = 0; i<lep_type->size(); i++)
        {
            std::cout<<std::left<<std::setw(3)<<lep_type->at(i)<<" ";
        }
        std::cout<<std::left<<std::setfill(' ')<<std::setw(19)<<"| Sum Lepton Type"<<" = "<<std::setw(2)<<sumLep_Type<<" | Cut Off = "<<typeCutOff<<std::endl;
    }
    return typeCutOff;
}

bool Cut_Lep_Charge(std::vector<int> *lep_charge, bool print=false)
{
    // sum all lepton charges
    // cut off all lepton groups that dont sum to a charge of 0
    int sumLep_Charge = std::accumulate(lep_charge->begin(), lep_charge->end(), 0);
    bool chargeCutOff = (sumLep_Charge != 0);
    if(print)
    {
        // printing out
        std::cout<<std::left<<std::setfill(' ')<<std::setw(7)<<"Charges"<<" = ";
        for (int i=0; i<lep_charge->size(); i++)
        {
            std::cout<<std::left<<std::setw(3)<<lep_charge->at(i)<<" ";
        }
        std::cout<<std::left<<std::setw(19)<<"| Sum Lepton Charge"<<" = "<<std::setw(2)<<sumLep_Charge<<" | Cut Off = "<<chargeCutOff<<std::endl;
    }
    return chargeCutOff;
}

double Calc_Invariant_Mass(UInt_t lep_n, std::vector<float> *pt, std::vector<float> *eta, std::vector<float> *phi, std::vector<float> *E)
{
    // pt -> momentum perpendicular to beam
    // eta -> pseudorapidity -> angle rel to beam dir
    // phi -> azimuthal angle
    // E -> energy of lepton
    std::vector<ROOT::Math::PtEtaPhiEVector> leptonsFour;
    for(int i=0; i<lep_n; i++)
    {
        ROOT::Math::PtEtaPhiEVector temp(pt->at(i), eta->at(i), phi->at(i), E->at(i));
        leptonsFour.push_back(temp);
    }
    // sum the 4 vectors of the leptons and get invariant mass of resulting system
    ROOT::Math::PtEtaPhiEVector sumFour;
    for (auto &l : leptonsFour) sumFour += l;

    return sumFour.M()/1000; // MeV to GeV converstion here
}

// MAIN METHOD
void fourleptonanalysis() {
    // Read Info File for MC Sims
    auto mcInfoFile = Load_MC_Info(MCINFOPATH);

    // define histograms
    TH1F *h_mc_Zee = new TH1F("h_mc_Zee", "Four-lepton invariant mass; m_{4l} [GeV]; Events", 36, 80, 250);
    TH1F *h_mass_signal = new TH1F("h_mass_signal", "Four-lepton invariant mass; m_{4l} [GeV]; Events", 36, 80, 250);

    // OPEN REAL DATA FILE ---- TODO: LOOP HERE OVER testLepA -------- for path in real data paths or something
    TFile *file = TFile::Open((DATADIRECTORY/testLepA).string().c_str());
    TTree *tree = file->Get<TTree>("mini");
    // check for sucessful location
    if(!tree){
        std::cerr << "Could not find tree 'mini' in file." <<std::endl;
        return;
    } else {
        std::cout<<"SUCCESSFULLY READ file "<<(DATADIRECTORY/testLepA).string()<<std::endl;
    }

    // TRACK NEEDED VARS, lifted from the jupyter notebook from the open data release

    UInt_t lep_n;
    //std::vector<bool>    *lep_truthMatched = nullptr;
    //std::vector<bool>    *lep_trigMatched = nullptr;
    std::vector<float>   *lep_pt = nullptr;
    std::vector<float>   *lep_eta = nullptr;
    std::vector<float>   *lep_phi = nullptr;
    std::vector<float>   *lep_E = nullptr;
    //std::vector<float>   *lep_z0 = nullptr;
    std::vector<int>     *lep_charge = nullptr;
    std::vector<unsigned int> *lep_type = nullptr;
    //std::vector<bool>    *lep_isTightID = nullptr;
    //std::vector<float>   *lep_ptcone30 = nullptr;
    //std::vector<float>   *lep_etcone20 = nullptr;
    //std::vector<float>   *lep_trackd0pvunbiased = nullptr;
    //std::vector<float>   *lep_tracksigd0pvunbiased = nullptr;

    tree->SetBranchAddress("lep_n", &lep_n);
    //tree->SetBranchAddress("lep_truthMatched", &lep_truthMatched);
    //tree->SetBranchAddress("lep_trigMatched", &lep_trigMatched);
    tree->SetBranchAddress("lep_pt",&lep_pt);
    tree->SetBranchAddress("lep_eta", &lep_eta);
    tree->SetBranchAddress("lep_phi", &lep_phi);
    tree->SetBranchAddress("lep_E", &lep_E);
    //tree->SetBranchAddress("lep_z0", &lep_z0);
    tree->SetBranchAddress("lep_charge", &lep_charge);
    tree->SetBranchAddress("lep_type", &lep_type);
    //tree->SetBranchAddress("lep_isTightID", &lep_isTightID);
    //tree->SetBranchAddress("lep_ptcone30", &lep_ptcone30);
    //tree->SetBranchAddress("lep_etcone20", &lep_etcone20);
    //tree->SetBranchAddress("lep_trackd0pvunbiased", &lep_trackd0pvunbiased);
    //tree->SetBranchAddress("lep_tracksigd0pvunbiased", &lep_tracksigd0pvunbiased);

    Long64_t nEntries = tree->GetEntries(); // get number of entries, 39 for file A
    // im passing this as a long64_t since I imagine for full data sets the number of entries can excede the size of a 32 bit integer but its likely not needed for this exact use case
    // pass number of entires to console to check all is expected
    std::cout << "Number of Entires in Tree = " << nEntries << std::endl;     

    // ANALYSIS
    for(int i=0;i<nEntries;i++)
    {
        tree->GetEntry(i);
        // TODO:
        // Check for low transverse momentum, tight_ID and if lepton is isolated outside a jet

        if (lep_n != 4) continue; // we are interested only in 4 leptons
        // cut off entries without eeee, uuuu, or eeuu signals
        bool typeCutOff = Cut_Lep_Type(lep_type, false);
        // cut off entries with leptons that dont add up to 0 total charge
        bool chargeCutOff = Cut_Lep_Charge(lep_charge, false);
        // continue if we are skipping
        if (typeCutOff || chargeCutOff)
        {
            continue;
        }
        // calculate COM energy here using ROOT library
        double invarMass = Calc_Invariant_Mass(lep_n, lep_pt, lep_eta, lep_phi, lep_E);
        //std::cout<<"Invariant mass of leptons = "<<invarMass<<" GeV"<<std::endl;
        // save result to histogram
        h_mass_signal->Fill(invarMass);
        //std::cout<<std::endl;
    }
    // END OF LOOP HERE for real data -------------------

    // TODO, go through sim data, assigning a color to each background type
    // MC ANALYSIS --- loop testMC -----------------------

    auto it = mcInfoFile.find(testMC.string()); // find entry for dsid
    if (it != mcInfoFile.end())
    {
        const MCInfo &s = it->second;
        std::string mcFilePath = "mc_"+std::to_string(s.DSID)+'.'+testMC.string()+'.'+"4lep"+'.'+"root";

        TFile *mcFile = TFile::Open((MCDIRECTORY/mcFilePath).string().c_str());
        TTree *mcTree = mcFile->Get<TTree>("mini"); // open file
        // check for sucessful location
        if(!mcTree){
            std::cerr << "Could not find tree 'mini' in file." <<std::endl;
            return;
        } else {
            std::cout<<"SUCCESSFULLY READ file "<<(MCDIRECTORY/mcFilePath).string()<<std::endl;
        }
        mcTree->Print();

        Float_t mcWeight;
        Float_t scaleFactor_PILEUP;
        Float_t scaleFactor_ELE;
        Float_t scaleFactor_MUON;
        Float_t scaleFactor_LepTRIGGER;
        // TODO: im not sure the "mc" prefix is entirely necessary as at this point in the script the originals should be no longer needed but I REALLY dont want to be messing with memory allocation while im just figuring out how this works.
        UInt_t mclep_n;
        std::vector<int>     *mclep_charge = nullptr;
        std::vector<unsigned int> *mclep_type = nullptr;
        std::vector<float>   *mclep_pt = nullptr;
        std::vector<float>   *mclep_eta = nullptr;
        std::vector<float>   *mclep_phi = nullptr;
        std::vector<float>   *mclep_E = nullptr;
        
        mcTree->SetBranchAddress("lep_n", &mclep_n);
        mcTree->SetBranchAddress("mcWeight", &mcWeight);
        mcTree->SetBranchAddress("scaleFactor_PILEUP", &scaleFactor_PILEUP);
        mcTree->SetBranchAddress("scaleFactor_ELE", &scaleFactor_ELE);
        mcTree->SetBranchAddress("scaleFactor_MUON", &scaleFactor_MUON);
        mcTree->SetBranchAddress("scaleFactor_LepTRIGGER", &scaleFactor_LepTRIGGER);
        // data stuff
        mcTree->SetBranchAddress("lep_charge", &mclep_charge);
        mcTree->SetBranchAddress("lep_type", &mclep_type);
        mcTree->SetBranchAddress("lep_pt",&mclep_pt);
        mcTree->SetBranchAddress("lep_eta", &mclep_eta);
        mcTree->SetBranchAddress("lep_phi", &mclep_phi);
        mcTree->SetBranchAddress("lep_E", &mclep_E);

        Long64_t nMCEntries = mcTree->GetEntries();
        std::cout << "Number of Entires in MCTree = " << nMCEntries << std::endl;
    
        float xsec_weight = (Luminosities["data_A"]*1000*s.xsec)/(s.red_eff*s.sumw); // pb-1
        std::cout<<"xsec_weight for MC file = "<<xsec_weight<<std::endl;

        for (int i=0; i<nMCEntries; i++)
        {
            mcTree->GetEntry(i);
            // do the same data cutoffs as in the real data
            if (lep_n != 4) continue; // good just to check incase of errors
            bool typeCutOff = Cut_Lep_Type(mclep_type, false);
            bool chargeCutOff = Cut_Lep_Charge(mclep_charge, false);
            if (typeCutOff || chargeCutOff)
            {
                continue;
            }

            float total_weight = xsec_weight*mcWeight*scaleFactor_PILEUP*scaleFactor_ELE*scaleFactor_MUON*scaleFactor_LepTRIGGER;
            //std::cout<<"Total Weight for MC Event = "<<total_weight<<std::endl;

            // calc mass as before
            double mcInvarMass = Calc_Invariant_Mass(mclep_n, mclep_pt, mclep_eta, mclep_phi, mclep_E);
            //std::cout<<"Invariant mass of leptons = "<<mcInvarMass<<" GeV"<<std::endl;
            // save result to histogram
            h_mc_Zee->Fill(mcInvarMass, total_weight);
            //std::cout<<std::endl;
        }
    }   

    THStack *hs = new THStack("hs", "Four-lepton invariant mass; m_{4l} [GeV]; Events");

    h_mc_Zee->SetFillColor(kAzure - 9);
    h_mass_signal->SetFillColor(kRed - 7);

    hs->Add(h_mc_Zee);
    hs->Add(h_mass_signal);

    // DRAW
    // TODO stop this opening a window for some reason its a little annoying
    TCanvas *c1 = new TCanvas("c1", "c1");
    hs->Draw("HIST");
    h_mass_signal->Draw("E SAME");
    c1->SaveAs("four_lepton_mass.png");

    file->Close();
    std::cout<<"Sucessful execution."<<std::endl;
}
