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
#include<TLegend.h>

// Variables for histograms
const double xmin = 80;
const double xmax = 250;
const int stepsize = 5;
const int nbinsx = (xmax - xmin)/5;

// enter data directory path

std::filesystem::path FILEPATH = __FILE__;
std::filesystem::path DATADIRECTORY = FILEPATH.parent_path()/"4lep/Data/";
std::filesystem::path MCDIRECTORY = FILEPATH.parent_path()/"4lep/MC/";
std::filesystem::path testLepA = "data_A.4lep.root";
std::filesystem::path testMC = "Zee";
std::string fileEnd = ".4lep.root";
std::filesystem::path MCINFOPATH = FILEPATH.parent_path()/"mcinfofile.json";

// any data specific data
float Lumi = 10; // fb-1, for sum of all data

// MAP OF DIFFERENT SIGNALS WE WANT TO ANALYSE
std::map<std::string, std::vector<std::string>> samples;
std::map<std::string, std::vector<std::string>> MCs;
std::map<std::string, Color_t *> colors;
samples["data"] = {"data_A", "data_B", "data_C", "data_D"};
samples["MC"] = {"Background Z,t#bar{t}", "Background ZZ^{*}", "Signal (m_{H} = 125 GeV)"};

MCs["Background Z,t#bar{t}"] = {"Zee","Zmumu","ttbar_lep"};
colors["Background Z,t#bar{t}"] = new Color_t(9);

MCs["Background ZZ^{*}"] = {"llll"};
colors["Background ZZ^{*}"] = new Color_t(2);

MCs["Signal (m_{H} = 125 GeV)"] = {"ggH125_ZZ4lep","VBFH125_ZZ4lep","WH125_ZZ4lep","ZH125_ZZ4lep"};
colors["Signal (m_{H} = 125 GeV)"] = new Color_t(4);

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
bool Cut_Lep_isTight(std::vector<bool> *lep_isTightID)
{
    if (std::find(lep_isTightID->begin(), lep_isTightID->end(), false) != lep_isTightID->end())
    {
        return true;
    } else {
        return false;
    }
}

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
    TH1F *h_data = nullptr;
    std::map<std::string, TH1F *> H_BACKGROUND; // need to map these so they can be accessed later
    H_BACKGROUND["Background Z,t#bar{t}"] = nullptr;
    H_BACKGROUND["Background ZZ^{*}"] = nullptr;
    H_BACKGROUND["Signal (m_{H} = 125 GeV)"] = nullptr;

    // ACTUAL DATA
    h_data = new TH1F("h_data", "Data; m_{4l} [GeV]; Events", nbinsx, xmin, xmax);
    h_data->SetFillColor(10);
    for (std::string &dataFile : samples["data"]) // loop through all files in data
    {
        TFile *file = TFile::Open((DATADIRECTORY/(dataFile+fileEnd)).string().c_str());
        TTree *tree = file->Get<TTree>("mini");
        // check for sucessful location
        if(!tree){
            std::cerr << "Could not find tree 'mini' in file." <<std::endl;
            return;
        } else {
            std::cout<<"SUCCESSFULLY READ file "<<(DATADIRECTORY/(dataFile+fileEnd)).string()<<std::endl;
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
        std::vector<bool>    *lep_isTightID = nullptr;
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
        tree->SetBranchAddress("lep_isTightID", &lep_isTightID);
        //tree->SetBranchAddress("lep_ptcone30", &lep_ptcone30);
        //tree->SetBranchAddress("lep_etcone20", &lep_etcone20);
        //tree->SetBranchAddress("lep_trackd0pvunbiased", &lep_trackd0pvunbiased);
        //tree->SetBranchAddress("lep_tracksigd0pvunbiased", &lep_tracksigd0pvunbiased);

        Long64_t nEntries = tree->GetEntries(); // get number of entries, 39 for file A
        std::cout << "Number of Entires in Tree = " << nEntries << std::endl;     

        // ANALYSIS
        for(int i=0;i<nEntries;i++)
        {
            tree->GetEntry(i);
            // TODO:
            // Check for low transverse momentum, tight_ID and if lepton is isolated outside a jet

            if (lep_n != 4) continue; // we are interested only in 4 leptons
            // cut off if lepton has a false tight ID
            //if (Cut_Lep_isTight(lep_isTightID)) continue;
            // cut off entries without eeee, uuuu, or eeuu signals
            if (Cut_Lep_Type(lep_type, false)) continue;
            // cut off entries with leptons that dont add up to 0 total charge
            if (Cut_Lep_Charge(lep_charge, false)) continue;

            // calculate COM energy here using ROOT library
            double invarMass = Calc_Invariant_Mass(lep_n, lep_pt, lep_eta, lep_phi, lep_E);
            // save result to histogram
            h_data->Fill(invarMass);
        }

        // close the current data file
        file->Close();
    }

    // SIMULATION DATA
    for (std::string &bgType : samples["MC"]) // loop over all types 
    {   
        // bgType is for assigning titles and color now
        H_BACKGROUND[bgType] = new TH1F(("h_"+bgType).c_str(), (bgType + "; m_{4l} [GeV]; Events").c_str(), nbinsx, xmin, xmax);
        H_BACKGROUND[bgType]->Sumw2(); // tells histogram to track sum of squared wieghts per bin
        H_BACKGROUND[bgType]->SetFillColor(*colors[bgType]);

        for (std::string &bgTypeFile : MCs[bgType]) // all input files per type
        {
            auto it = mcInfoFile.find(bgTypeFile); // find entry for dsid
            if (it != mcInfoFile.end())
            {
                const MCInfo &s = it->second;
                std::string mcFilePath = "mc_"+std::to_string(s.DSID)+'.'+bgTypeFile+fileEnd;

                TFile *mcFile = TFile::Open((MCDIRECTORY/mcFilePath).string().c_str());
                TTree *mcTree = mcFile->Get<TTree>("mini"); // open file
                // check for sucessful location
                if(!mcTree){
                    std::cerr << "Could not find tree 'mini' in file." <<std::endl;
                    return;
                } else {
                    std::cout<<"SUCCESSFULLY READ file "<<(MCDIRECTORY/mcFilePath).string()<<std::endl;
                }

                Float_t mcWeight;
                Float_t scaleFactor_PILEUP;
                Float_t scaleFactor_ELE;
                Float_t scaleFactor_MUON;
                Float_t scaleFactor_LepTRIGGER;
                // NOTE: im not sure the "mc" prefix is entirely necessary as at this point in the script the originals should be no longer needed but I REALLY dont want to be messing with memory allocation while im just figuring out how this works.
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
                mcTree->SetBranchAddress("lep_charge", &mclep_charge);
                mcTree->SetBranchAddress("lep_type", &mclep_type);
                mcTree->SetBranchAddress("lep_pt",&mclep_pt);
                mcTree->SetBranchAddress("lep_eta", &mclep_eta);
                mcTree->SetBranchAddress("lep_phi", &mclep_phi);
                mcTree->SetBranchAddress("lep_E", &mclep_E);

                Long64_t nMCEntries = mcTree->GetEntries();
                std::cout << "Number of Entires in MCTree = " << nMCEntries << std::endl;
            
                float xsec_weight = (Lumi*1000*s.xsec)/(s.red_eff*s.sumw); // pb-1
                std::cout<<"xsec_weight for MC file = "<<xsec_weight<<std::endl;

                for (int i=0; i<nMCEntries; i++)
                {
                    mcTree->GetEntry(i);
                    // do the same data cutoffs as in the real data
                    if (mclep_n != 4) continue; // good just to check incase of errors
                    bool typeCutOff = Cut_Lep_Type(mclep_type, false);
                    bool chargeCutOff = Cut_Lep_Charge(mclep_charge, false);
                    if (typeCutOff || chargeCutOff)
                    {
                        continue;
                    }
                    // TODO add in other stuff once done for real data ie truth matching.

                    float total_weight = xsec_weight*mcWeight*scaleFactor_PILEUP*scaleFactor_ELE*scaleFactor_MUON*scaleFactor_LepTRIGGER;
                    // calc mass as before
                    double mcInvarMass = Calc_Invariant_Mass(mclep_n, mclep_pt, mclep_eta, mclep_phi, mclep_E);

                    // save result to histogram
                    H_BACKGROUND[bgType]->Fill(mcInvarMass, total_weight);
                }
                mcFile->Close();
            }   
        }
    }

    // calculate signal significance = N_signal/root(N_otherbackground)
    // from visual inspection, Higgs signal is significnt between 115 GeV to 130 GeV (under these hist settings)
    float higher = 130;
    float lower = 115; // replace with function params
    int lowerbin = (int)((lower-xmin)/stepsize) + 1;
    int higherbin = (int)((higher-xmin)/stepsize + 0.5) + 1;
    float N_sig = 0;
    float N_bg = 0;
    for (std::string &bgtype : samples["MC"])
    {
        for (int i = lowerbin; i<higherbin; i++)
        {
            if (bgtype != "Signal (m_{H} = 125 GeV)") N_bg += H_BACKGROUND[bgtype]->At(i);
            N_sig += H_BACKGROUND[bgtype]->At(i);
        }
    }
    float stat_sig = N_sig / std::sqrt(N_bg + (0.3*std::pow(N_bg,2))); // extra term for unaccounted uncert
    std::cout<<"Statistical Significance = "<<stat_sig<<std::endl;

    // store total background stacked
    THStack *hs = new THStack("hs", "Four-lepton invariant mass; m_{4l} [GeV]; Events");
    // store background statistical uncertainties
    TH1F *h_mc_total = new TH1F("h_mc_total", "Syst Uncert; m_{4l} [GeV]; Events", nbinsx, xmin, xmax);

    h_data->SetFillColor(kRed - 7);
    h_data->SetMarkerStyle(20); 
    h_data->SetMarkerSize(1.0);
    h_data->SetMarkerColor(kBlack);
    h_data->SetLineColor(kBlack);

    for (std::string &bgType : samples["MC"]) // loop over all types to stack bgs
    {
        hs->Add(H_BACKGROUND[bgType]);
        h_mc_total->Add(H_BACKGROUND[bgType]);
    }

    // DRAWING 

    // ensure range of hist does not cut off any data
    double maxY = std::max(hs->GetMaximum(), h_data->GetMaximum());
    hs->SetMaximum(maxY * 1.2);
    // stat uncert styling
    h_mc_total->SetFillColor(kGray + 2);
    h_mc_total->SetFillStyle(3345);   // hatched pattern
    h_mc_total->SetMarkerSize(0); 

    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 800);

    hs->Draw("HIST");
    h_mc_total->Draw("E2 SAME");
    h_data->Draw("E SAME");

    TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->AddEntry(h_data, "data", "lp");
    leg->AddEntry(h_mc_total, "Stat. uncert", "f");
    for (std::string &bgType : samples["MC"]) // loop over all types 
    {
        leg->AddEntry(("h_"+bgType).c_str(), bgType.c_str(), "f");
    }
    leg->Draw();

    // draw text
    TLatex Text;
    Text.SetNDC();              // coordinates as fractions of the pad (0 to 1), not data units
    Text.SetTextSize(0.04);
    Text.DrawLatex(0.6, 0.65, "#sqrt{s} = 13 TeV, #int L = 10.0 fb^{-1}");
    
    c1->SaveAs("four_lepton_mass.png");
    std::cout<<"Sucessful execution."<<std::endl;
}
