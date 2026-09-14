// An analysis program designed to search for the higgs boson through the decay of two Z bosons into 4 leptons

#include<filesystem>
#include<TFile.h>
#include<TTree.h>
#include<TH1F.h>
#include<Math/Vector4D.h>

// enter data directory path
std::filesystem::path DATADIRECTORY = "/Users/ecd/Desktop/Academic Work/ATLAS-OPEN-2020-13TEV/4lep";
std::filesystem::path lepA = "Data/data_A.4lep.root";

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
    // OPEN FILE

    // start with the A real data file
    TFile *file = TFile::Open((DATADIRECTORY/lepA).string().c_str());
    TTree *tree = file->Get<TTree>("mini");
    // check for sucessful location
    if(!tree){
        std::cerr << "Could not find tree 'mini' in file." <<std::endl;
        return;
    } else {
        std::cout<<"SUCCESSFULLY READ file "<<(DATADIRECTORY/lepA).string()<<std::endl;
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

    TH1F *h_mass = new TH1F("h_mass", "Four-lepton invariant mass; m_{4l} [GeV]; Events", 36, 80, 250);

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
        // this is already the case for the data set im using but worth adding in incase i switch to others

        // cut off entries without eeee, uuuu, or eeuu signals
        bool typeCutOff = Cut_Lep_Type(lep_type, true);
        // cut off entries with leptons that dont add up to 0 total charge
        bool chargeCutOff = Cut_Lep_Charge(lep_charge, true);
        // continue if we are skipping
        if (typeCutOff || chargeCutOff)
        {
            std::cout<<std::endl;
            continue;
        }
        // calculate COM energy here using ROOT library
        double invarMass = Calc_Invariant_Mass(lep_n, lep_pt, lep_eta, lep_phi, lep_E);
        std::cout<<"Invariant mass of leptons = "<<invarMass<<" GeV"<<std::endl;

        // save result to histogram
        h_mass->Fill(invarMass);
        
        std::cout<<std::endl;
    }

    // lets draw the data
    // TODO stop this opening a window for some reason its a little annoying
    TCanvas *c1 = new TCanvas("c1", "c1");
    h_mass->Draw("E");
    c1->SaveAs("four_lepton_mass.png");

    file->Close();
    std::cout<<"Sucessful execution."<<std::endl;
}
