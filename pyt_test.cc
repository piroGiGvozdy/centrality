#include <iostream>
#include "TFile.h"
#include <cmath>
#include "TH1.h"
#include "TTree.h"

#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"

using namespace Pythia8;
using namespace std;

int main()
{
    int counter=0, evNumber=0;
    double b,p0,px,py,pz,m,sigma,ityp,chg;
    bool isFinalChecked;

    TFile *output = new TFile("Pythia_input.root", "recreate");
    TTree *tree = new TTree("particles", "particles");

    tree->Branch("evNumber", &evNumber, "evNumber/I");
    tree->Branch("b", &b, "b/D");
    tree->Branch("sigma", &sigma, "sigma/D");
    tree->Branch("p0", &p0, "p0/D");
    tree->Branch("px", &px, "px/D");
    tree->Branch("py", &py, "py/D");
    tree->Branch("pz", &pz, "pz/D");
    tree->Branch("m", &m, "m/D");
    tree->Branch("chg", &chg, "chg/D");
    tree->Branch("ityp", &ityp, "ityp/D");

    Pythia8::Pythia pythia;

    pythia.readString("Beams:idA = 1000130270");
    pythia.readString("Beams:idB = 1000130270");
    pythia.readString("Beams:eCM = 10");
    pythia.readString("SoftQCD:all = on");
    pythia.readString("HeavyIon:SigFitDefPar = 8.16,22.66,0.29");
    pythia.readString("HeavyIon:SigFitDefAvNDb = 0.61");
    pythia.readString("HeavyIon:SigFitNGen = 0");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 1000");

    pythia.init();

    while(counter!=100000)
    {
        if(!pythia.next()) continue;
        isFinalChecked = true;

        if (pythia.info.hiInfo->b()<=7)
        {
            int entries = pythia.event.size();

            for (int j=0; j<entries; j++)
            {
                if (pythia.event[j].isFinal())
                {
                    if (isFinalChecked)
                    {
                        counter++;
                        isFinalChecked = false;
                    }

                    evNumber=counter;
                    b=pythia.info.hiInfo->b();
                    sigma=pythia.info.sigmaGen();
                    p0=pythia.event[j].e();
                    px=pythia.event[j].px();
                    py=pythia.event[j].py();
                    pz=pythia.event[j].pz();
                    m=pythia.event[j].m();
                    chg=pythia.event[j].charge();
                    ityp=pythia.event[j].id();

                    tree->Fill();
                }
            }
        }
        else pythia.next();
    }

    output->Delete("particles;12");
    tree->Write();
    output->Close();

    return 0;
}