#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"



using namespace std;

vector<vector<double>> extractAllImpactParameters(int& k, int startValue) 
{
    string filename = Form("./data/data%i.f14", k);

    vector<double> params;
    vector<vector<double>> results;
    int counter=0;
    streampos previous_pos;
    bool FirstParticle = true, Mch_checked, ComeBack = false;

    ifstream file(filename);

    string line;
    string prefix="impact_parameter_real/min/max(fm):";

    int Mch=0;
    double b_old=0, sigma_old;

    while (true) 
    {
        if (!getline(file, line)) 
        {
            break;
        }

        double b, min, max, sigma;
        
        if (line.find(prefix) == 0) 
        {
            size_t colon_pos = line.find(':');
            if (colon_pos != string::npos) 
            {
                string b_value = line.substr(colon_pos + 1);
                size_t colon_pos1 = b_value.find(':');
                if (colon_pos1 != string::npos) 
                {
                    string sigma_value = b_value.substr(colon_pos1 + 1);
                    b_value = b_value.substr(0, colon_pos1);
                    istringstream iss(sigma_value);
                    iss >> sigma;
                }


                istringstream iss1(b_value);
                iss1 >> b >> min >> max;
                
                Mch_checked = false;
                FirstParticle = true;

                if(Mch!=0)
                {
                    counter++;
                    // cout << startValue+counter << endl;
                    Mch=0;
                    file.seekg(previous_pos);
                    getline(file, line);
                    Mch_checked = true;
                    ComeBack = true;
                    continue;
                }

                if(!ComeBack)
                {
                    previous_pos = file.tellg();
                    b_old = b;
                    sigma_old = sigma;
                }
            }
        }
        else
        {
            double r0,rx,ry,rz,p0,px,py,pz,m,ityp,i3,chg,lcl,ncl,ora; 
            istringstream iss2(line);
            if (iss2 >> r0 >> rx >> ry >> rz >> p0 >> px >> py >> pz >> m >> ityp >> i3 >> chg >> lcl >> ncl >> ora)
            {
                if (!Mch_checked)
                {
                    if (chg!=0 && abs(atanh(pz/sqrt(px*px+py*py+pz*pz)))<2.266) 
                    {
                        Mch++;
                    }
                }
                else
                {

                    params.push_back(startValue+counter);
                    params.push_back(b_old);
                    params.push_back(sigma_old);
                    params.push_back(r0);
                    params.push_back(rx);
                    params.push_back(ry);
                    params.push_back(rz);
                    params.push_back(p0);
                    params.push_back(px);
                    params.push_back(py);
                    params.push_back(pz);
                    params.push_back(m);
                    params.push_back(ityp);
                    params.push_back(i3);
                    params.push_back(chg);
                    params.push_back(lcl);
                    params.push_back(ncl);
                    params.push_back(ora);

                    results.push_back(params);
                    params.clear();

                    ComeBack = false;
                    
                }
            }
        }
        if (startValue+counter == 100001)
        {
            cout << startValue+counter << " complete" << endl;
            break;
        } 
        // cout << startValue+counter << endl;
    }

    if (results.empty()) 
    {
        cerr << "Error: Could not find any lines starting with '" << prefix << "' in file " << filename << endl;
    }

    return results;
}  

void test() 
{
    TFile *output = new TFile("UrQMD_input.root", "recreate");
    TTree *tree = new TTree("particles", "particles");

    int nThreads = 5, startValue, evNumber;
    

    double b,r0,rx,ry,rz,p0,px,py,pz,m,ityp,i3,chg,lcl,ncl,ora,sigma;

    tree->Branch("evNumber", &evNumber, "evNumber/I");
    tree->Branch("b", &b, "b/D");
    tree->Branch("sigma", &sigma, "sigma/D");
    tree->Branch("r0", &r0, "r0/D");
    tree->Branch("rx", &rx, "rx/D");
    tree->Branch("ry", &ry, "ry/D");
    tree->Branch("rz", &rz, "rz/D");
    tree->Branch("p0", &p0, "p0/D");
    tree->Branch("px", &px, "px/D");
    tree->Branch("py", &py, "py/D");
    tree->Branch("pz", &pz, "pz/D");
    tree->Branch("m", &m, "m/D");
    tree->Branch("ityp", &ityp, "ityp/D");
    tree->Branch("i3", &i3, "i3/D");
    tree->Branch("chg", &chg, "chg/D");
    tree->Branch("lcl", &lcl, "lcl/D");
    tree->Branch("ncl", &ncl, "ncl/D");
    tree->Branch("ora", &ora, "ora/D");

    vector<vector<double>> real_params;

    for (int k=0; k<nThreads; k++)
    {
        if (k==0)
        {
            real_params = extractAllImpactParameters(k, 0);
            startValue = real_params[real_params.size()-1][0];
            // cout << startValue << endl;
        }
        else
        {
            real_params = extractAllImpactParameters(k, startValue);
            startValue = real_params[real_params.size()-1][0];
        }
        

        int counter=0;

        for (int j = 0; j < real_params.size(); j++) 
        {   
            if(tree)
            {
                evNumber=int(real_params[j][0]);
                b=real_params[j][1];
                sigma=real_params[j][2];
                r0=real_params[j][3];
                rx=real_params[j][4];
                ry=real_params[j][5];
                rz=real_params[j][6];
                p0=real_params[j][7];
                px=real_params[j][8];
                py=real_params[j][9];
                pz=real_params[j][10];
                m=real_params[j][11];
                ityp=real_params[j][12];
                i3=real_params[j][13];
                chg=real_params[j][14];
                lcl=real_params[j][15];
                ncl=real_params[j][16];
                ora=real_params[j][17];

                tree->Fill(); 
            }
        }
        cout << real_params.size() << endl;
        real_params.clear();
    }

    tree->Write();
    output->Delete("particles;22");
    output->Close();
}