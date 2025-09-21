#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include <cmath>
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"
#include "TPaveText.h"

using namespace std;

class Generator
{
    public:
        string genName;

        TFile *input;
        TTree *tree;

        int b_max;
        double sigmaTotal, sigmaTotal_Prev;

        void GetGraphics(string genName, TFile& output);
        void GetProjections(TFile& output, string genName, TH2D& Graph);
        void GetTH1(string GraphName, TFile& output, string log);
        void GetGraphs(string er_GraphName, TFile& output);
        void GetTH2(string GraphName, TFile& output);
};

void Generator::GetGraphics(string genName, TFile& output)
{
    int Mch_max, Pt_max, Pt2_max, Pt24_max, Pt46_max, E_max;
    double b,p0,px,py,pz,m,chg,P_T=0, P_T_2=0, P_T_24=0, P_T_46=0, b_old=0, ityp, sigma, sigmaAvg=0, E_sum=0;
    int counter=0, nEvents=1, evNumber=0, error=0, Mch2=0, Mch46=0;

    if (genName == "Pythia")
    {
        Mch_max=200, b_max=7, Pt_max=80, Pt2_max=8, Pt24_max=25, Pt46_max=15, E_max=30;
    }
    else if(genName == "UrQMD")
    {
        Mch_max=200, b_max=7, Pt_max=80, Pt2_max=8, Pt24_max=25, Pt46_max=15, E_max=300;
    }

    input = new TFile(Form("%s_input.root", genName.c_str()), "read");
    tree = (TTree*)input->Get("particles");

    tree->SetBranchAddress("evNumber", &evNumber);
    tree->SetBranchAddress("b", &b);
    tree->SetBranchAddress("p0", &p0);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("m", &m);
    tree->SetBranchAddress("chg", &chg);
    tree->SetBranchAddress("ityp", &ityp);
    tree->SetBranchAddress("sigma", &sigma);

    int nEntries = tree->GetEntries(); 
    TH1D Eta_b[7];
    for (int i=0; i<7; i++)
    {
        Eta_b[i] = TH1D(Form("%s_Eta_%i", genName.c_str(), i), Form("%s_Eta %i fm", genName.c_str(), i+1), 100, -10, 10);
    }

    TH1D Eta_total(Form("%s_Eta_total", genName.c_str()), "Eta_total", 100, -10, 10);
    TH1D Teta_total(Form("%s_Teta_total", genName.c_str()), "Teta_total", 180, 0, 180);
    TH1D Eta_Straw(Form("%s_Eta_Straw", genName.c_str()),"Eta_Straw", 100, -3, 3);

    TH1D Mch(Form("%s_Mch", genName.c_str()), "Mch", Mch_max, 0, Mch_max);
    TH2D Mch_b(Form("%s_Mch_b", genName.c_str()), "Mch_b", Mch_max, 0, Mch_max, 50*b_max, 0, b_max);

    TH1D Espec(Form("%s_Espec", genName.c_str()), "Espec", 10*E_max, 0, E_max);
    TH2D Espec_b(Form("%s_Espec_b", genName.c_str()), "Espec_b", 10*E_max, 0, E_max, 50*b_max, 0, b_max);

    TH1D Pt_total(Form("%s_Pt_total", genName.c_str()), "Pt_total", Pt_max*10, 0, Pt_max);
    TH2D Pt_b(Form("%s_Pt_b", genName.c_str()), "Pt_b", Pt_max*10, 0, Pt_max, 50*b_max, 0, b_max);

    TH1D Pt_2(Form("%s_Pt_2", genName.c_str()), "Pt_2", 10*Pt2_max, 0, Pt2_max);
    TH2D Pt_b_2(Form("%s_Pt_b_2", genName.c_str()), "Pt_b_2", 10*Pt2_max, 0, Pt2_max, 50*b_max, 0, b_max);

    TH1D Pt_24(Form("%s_Pt_24", genName.c_str()), "Pt_24", 10*Pt24_max, 0, Pt24_max);
    TH2D Pt_b_24(Form("%s_Pt_b_24", genName.c_str()), "Pt_b_24", 10*Pt24_max, 0, Pt24_max, 50*b_max, 0, b_max);

    TH1D Pt_46(Form("%s_Pt_46", genName.c_str()), "Pt_46", 10*Pt46_max, 0, Pt46_max);
    TH2D Pt_b_46(Form("%s_Pt_b_46", genName.c_str()), "Pt_b_46", 10*Pt46_max, 0, Pt46_max, 50*b_max, 0, b_max);

    bool FirstParticle = true;
    array b_counter = {0, 0, 0, 0, 0, 0, 0};
    double P_module, eta, teta, p_t;

    ifstream fin("output.txt");
    ofstream csv_fout("ml/train.csv");
    csv_fout << "Multiplicity,Summarized Transverse Momentum,Spectators Energy Deposition,Multiplicity of soft particles,Multiplicity of hard particles,Impact Parameter" << endl;

    for (int j=0; j<nEntries; j++)
    {
        tree->GetEntry(j);

        if (evNumber == nEvents)
        {
            b_old=b;
            P_module = sqrt(px*px+py*py+pz*pz);
            eta = atanh(pz/P_module);
            teta = acos(pz/P_module)*180/M_PI;
            p_t = sqrt(px*px+py*py);

            if (!fin)
            {
                if (FirstParticle)
                {
                    sigmaAvg+=sigma;
                    FirstParticle = false;
                    for (int i=0; i<b_counter.size(); i++)
                    {
                        if (i == 0)
                        {
                            if (b_old >= double(i) && b_old <= double(i+1))
                            {
                                b_counter[i]++;
                            }
                        }
                        else
                        {
                            if (b_old > double(i) && b_old <= double(i+1))
                            {
                                b_counter[i]++;
                            }
                        }
                    }
                }

                for (int i=0; i<7; i++)
                {
                    if ((genName == "Pythia" && (ityp == 2212 || ityp == 2112)) || (genName == "UrQMD" && ityp == 1))
                    {
                        if (i == 0)
                        {
                            if (b_old >= double(i) && b_old <= double(i+1))
                            {
                                Eta_b[i].Fill(eta);
                            }
                        }
                        else
                        {
                            if (b_old > double(i) && b_old <= double(i+1))
                            {
                                Eta_b[i].Fill(eta);
                            }
                        }
                    }
                }

                // if (genName == "Pythia" && ityp!=2212 && chg!=0 && abs(eta)<2.266)
                // {
                    Eta_total.Fill(eta);
                    Teta_total.Fill(teta);
                // }
                // else if (genName == "UrQMD" && ityp!=1 && chg!=0 && abs(eta)<2.266)
                // {
                //     Eta_total.Fill(eta);
                //     Teta_total.Fill(teta);
                // }
            }

            if (chg !=0 && abs(eta)<2.266)
            {
                counter++;
                P_T += p_t;
                if (!fin)
                {
                    if (p_t<0.2)
                    {
                        P_T_2 += p_t;
                        Mch2++;
                    }

                    if (p_t>0.2 && p_t<0.4)
                    {
                        P_T_24 += p_t;
                    }

                    if (p_t>0.4 && p_t<0.6)
                    {
                        P_T_46 += p_t;
                        Mch46++;
                    }
                    Eta_Straw.Fill(eta);
                }
            }

            if (abs(eta)>3.5)
            {
                E_sum+=p0;
                // if (nEvents<=100) cout << genName << " " << nEvents << " " << b_old << " " << p0 << " " << E_sum << endl;
            }
        }
        else 
        {
            if(!fin)
            {
                FirstParticle = true;
                
                if (counter>1)
                {
                    Mch_b.Fill(counter, b_old);
                    Mch.Fill(counter);
                }
                if(P_T_2!=0)
                {
                    Pt_2.Fill(P_T_2);
                    Pt_b_2.Fill(P_T_2, b_old);
                }

                if(P_T_24!=0)
                {
                    Pt_24.Fill(P_T_24);
                    Pt_b_24.Fill(P_T_24, b_old);
                }

                if(P_T_46!=0)
                {
                    Pt_46.Fill(P_T_46);
                    Pt_b_46.Fill(P_T_46, b_old);
                }

                if(P_T!=0)
                {
                    Pt_total.Fill(P_T);
                    Pt_b.Fill(P_T, b_old);
                }

                if(E_sum!=0)
                {
                    Espec.Fill(E_sum);
                    Espec_b.Fill(E_sum, b_old);
                    if (nEvents<=100)
                    {
                        // cout << genName << " " << nEvents << " " << E_sum << endl;
                    }
                }
   
                P_T_2=0;
                P_T_24=0;
                P_T_46=0;
            }
            else
            {
                ifstream fin("output.txt");
                string line;
                double params[4][6] = {{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0}};
                string names[2];
                double sigmaP, sigmaU, a1, b1, c1, d1;
                if(fin.is_open())
                {
                    for (int i=0; i<4; i++)
                    {
                        getline(fin, line);
                        size_t colon_pos = line.find(':');
                        if (colon_pos != string::npos) 
                        {
                            if (i == 0) names[0] = line.substr(0, colon_pos);
                            if (i == 2) names[1] = line.substr(0, colon_pos);
                            line = line.substr(colon_pos + 1);
                            
                            std::istringstream iss(line);
                            iss >> sigmaP >> a1 >> b1 >> sigmaU >> c1 >> d1;
                            params[i][0] = sigmaP;
                            params[i][1] = a1;
                            params[i][2] = b1;
                            params[i][3] = sigmaU;
                            params[i][4] = c1;
                            params[i][5] = d1;
                        }
                    }
                }
                fin.close();

                int integralMch=0, integralPt=0;
                for (int j=0; j<2; j++)
                {
                    ofstream fout(Form("%s_centr_%s.txt", genName.c_str(), names[j].c_str()), ios::app);
                    
                    if (j==0)
                    {
                        if (genName == "Pythia")
                        {
                            TF1 *Pythia_centr_Mch= new TF1("Pythia_centr_Mch", "[0]*x+[1]", -params[0][2]/params[0][1], counter);
                            Pythia_centr_Mch->FixParameter(0, params[1][1]);
                            Pythia_centr_Mch->FixParameter(1, params[1][2]);
                            integralMch = int(100*params[0][1]/params[0][0]*Pythia_centr_Mch->Integral(-params[0][2]/params[0][1], counter));
                            if (integralMch <=0) integralMch = 0;
                            if (integralMch >=100) integralMch = 99;
                        }
                        else
                        {
                            TF1 *UrQMD_centr_Mch = new TF1("UrQMD_centr_Mch", "[0]*x+[1]", -params[0][5]/params[0][4], counter);
                            UrQMD_centr_Mch->FixParameter(0, params[1][4]);
                            UrQMD_centr_Mch->FixParameter(1, params[1][5]);
                            integralMch = int(100*params[0][4]/params[0][3]*UrQMD_centr_Mch->Integral(-params[0][5]/params[0][4], counter));
                            if (integralMch <=0) integralMch = 0;
                            if (integralMch >=100) integralMch = 99;
                        }
                        fout << Form("%s: #%i Mch: %i Centrality: %i Class: %i", genName.c_str(), nEvents, counter, integralMch, int(integralMch/20 +1)) << endl;
                    }
                    else
                    {
                        if (genName == "Pythia")
                        {
                            TF1 *Pythia_centr_Pt= new TF1("Pythia_centr_Pt", "[0]*x+[1]", -params[2][2]/params[2][1], P_T);
                            Pythia_centr_Pt->FixParameter(0, params[3][1]);
                            Pythia_centr_Pt->FixParameter(1, params[3][2]);
                            integralPt = int(100*params[2][1]/params[2][0]*Pythia_centr_Pt->Integral(-params[2][2]/params[2][1], P_T));
                            if (integralPt <=0) integralPt = 0;
                            if (integralPt >=100) integralPt = 99; 
                        }
                        else
                        {
                            TF1 *UrQMD_centr_Pt = new TF1("UrQMD_centr_Pt", "[0]*x+[1]", -params[2][5]/params[2][4], P_T);
                            UrQMD_centr_Pt->FixParameter(0, params[3][4]);
                            UrQMD_centr_Pt->FixParameter(1, params[3][5]);
                            integralPt = int(100*params[2][4]/params[2][3]*UrQMD_centr_Pt->Integral(-params[2][5]/params[2][4], P_T));
                            if (integralPt <=0) integralPt = 0;
                            if (integralPt >=100) integralPt = 99;
                        }
                        fout << Form("%s: #%i P_T: %f Centrality: %i Class: %i", genName.c_str(), nEvents, P_T, integralPt, int(integralPt/20 +1)) << endl;
                    }
                }
                if (int(integralMch/20 +1) != int(integralPt/20 +1))
                {
                    error++;
                    // cout << genName << " " << evNumber << " " << int(integralMch/20 +1) << " " << int(integralPt/20 +1) << endl;
                }
            }
            j--;
            nEvents++;
            if (int(b_old/(7/3))<2) csv_fout << counter << "," << P_T << "," << E_sum << "," << Mch2 << "," << Mch46 << "," << int(b_old/(7/3)) << endl;
            else csv_fout << counter << "," << P_T << "," << E_sum << "," << Mch2 << "," << Mch46 << "," << 2 << endl;
            
            Mch2=0;
            Mch46=0;
            P_T=0;
            counter=0;
            E_sum=0;
        }
    }

    if (fin)
    {
        cout << genName << " " << 100*error/double(evNumber) << endl;
    }
    else
    {
        sigmaTotal = sigmaAvg/(evNumber*10);

        GetProjections(output, genName, Mch_b);

        TGraphErrors *Hist1 = (TGraphErrors*)output.Get(Form("%s_Mch_b_Graph", genName.c_str()));
        TGraphErrors *Hist2 = new TGraphErrors();
        for (int i=0; i<6; i++)
        {
            Hist2->SetPoint(i, Hist1->GetPointX(i), sigmaTotal/evNumber*b_counter[i]);            
            Hist2->SetPointError(i, Hist1->GetErrorX(i), sqrt(sigmaTotal)/evNumber*b_counter[i]);
        }
        Hist2->SetName(Form("%s_Sigma_Mch_Graph", genName.c_str()));
        Hist2->Write();
        GetGraphs(Hist2->GetName(), output);

        GetProjections(output, genName, Pt_b);
        TGraphErrors *Hist3 = (TGraphErrors*)output.Get(Form("%s_Pt_b_Graph", genName.c_str()));
        TGraphErrors *Hist4 = new TGraphErrors();
        for (int i=0; i<6; i++)
        {
            Hist4->SetPoint(i, Hist3->GetPointX(i), sigmaTotal/evNumber*b_counter[i]);            
            Hist4->SetPointError(i, Hist3->GetErrorX(i), sqrt(sigmaTotal)/evNumber*b_counter[i]);
        }
        Hist4->SetName(Form("%s_Sigma_Pt_Graph", genName.c_str()));
        Hist4->Write();
        GetGraphs(Hist4->GetName(), output);

        GetProjections(output, genName, Pt_b_2);
        GetProjections(output, genName, Pt_b_24);
        GetProjections(output, genName, Pt_b_46);
        GetProjections(output, genName, Espec_b);

        for (int i=0; i<7; i++)
        {
            Eta_b[i].Write();
            GetTH1(Eta_b[i].GetName(), output, "log");
        }

        Mch.Write();
        GetTH1(Mch.GetName(), output, "log");
        Mch_b.Write();
        GetTH2(Mch_b.GetName(), output);

        Espec.Write();
        GetTH1(Espec.GetName(), output, "log");
        Espec_b.Write();
        GetTH2(Espec_b.GetName(), output);
        
        Pt_total.Write();
        GetTH1(Pt_total.GetName(), output, "log");
        Pt_b.Write();
        GetTH2(Pt_b.GetName(), output);
        Pt_2.Write();
        GetTH1(Pt_2.GetName(), output, "log");
        Pt_b_2.Write();
        GetTH2(Pt_b_2.GetName(), output);
        Pt_24.Write();
        GetTH1(Pt_24.GetName(), output, "log");
        Pt_b_24.Write();
        GetTH2(Pt_b_24.GetName(), output);
        Pt_46.Write();
        GetTH1(Pt_46.GetName(), output, "log");
        Pt_b_46.Write();
        GetTH2(Pt_b_46.GetName(), output);
        Teta_total.Write();
        GetTH1(Teta_total.GetName(), output, "log");
        Eta_total.Write();
        GetTH1(Eta_total.GetName(), output, "");
        Eta_Straw.Write();
        GetTH1(Eta_Straw.GetName(), output, "");

        input->Close();
    }
    
    csv_fout.close();
    };

void Generator::GetProjections(TFile& output, string genName, TH2D& Graph)
{
    output.cd();

    TH1D *projX[7];
    double Xmax=0;
    TGraphErrors *er_Graph = new TGraphErrors();

    for (int i=0; i<7; i++)
    {
        projX[i] = Graph.ProjectionX(Form("%s %s %i-%i fm", genName.c_str(), Graph.GetName(), i, i+1), i*50, (i+1)*50);

        er_Graph->SetPoint(i, projX[i]->GetMean(), i+0.5);
        Xmax=projX[i]->GetMean();
        er_Graph->SetPointError(i, projX[i]->GetStdDev(), 0.5);
        projX[i]->Write();
    }
    er_Graph->SetName(Form("%s_Graph", Graph.GetName()));
    er_Graph->Write();
    GetGraphs(er_Graph->GetName(), output);
};

void Generator::GetTH1(string GraphName, TFile& output, string log)
{
    gStyle->SetOptStat(kFALSE);

    string GenName, GenName_Prev="Pythia", HistName, ParamName, ExtraInfo;
    size_t pos = GraphName.find('_');
    if(pos!= string::npos)
    {
        GenName = GraphName.substr(0, pos);
        HistName = GraphName.substr(pos+1);
        size_t pos1 = HistName.find('_');
        if(pos1!= string::npos)
        {
            ExtraInfo = HistName.substr(pos1+1);
            ParamName = HistName.substr(0, pos1);
        }
        else
        {
            ParamName = HistName;
        }
    }

    if(GenName != "Pythia")
    {
        TCanvas *C = new TCanvas();
        if(log == "log")
        {
            C->SetLogy();
        }

        TH1D *Hist1 = (TH1D*)output.Get(Form("%s_%s", GenName.c_str(), HistName.c_str()));

        Hist1->GetXaxis()->SetTitle(ParamName.c_str());
        Hist1->SetLineColor(kRed);
        Hist1->SetLineWidth(2);
        Hist1->Draw();

        TH1D *Hist2 = (TH1D*)output.Get(Form("%s_%s", GenName_Prev.c_str(), HistName.c_str()));
        if (Hist2)
        {
            Hist2->SetLineColor(kBlue);
            Hist2->SetLineWidth(2);
            Hist2->Draw("same");
        }

        TLegend *legend = new TLegend(0.15, 0.15, 0.3, 0.25);
        if (Hist2) legend->AddEntry(Form("%s_%s", GenName_Prev.c_str(), HistName.c_str()), "Pythia", "l");
        legend->AddEntry(Form("%s_%s", GenName.c_str(), HistName.c_str()), GenName.c_str(), "l");
        legend->SetLineWidth(0);
        legend->SetTextSize(0.03);
        legend->Draw("same");

        TPaveText *pt = new TPaveText(0.65, 0.75, 0.89, 0.89, "NDC");  // Координаты NDC
        pt->AddText("Al+Al #sqrt{s_{NN}} = 10 GeV");
        pt->SetTextFont(42);
        pt->SetFillColor(kWhite);
        pt->SetTextSize(0.04);
        pt->SetBorderSize(0);
        // pt->Draw("same");

        C->SaveAs(Form("graphs/%s.pdf", HistName.c_str()));
    }
};

void Generator::GetGraphs(string er_GraphName, TFile& output)
{
    string GenName, GenName_Prev="Pythia", HistName, FirstParamName, SecondParamName, ExtraInfo;
    size_t pos = er_GraphName.find('_');
    if(pos!= string::npos)
    {
        GenName = er_GraphName.substr(0, pos);
        HistName = er_GraphName.substr(pos+1);
        size_t pos1 = HistName.find('_');
        if(pos1!= string::npos)
        {
            FirstParamName = HistName.substr(0, pos1);
            SecondParamName = HistName.substr(pos1+1);
            size_t pos2 = SecondParamName.find('_');
            if(pos2!= string::npos)
            {
                ExtraInfo = SecondParamName.substr(pos2+1);
                SecondParamName = SecondParamName.substr(0, pos2);
            }
        }
    }

    if(GenName == "Pythia")
    {
        sigmaTotal_Prev = sigmaTotal;
    }
    else
    {
        TCanvas *C = new TCanvas();
        TGraphErrors *Hist1 = (TGraphErrors*)output.Get(Form("%s_%s", GenName_Prev.c_str(), HistName.c_str()));
        if (Hist1)
        {
            double Xmax_Prev = Hist1->GetPointX(0);
            Hist1->SetTitle(HistName.c_str());
            if (FirstParamName == "Sigma")
            {
                Hist1->GetYaxis()->SetTitle(FirstParamName.c_str());
                Hist1->GetXaxis()->SetTitle(SecondParamName.c_str());
            }
            else
            {
                Hist1->GetXaxis()->SetTitle(FirstParamName.c_str());
                Hist1->GetYaxis()->SetTitle(SecondParamName.c_str());
            }
            Hist1->SetMarkerStyle(20);
            Hist1->SetMarkerColor(kBlue);
            Hist1->SetMarkerSize(1);
            Hist1->SetLineColor(kBlue);
            Hist1->SetLineWidth(1);
            Hist1->Draw("APE");

            TF1 *approxGraph1 = new TF1("PythiaApprox", "[0]*x + [1]", 0, Xmax_Prev);

            approxGraph1->SetParLimits(0, -10, 0.);
            approxGraph1->SetParLimits(1, 0, 1000);
            approxGraph1->SetLineColor(kBlue);
            // approxGraph1->SetLineStyle(7);
            approxGraph1->SetNpx(50);

            Hist1->Fit(approxGraph1, "ROB");
            approxGraph1->Draw("same");
        }

        TGraphErrors *Hist2 = (TGraphErrors*)output.Get(Form("%s_%s", GenName.c_str(), HistName.c_str()));
        double Xmax = Hist2->GetPointX(0);
        Hist2->SetTitle(HistName.c_str());
        Hist2->SetMarkerStyle(20);
        Hist2->SetMarkerColor(kRed);
        Hist2->SetMarkerSize(1);
        Hist2->SetLineColor(kRed);
        Hist2->SetLineWidth(1);
         if (Hist1) Hist2->Draw("PE same");
         else Hist2->Draw("APE");

        TF1 *approxGraph2 = new TF1("UrQMDApprox", "[0]*x + [1]", 0, Xmax);

        if (FirstParamName != "Espec") 
        {
            approxGraph2->SetParLimits(0, -1000, 0.);
            approxGraph2->SetParLimits(1, 0, 1000);
        }
        else
        {
            approxGraph2->SetParLimits(0, 0, 1000);
            approxGraph2->SetParLimits(1, -1000, 0);
        }
        approxGraph2->SetLineColor(kRed);
        // approxGraph2->SetLineStyle(7);
        approxGraph2->SetNpx(50);

        Hist2->Fit(approxGraph2, "ROB");
        approxGraph2->Draw("same");

        TLegend *legend;
        if (FirstParamName != "Espec") {legend = new TLegend(0.15, 0.15, 0.3, 0.35);}
        else {legend = new TLegend(0.15, 0.6, 0.3, 0.75);}
        if (Hist1) legend->AddEntry(Form("%s_%s", GenName_Prev.c_str(), HistName.c_str()), GenName_Prev.c_str(), "pl");
        legend->AddEntry(Form("%s_%s", GenName.c_str(), HistName.c_str()), GenName.c_str(), "pl");
        if (Hist1) legend->AddEntry(Form("%sApprox", GenName_Prev.c_str()), Form("%s_{approx}", GenName_Prev.c_str()), "l");
        legend->AddEntry(Form("%sApprox", GenName.c_str()), Form("%s_{approx}", GenName.c_str()), "l");
        legend->SetLineWidth(0);
        legend->SetTextSize(0.03);
        legend->Draw("same");

        TPaveText *pt;
        if (FirstParamName != "Espec") pt = new TPaveText(0.65, 0.75, 0.89, 0.89, "NDC");
        else pt = new TPaveText(0.15, 0.75, 0.39, 0.89, "NDC");
        pt->AddText("Al+Al #sqrt{s_{NN}} = 10 GeV");
        pt->SetTextFont(42);
        pt->SetFillColor(kWhite);
        pt->SetTextSize(0.04);
        pt->SetBorderSize(0);
        pt->Draw("same");

        ofstream fout("output.txt", ios::app);

        fout << Form("%s: %f %f %f", FirstParamName.c_str(), sigmaTotal, approxGraph2->GetParameter(0), approxGraph2->GetParameter(1)) << endl;

        fout.close();
        C->SaveAs(Form("graphs/%s.pdf", HistName.c_str()));
    }
};

void Generator::GetTH2(string GraphName, TFile& output)
{
    gStyle->SetOptStat(kFALSE);

    string GenName, HistName, FirstParamName, SecondParamName, ExtraInfo;
    size_t pos = GraphName.find('_');
    if(pos!= string::npos)
    {
        GenName = GraphName.substr(0, pos);
        HistName = GraphName.substr(pos+1);
        size_t pos1 = HistName.find('_');
        if(pos1!= string::npos)
        {
            FirstParamName = HistName.substr(0, pos1);
            SecondParamName = HistName.substr(pos1+1);
            size_t pos2 = SecondParamName.find('_');
            if(pos2!= string::npos)
            {
                ExtraInfo = SecondParamName.substr(pos2+1);
                SecondParamName = SecondParamName.substr(0, pos2);
            }
        }
    }

    TCanvas *C = new TCanvas();
    TH1D *Hist1 = (TH1D*)output.Get(Form("%s_%s", GenName.c_str(), HistName.c_str()));
    Hist1->GetXaxis()->SetTitle(FirstParamName.c_str());
    Hist1->GetYaxis()->SetTitle(SecondParamName.c_str());
    Hist1->Draw();

    TPaveText *pt = new TPaveText(0.65, 0.75, 0.89, 0.89, "NDC");  // Координаты NDC
    pt->AddText("Al+Al #sqrt{s_{NN}} = 10 GeV");
    pt->SetTextFont(42);
    pt->SetFillColor(kWhite);
    pt->SetTextSize(0.04);
    pt->SetBorderSize(0);
    pt->Draw("same");

    C->SaveAs(Form("graphs/%s.pdf", GraphName.c_str()));
};


void analyze()
{
    TFile output("output.root", "recreate");

    // Generator Pythia;
    // Pythia.GetGraphics("Pythia", output);

    Generator UrQMD;
    UrQMD.GetGraphics("UrQMD", output);

    output.Close();
}