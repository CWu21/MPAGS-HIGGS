/* ----------------------------------------------
----------------HIGGS HOMEWORK ------------------
--------------------P4---------------------------
-------------------------------------------------*/

#include "TCanvas.h"
#include "TMath.h"
#include "TLatex.h"
#include "TF1.h"
#include "TGraph.h"
#include "TList.h"
#include "TObject.h"
#include "TAttLine.h"
#include <TMarker.h>
#include <vector>

TRandom3*rndm=new TRandom3(0);

vector<Double_t> c{ 1, 0.05429, 0.008939, 0.0000890, 0.000161, 1.070, 0.5256, 0.0678, 0.00179, 0.0000659, 0.0737, 114.9};
vector<Double_t> d{ 0.2312527, 0.0004729, 0.0000207, 0.00000385, -0.00000185, 0.0207, -0.002851, 0.000182, -0.00000976, 0.000398, -0.655};

Double_t dH(Double_t mH) {return log(mH/100);}
Double_t dh(Double_t mH) {return pow(mH/100,2);}
Double_t dt(Double_t mt) {return pow(mt/174.3,2)-1;}

Double_t LH(Double_t mH) {return log(mH/100);}
Double_t deH(Double_t mH) {return mH/100;}
Double_t det(Double_t mt) {return pow(mt/178.0,2)-1;}

Double_t mWini = 80.3799;

Double_t mW(Double_t mH, Double_t mt){

    Double_t m;
    m = mWini - c[1]*dH(mH) - c[2]*dH(mH)*dH(mH) +c[3]*dH(mH)*dH(mH)*dH(mH)*dH(mH) + c[4]*(dh(mH)-1) + c[6]*dt(mt) -c[7]*dt(mt)*dt(mt)-c[8]*dH(mH)*dt(mt) + c[9]*dh(mH)*dt(mt);

    return m;
}

Double_t sin2(Double_t mH, Double_t mt){

    Double_t s;
    s = d[0] + d[1]*LH(mH) + d[2]*LH(mH)*LH(mH) + d[3]*LH(mH)*LH(mH)*LH(mH)*LH(mH) + d[4]*(deH(mH)*deH(mH)-1) + d[6]*det(mt) + d[7]*dt(mt)*dt(mt) + d[8]*dt(mt)*(deH(mH)-1) + d[9]*(0.1192/0.117-1) + d[10]*(91.1875/91.1876-1);

    return s;
}

Double_t myFunc1(Double_t x) {return mW(x, 172.4);} 
Double_t myFunc2(Double_t x) {return sin2(x, 172.4);}

Double_t chi2mh(Double_t x){

    Double_t chi2 = 0;
    Double_t mmeas = 80.399;
    chi2 = pow((myFunc1(x)-mmeas),2)/pow(mmeas*0.0001,2);  //0.01, 0.001, 0.0001

    return chi2;
    
}

Double_t chi2mh_mt(Double_t x, Double_t y){
    
    Double_t chi2 = 0;
    Double_t mwmeas = 80.399; 
    Double_t mtmeas = 172.4;

    chi2 = pow((mW(x,y)-mwmeas),2)/pow(0.025,2) + pow(y-mtmeas,2)/pow(mtmeas*0.001,2); // 0.1, 0.01, 0.001

    return chi2;
}

Double_t chi2sin2(Double_t x){

    Double_t chi2 = 0;
    Double_t smeas = 0.2324;
    chi2 = pow((myFunc2(x)-smeas),2)/pow(smeas*0.00001,2);  //0.001, 0.0001, 0.00001

    return chi2;

}


void DrawMW(){
    TCanvas *c = new TCanvas();

    TF1 *f1 = new TF1("f1","myFunc1(x)",10,1000);

    f1->GetXaxis()->SetTitle("M_{H} (GeV)");
    f1->GetYaxis()->SetTitle("M_{W} (GeV)");
    f1->Draw();
    c->SaveAs("P4/1.png");

    TF1 *f2 = new TF1("f2","1.0*ROOT::Math::normal_pdf(x, 0.025, 80.399)", 80.2, 80.6);
    f2->GetXaxis()->SetTitle("M_{W} (GeV)");
    f2->SetTitle("M_{W}");
    f2->Draw();
    c->SaveAs("P4/mW.png");

}

void Drawsin2(){
    
    TCanvas *c = new TCanvas();
    TF1 *f1 = new TF1("f1","myFunc2(x)",10,1000);

    f1->GetXaxis()->SetTitle("M_{H} (GeV)");
    f1->GetYaxis()->SetTitle("sin^{2} #theta_{eff}^{l}");
    f1->SetTitle("sin^{2} #theta_{eff}^{l}");
    f1->Draw();
    c->SaveAs("P4/2.png");
}


void DrawChi2mh(){
    
    TCanvas *c = new TCanvas();
    TF1 *f1 = new TF1("f1","chi2mh(x)",10,1000);

    double x0 = f1->GetMinimumX(10,200);
    std::cout << "x0 = " << x0<< std::endl;
    double xmin = f1->GetX(1,10,x0);
    double xmax = f1->GetX(1,x0,1000);
    std::cout << "xmin = " << xmin<< ", xmax = " << xmax<< std::endl;

    f1->GetXaxis()->SetTitle("M_{H} (GeV)");
    f1->GetYaxis()->SetTitle("#chi^{2}");
    f1->Draw();
    c->SaveAs("P4/chi2mh.png");

}

void DrawChi2sin2(){
    
    TCanvas *c = new TCanvas();
    TF1 *f1 = new TF1("f1","chi2sin2(x)",10,1000);

    double x0 = f1->GetMinimumX(10,1000);
    std::cout << "x0 = " << x0<< std::endl;
    double xmin = f1->GetX(1,10,x0);
    double xmax = f1->GetX(1,x0,1000);
    std::cout << "xmin = " << xmin<< ", xmax = " << xmax<< std::endl;

    f1->GetXaxis()->SetTitle("M_{H} (GeV)");
    f1->GetYaxis()->SetTitle("#chi^{2}");
    f1->Draw();
    c->SaveAs("P4/chi2sin2.png");

}

void DrawChi2mh_mt(){
    
    TCanvas *c = new TCanvas();
    TF2 *f1 = new TF2("f1","chi2mh_mt(x,y)",10,1000,100,200);
    
    double x0, y0;
    f1->GetMinimumXY(x0, y0);
    std::cout << "x0 = " << x0<< ", y0 = " << y0<< std::endl;

    f1->GetXaxis()->SetTitle("M_{H} (GeV)");
    f1->GetYaxis()->SetTitle("m_{t} (GeV)");
    f1->SetTitle("#chi^{2}");
    f1->Draw("CONT4Z");
    c->SaveAs("P4/chi2mh_mt.png");

}

Double_t Combine(Double_t x, Double_t y){

    Double_t chi2 = 0;
    Double_t mwmeas = 80.399; 
    Double_t mtmeas = 172.4;
    Double_t smeas = 0.2324;

    chi2 = pow((mW(x,y)-mwmeas),2)/pow(0.025,2) + pow(y-mtmeas,2)/pow(mtmeas*0.01,2) + pow((sin2(x,y)-smeas),2)/pow(0.0012,2);

    return chi2;

}

void DrawCombine(){

    TCanvas *c = new TCanvas();
    TF2 *f1 = new TF2("f1","chi2mh_mt(x,172.4)",10,1000,100,200);
    
    double x0, y0;
    f1->GetMinimumXY(x0, y0);
    std::cout << "x0 = " << x0<< ", y0 = " << y0<< std::endl;
/*
    double x0 = f1->GetMinimumX(10,1000);
    std::cout << "x0 = " << x0<< std::endl;
    double xmin = f1->GetX(1,10,x0);
    double xmax = f1->GetX(1,x0,1000);
    std::cout << "xmin = " << xmin<< ", xmax = " << xmax<< std::endl;*/


    f1->GetXaxis()->SetTitle("M_{H} (GeV)");
    f1->GetYaxis()->SetTitle("m_{t} (GeV)");
    f1->SetTitle("#chi^{2}");
    f1->Draw("CONT4Z");
    c->SaveAs("P4/chi2combine.png");


}


void Analysis4(){

    //DrawMW();
    //Drawsin2();

    //DrawChi2mh();
    //DrawChi2mh_mt();

    //DrawChi2sin2();

    DrawCombine();




}
