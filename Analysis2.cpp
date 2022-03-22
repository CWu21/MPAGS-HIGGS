/* ----------------------------------------------
----------------HIGGS HOMEWORK ------------------
--------------------P2---------------------------
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

Double_t myFunc1(double x) {return -18.3*x+24*log(16.2*x+5.2)+8*log(2.1*x+0.9)-6.1-log(TMath::Factorial(24)*TMath::Factorial(8));}
Double_t myFunc2(double x) {return -2*myFunc1(x)-2*6.4757;}
Double_t myFunc3(double x) {return -5.1*x+24*log(0.9*x+5.2)+8*log(4.2*x+0.9)-6.1-log(TMath::Factorial(24)*TMath::Factorial(8));}
Double_t myFunc4(double x) {return -2*myFunc3(x)-2*16.356;}
Double_t myFunc5(double x, double y) {return -18.3*x-5.1*y-6.1+24*log(16.2*x+0.9*y+5.2)+8*log(2.1*x+4.2*y+0.9)-log(TMath::Factorial(24)*TMath::Factorial(8));}
Double_t myFunc6(double x, double y) {return -2*myFunc5(x,y)-2*6.4757;}

Double_t myFunc7(double x) {return -18.3*x+22.3*log(16.2*x+5.2)+7.2*log(2.1*x+0.9)-6.1;}
Double_t myFunc8(double x) {return -2*myFunc7(x)+2*myFunc7(1.26606);}
Double_t myFunc9(double x) {return -5.1*x+22.3*log(0.9*x+5.2)+7.2*log(4.2*x+0.9)-6.1;}
Double_t myFunc10(double x) {return -2*myFunc9(x)+2*myFunc9(2.70086);}
Double_t myFunc11(double x, double y) {return -18.3*x-5.1*y-6.1+22.3*log(16.2*x+0.9*y+5.2)+7.2*log(2.1*x+4.2*y+0.9);}
Double_t myFunc12(double x, double y) {return -2*myFunc11(x,y)+2*myFunc11(1,1);}


void ggF(){
        
    TCanvas *c = new TCanvas();
    std::cout << "--------------ggF----------------- "<< std::endl;

    TF1 *f1 = new TF1("f1","-myFunc1(x)",0.75,2.25);
    TF1 *f2 = new TF1("f2","myFunc2(x)",0.75,2.25);

    double x0 = f1->GetMinimumX(0.75,2.25);
    double ymin = -myFunc1(x0);
    std::cout << "x0 = " << x0<< ", ymin = " << ymin<< std::endl;
    double xmin = f1->GetX(ymin+0.5,0.75,x0);
    double xmax = f1->GetX(ymin+0.5,x0,2.25);
    std::cout << "xmin = " << xmin<< ", xmax = " << xmax<< std::endl;
    std::cout << "negative uncertainty = " << x0-xmin<< ", positive uncertainty = " << xmax-x0<< std::endl;


    f1->GetXaxis()->SetTitle("#mu_{ggF}");
    f1->GetYaxis()->SetTitle("-ln L_{ggF}");
    f1->Draw();
    c->SaveAs("P2/1.png");

    double xmin1 = f2->GetX(1,0.75,x0);
    double xmax1 = f2->GetX(1,x0,2.25);
    std::cout << "xmin = " << xmin1<< ", xmax = " << xmax1<< std::endl;

    f2->GetXaxis()->SetTitle("#mu_{ggF}");
    f2->GetYaxis()->SetTitle("-2ln #Lambda_{ggF}");
    f2->Draw();
    c->SaveAs("P2/2.png");
}

void VBF(){
        
    TCanvas *c = new TCanvas();
    std::cout << "--------------VBF----------------- "<< std::endl;

    TF1 *f3 = new TF1("f3","-myFunc3(x)",2,5);
    TF1 *f4 = new TF1("f4","myFunc4(x)",2,5);

    double x0 = f3->GetMinimumX(2,5);
    double ymin = -myFunc3(x0);
    std::cout << "x0 = " << x0<< ", ymin = " << ymin<< std::endl;
    double xmin = f3->GetX(ymin+0.5,2,x0);
    double xmax = f3->GetX(ymin+0.5,x0,5);
    std::cout << "xmin = " << xmin<< ", xmax = " << xmax<< std::endl;
    std::cout << "negative uncertainty = " << x0-xmin<< ", positive uncertainty = " << xmax-x0<< std::endl;


    f3->GetXaxis()->SetTitle("#mu_{VBF}");
    f3->GetYaxis()->SetTitle("-ln L_{VBF}");
    f3->Draw();
    c->SaveAs("P2/3.png");

    double xmin1 = f4->GetX(1,2,x0);
    double xmax1 = f4->GetX(1,x0,5);
    std::cout << "xmin = " << xmin1<< ", xmax = " << xmax1<< std::endl;

    f4->GetXaxis()->SetTitle("#mu_{VBF}");
    f4->GetYaxis()->SetTitle("-2ln #Lambda_{VBF}");
    f4->Draw();
    c->SaveAs("P2/4.png");
}

void ggFVBF(){
    
    TCanvas *c = new TCanvas();
    std::cout << "--------------ggF+VBF----------------- "<< std::endl;

    TF2 *f5 = new TF2("f5","-myFunc5(x,y)",0,3,-0.3,5);
    TF2 *f6 = new TF2("f6","myFunc6(x,y)",0,3,-0.3,5);

    double zmin = f5->GetMinimum();
    double x0, y0;
    f5->GetMinimumXY(x0, y0);
    std::cout << "zmin = " << zmin<< std::endl;
    std::cout << "x0 = " << x0<< ", y0 = " << y0<< std::endl;

    f5->GetXaxis()->SetTitle("#mu_{ggF}");
    f5->GetYaxis()->SetTitle("#mu_{VBF}");
    f5->SetTitle("-ln L_{ggF,VBF}");
    //f5->Draw("lego");
    f5->Draw("CONT4Z LIST");
    c->SaveAs("P2/5.png");

    Double_t contours[2];
    contours[0] = 1;
    contours[1] = 4;
    f6-> SetContour(2, contours);
    
    f6->GetXaxis()->SetTitle("#mu_{ggF}");
    f6->GetYaxis()->SetTitle("#mu_{VBF}");
    f6->SetTitle("-2ln #Lambda");
    f6->Draw("CONT 4 LIST");
    c->Update();
    //c->SaveAs("P2/6.png");

 /*   TH2D *h1 = new TH2D("h1","-2ln #Lambda",  100, 0, 3, 100, 0, 5);
    for (Int_t i = 0; i < 100; i++) {
        for(Int_t j = 0; j < 100; j++){
            Double_t x = i*3.0/100;
            Double_t y = j*5.0/100;
            h1->SetBinContent(i,j, myFunc6(x,y));
        }
    }
    
    h1-> SetContour(2, contours);
    gStyle->SetOptStat(0);
    h1->Draw("CONT Z LIST");
    c->Update();
*/


    TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    TList* contLevel = NULL;
    TGraph* curv     = NULL;
    TGraph* gc       = NULL;
 
    Int_t nGraphs    = 0;
    Int_t TotalConts = 0;
 
    if (conts == NULL){
        TotalConts = 0;
        std::cout<<"TotalConts : "<<TotalConts<<std::endl;
        return 0;
    } else {
        TotalConts = conts->GetSize();
        std::cout<<"TotalConts : "<<TotalConts<<std::endl;
    }
 
    TCanvas* c1 = new TCanvas("c1","Contour List");
    c1->SetTopMargin(0.15);
    TH2F *h2 = new TH2F("h2","-2ln #Lambda; #mu_{ggF}; #mu_{VBF}", 1, 0, 3, 1, -0.3, 5); 

    gStyle->SetOptStat(0);
    h2->Draw();
    Double_t xval0, yval0, zval0;
    TLegend *leg = new TLegend(0.8,0.9,0.95,0.7);

    for(int i = 0; i < TotalConts; i++){
        contLevel = (TList*)conts->At(i);
        nGraphs += contLevel->GetSize();
        zval0 = contours[i];

        curv = (TGraph*)contLevel->First();
        for(int j = 0; j < contLevel->GetSize(); j++){
        curv->GetPoint(0, xval0, yval0);
        if (zval0 == 1) {
            curv-> SetLineStyle(kSolid);
            //curv-> SetTitle("68% contour");
            leg->AddEntry(curv, "68% contour","l");
        }
        if (zval0 == 4) {
            curv-> SetLineStyle(kDashed);
            //curv-> SetTitle("95% contour");
            leg->AddEntry(curv, "95% contour","l");
        }
        gc = (TGraph*)curv->Clone();
        gc->Draw("L");
        curv = (TGraph*)contLevel->After(curv);
        }

    }

    auto m1 = new TMarker(1,1,5);
    m1->Draw();

    auto m2 = new TMarker(x0,y0,2);
    //m2->SetMarkerColor(2);
    m2->Draw();

    leg->AddEntry(m1, "SM prediction", "P");
    leg->AddEntry(m2, "Best Fit", "P");
    leg->Draw();



    c1->Update();
    //c1->BuildLegend(0.8,0.9,0.95,0.7);
    c1->SaveAs("P2/6.png");

}

void ggF2(){
        
    TCanvas *c = new TCanvas();
    std::cout << "--------------GGF----------------- "<< std::endl;

    TF1 *f1 = new TF1("f1","-myFunc7(x)",0.75,2.25);
    TF1 *f2 = new TF1("f2","myFunc8(x)",0.75,2.25);

    double x0 = f1->GetMinimumX(0.75,2.25);
    double ymin = -myFunc7(x0);
    std::cout << "x0 = " << x0<< ", ymin = " << ymin<< std::endl;
    double xmin1 = f2->GetX(1,0.75,x0);
    double xmax1 = f2->GetX(1,x0,2.25);
    std::cout << "xmin = " << xmin1<< ", xmax = " << xmax1<< std::endl;
    std::cout << "negative uncertainty = " << x0-xmin1<< ", positive uncertainty = " << xmax1-x0<< std::endl;

    f2->GetXaxis()->SetTitle("#mu_{ggF}");
    f2->GetYaxis()->SetTitle("-2ln #Lambda_{ggF}");
    f2->Draw();
    c->SaveAs("P2/7.png");
}

void VBF2(){
        
    TCanvas *c = new TCanvas();
    std::cout << "--------------VBF----------------- "<< std::endl;

    TF1 *f1 = new TF1("f1","-myFunc9(x)",1,5);
    TF1 *f2 = new TF1("f2","myFunc10(x)",1,5);

    double x0 = f1->GetMinimumX(0,5);
    double ymin = -myFunc9(x0);
    
    std::cout << "x0 = " << x0<< ", ymin = " << ymin<< std::endl;
    double xmin1 = f2->GetX(1,1,x0);
    double xmax1 = f2->GetX(1,x0,5);
    std::cout << "xmin = " << xmin1<< ", xmax = " << xmax1<< std::endl;
    std::cout << "negative uncertainty = " << x0-xmin1<< ", positive uncertainty = " << xmax1-x0<< std::endl;

    f2->GetXaxis()->SetTitle("#mu_{VBF}");
    f2->GetYaxis()->SetTitle("-2ln #Lambda_{VBF}");
    f2->Draw();
    c->SaveAs("P2/8.png");
}

void ggFVBF2(){
    
    TCanvas *c = new TCanvas();
    std::cout << "--------------ggF+VBF----------------- "<< std::endl;

    TF2 *f5 = new TF2("f5","-myFunc11(x,y)",0,3,-0.3,5);
    TF2 *f6 = new TF2("f6","myFunc12(x,y)",0,3,-0.3,5);

    double zmin = f5->GetMinimum();
    double x0, y0;
    f5->GetMinimumXY(x0, y0);
    std::cout << "zmin = " << zmin<< std::endl;
    std::cout << "x0 = " << x0<< ", y0 = " << y0<< std::endl;

    f5->GetXaxis()->SetTitle("#mu_{ggF}");
    f5->GetYaxis()->SetTitle("#mu_{VBF}");
    f5->SetTitle("-ln L_{ggF,VBF}");
    f5->Draw("lego");
    //f5->Draw("CONT4Z LIST");
    c->SaveAs("P2/10.png");

    Double_t contours[2];
    contours[0] = 1;
    contours[1] = 4;
    f6-> SetContour(2, contours);
    
    f6->GetXaxis()->SetTitle("#mu_{ggF}");
    f6->GetYaxis()->SetTitle("#mu_{VBF}");
    f6->SetTitle("-2ln #Lambda");
    f6->Draw("CONT 4 LIST");
    c->Update();

    TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    TList* contLevel = NULL;
    TGraph* curv     = NULL;
    TGraph* gc       = NULL;
 
    Int_t nGraphs    = 0;
    Int_t TotalConts = 0;
 
    if (conts == NULL){
        TotalConts = 0;
        std::cout<<"TotalConts : "<<TotalConts<<std::endl;
        return 0;
    } else {
        TotalConts = conts->GetSize();
        std::cout<<"TotalConts : "<<TotalConts<<std::endl;
    }
 
    TCanvas* c1 = new TCanvas("c1","Contour List");
    c1->SetTopMargin(0.15);
    TH2F *h2 = new TH2F("h2","-2ln #Lambda; #mu_{ggF}; #mu_{VBF}", 1, 0, 3, 1, -0.3, 5); 

    gStyle->SetOptStat(0);
    h2->Draw();
    Double_t xval0, yval0, zval0;
    TLegend *leg = new TLegend(0.8,0.9,0.95,0.7);

    for(int i = 0; i < TotalConts; i++){
        contLevel = (TList*)conts->At(i);
        nGraphs += contLevel->GetSize();
        zval0 = contours[i];

        curv = (TGraph*)contLevel->First();
        for(int j = 0; j < contLevel->GetSize(); j++){
        curv->GetPoint(0, xval0, yval0);
        if (zval0 == 1) {
            curv-> SetLineStyle(kSolid);
            //curv-> SetTitle("68% contour");
            leg->AddEntry(curv, "68% contour","l");
        }
        if (zval0 == 4) {
            curv-> SetLineStyle(kDashed);
            //curv-> SetTitle("95% contour");
            leg->AddEntry(curv, "95% contour","l");
        }
        gc = (TGraph*)curv->Clone();
        gc->Draw("L");
        curv = (TGraph*)contLevel->After(curv);
        }

    }

    auto m1 = new TMarker(1,1,5);
    m1->Draw();

    auto m2 = new TMarker(x0,y0,2);
    //m2->SetMarkerColor(2);
    m2->Draw();

    leg->AddEntry(m1, "SM prediction", "P");
    leg->AddEntry(m2, "Best Fit", "P");
    leg->Draw();



    c1->Update();
    c1->SaveAs("P2/9.png");

}


void Analysis2(){

    ggF();
    VBF();
    ggFVBF();
std::cout << "------------------------------------------------------ "<< std::endl;
std::cout << "--------------experiment design level----------------- "<< std::endl;
std::cout << "------------------------------------------------------ "<< std::endl;
    ggF2();
    VBF2();
    ggFVBF2();
    



}