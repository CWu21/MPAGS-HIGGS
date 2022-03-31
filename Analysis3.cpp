/* ----------------------------------------------
----------------HIGGS HOMEWORK ------------------
--------------------P3---------------------------
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

Double_t pi = TMath::Pi();
Double_t e = TMath::E();
Double_t alpha = (1E-11)*(0.511*0.511*1.16637*sqrt(2))/(4*pi);

Double_t lifetime(Double_t x){
    
    Double_t tau = 1.0/(0.5*alpha*x*pow((1-pow((2*0.511*0.001/x),2)),1.5));
    tau = tau*6.8*1E-25;
    return tau;
}

Double_t beta_gamma(Double_t tau){
    
    TF1 *f1 = new TF1("f1","lifetime(x)",0.001,0.06);
    Double_t m = f1->GetX(tau, 0.001, 0.006);
    Double_t beta_gamma = (sqrt(1.6*1.6-m*m))/m;

    return beta_gamma;
}


Double_t fiction1(Double_t x){

    Double_t f = 1 - pow(e, -2.0/(beta_gamma(x)*3*1E8*x));
    return f;
}

Double_t fiction2(Double_t x){

    Double_t bg =  (sqrt(1.6*1.6-x*x))/x;
    Double_t f = 1 - pow(e, -2.0/(bg*3*1E8*lifetime(x)));
    return f;
}

Double_t dsigma(Double_t mh, Double_t z){

    Int_t Z=82;

    Double_t f =  mh*mh*(1-z)/(0.000511*0.000511*z*z);
    Double_t F1 = log(184*pow(Z,-1/3));
    Double_t F2 = log(2*1.6*(1-z)/(0.000511*z*sqrt(1+f)))-0.5;

    Double_t dsigma = 9.405*1E-7*Z*Z*z*(1+f*2/3)*(F1)/pow((1+f),2);

    return dsigma;

}

void Drawlifetime(){

    TCanvas *c = new TCanvas();
    TF1 *f1 = new TF1("f1","lifetime(x)",0.001,0.06);

    f1->GetXaxis()->SetTitle("M_{H} (GeV)");
    f1->GetYaxis()->SetTitle("#tau_{H} (s)");
    f1->Draw();
    c->SetLogy();
    c->SaveAs("P3/1.png");
}

void fiction(){

    TCanvas *c = new TCanvas();
    TF1 *f1 = new TF1("f1","fiction1(x)",5*1E-11,1E-8);
    TF1 *f2 = new TF1("f2","fiction2(x)",0.001,0.05);

    f1->GetXaxis()->SetTitle("#tau_{H} (s)");
    f1->GetYaxis()->SetTitle("f");
    f1->Draw();
    c->SetLogx();
    c->SaveAs("P3/2.png");

    TCanvas *c1 = new TCanvas();
    f2->GetXaxis()->SetTitle("m_{H} (GeV)");
    f2->GetYaxis()->SetTitle("f");
    f2->Draw();
    c1->SaveAs("P3/3.png");
}

Double_t sigma(Double_t mh){

    Double_t sigma = 0;
    Double_t z = 0;
    Int_t n=1000;

    for(Int_t i=1;i<n;i++){

        z = i*1.0/n;
        sigma += dsigma(mh,z)*1.0/n;
    }


    return sigma;
}


void Crosssection(){
    
    TCanvas *c = new TCanvas();

/*    TF1 *f2 = new TF1("f2","dsigma(0.05, x)",0,1);

    f2->GetXaxis()->SetTitle("z");
    f2->GetYaxis()->SetTitle();
    f2->Draw();
    //c->SetLogy();
    c->SaveAs("P3/dsigma.png");*/
    TF1 *f1 = new TF1("f1","sigma(x)*0.389379*1E9",0.001,0.05);

    f1->GetXaxis()->SetTitle("m_{H} (GeV)");
    f1->GetYaxis()->SetTitle("#sigma (pb)");
    f1->Draw();
    c->SetLogy();
    c->SaveAs("P3/4.png");

}

void Analysis3(){
    
    //Drawlifetime();
    //fiction();
    Crosssection();


}
