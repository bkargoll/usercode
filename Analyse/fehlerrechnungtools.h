//#include <math.h>
#include "TMath.h"
#include <fehlerrechnung.h>
//print für fehlerrechnung

void cprint(char *txt, fehlerrechnung w, int digits=6){
   double fac = TMath::Power(10.,digits);
   double mittel = TMath::Nint(fac*w.mittel())/fac;
   double emittel = TMath::Nint(fac*w.emittel())/fac;
  cout<<txt<<setw(15)<<setprecision(15)<<mittel<<" +/- "<<emittel<<endl;
}

void cprint(TString txt, fehlerrechnung w, int digits=6){
   double fac = TMath::Power(10.,digits);
   double mittel = TMath::Nint(fac*w.mittel())/fac;
   double emittel = TMath::Nint(fac*w.emittel())/fac;
  cout<<txt<<setw(15)<<setprecision(15)<<mittel<<" +/- "<<emittel<<endl;
}

void cprint(fehlerrechnung w, int digits){
   double fac = TMath::Power(10.,digits);
   double mittel = TMath::Nint(fac*w.mittel())/fac;
   double emittel = TMath::Nint(fac*w.emittel())/fac;
  cout<<setw(12)<<setprecision(15)<<mittel<<" +/-  "<<emittel<<endl;
}

void cprint(fehlerrechnung w){
   double mittel = w.mittel();
   double emittel = w.emittel();
  cout<<setw(12)<<setprecision(15)<<mittel<<" +/- "<<emittel<<endl;
}

// Einlesefunktion
void cread(){
   double mittel, fehler;
  cout<<"Bitte Mittelwert eingeben:";
  cin>>mittel;
  cout<<"Fehler eingeben:";
  cin>>fehler;
  fehlerrechnung z1 = fehlerrechnung(mittel,fehler);
  cout<<" z1 = ";
  cprint(z1);
  cout<<" wurde definiert."<<endl;
}



