#include "TROOT.h"
#include "TCanvas.h"
#include "fehlerrechnungtools.h"
#include "TTree.h"
#include "TFile.h"


void WQVergleich(){
    
   using namespace std;

   cout << "sigma(10 TeV)/sigma(7 TeV):" << endl;
   cout << "Ttbar            :" << 385.      / 165. << endl;
   cout << "WW               :" <<  74.      / 43. << endl;
   cout <<  "ZZ               :" << 10.5      / 5.9  << endl;
   cout << "WZ               :" << 32.       / 18. << endl;
   cout << "Zee              :" << 1944.     / 1300. << endl;
   cout << "Zmumu            :" << 1944.     / 1300. << endl;
   cout << "Ztautau          :" << 1944.     / 1300. << endl;
   cout << "InclusiveMu15    :" << 0.5091e9 / 0.2969e9 << endl;
   cout << "BCtoE20to30      :" << 0.4e9     / 0.2355e9 << endl;
   cout << "BCtoE30to80      :" << 0.1e9     / 0.0593e9 << endl;
   cout << "BCtoE80to170     :" << 1.9e6     / 0.906e6 << endl;
   cout << "EMenriched20to30 :" << 0.4e9     / 0.2355e9 << endl;
   cout << "EMenriched30to80 :" << 0.1e9     / 0.0593e9 << endl;
   cout << "EMenriched80to170:" << 1.9e6     / 0.906e6 << endl;
   cout << "Wenu             :" << 42800./3  / (28000./3) << endl;
   cout << "Wmunu            :" << 42800./3  / (28000./3) << endl;
   cout << "SingleTops       :" << 5.        / 4.6 << endl;
   cout << "SingleTopt       :" << 130.      / 63. << endl;
   cout << "SingleToptW      :" << 29.       / 10.6 << endl;
   cout << endl;
   cout << "genauere Betrachtung von W+Jets: Cut-Flow des Verhaeltnisses:" << endl;
   cout << "Der Cut-Flow ist: Alle -- nach Trigger -- nach 2 Leptonen -- nach 2 Jets -- nach versch. Ladung -- nach SameFlavour-Schnitten" << endl;
   cout << "Wenu:" << endl;
   cout <<105288 /  72706.667 << endl;
   cout << 85943/    59840.487 << endl;
   cout << 8.293/    4.863<< endl;
   cout << 0.74/    0.35 << endl;
   cout << 0.494/    0.07 << endl;
   cout << 0.247/    0.07 << endl;
     cout << "Wmunu:" <<  endl;
   cout << 98582.667/  69253.333<< endl;
   cout << 71604.645/    50835.658 << endl;
   cout << 8.473/    6.051 << endl;
   cout << 0.568/    0.372 << endl;
   cout << 0.521/    0.203 << endl;
   cout << 0.521/    0.203 << endl;

}
