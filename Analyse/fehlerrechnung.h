#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <fstream>
//#include <fehlerrechnungtools.h>
using namespace std;


class fehlerrechnung{
private: 
   double m;
   double em;
   
public:

 double &mittel(){return m;} 
 double &emittel(){return em;}
 double signifikanz(){return m/em;}

//########## konstruktoren
fehlerrechnung(  double mittelwert=0,  double fehler=0){m=mittelwert;em=fehler;}
fehlerrechnung( double mittelwert, std::string name):m(0),em(0){
  if(name == std::string("poisson")){
    m = mittelwert;
    em = sqrt(mittelwert);
    //cout << "Der Konstruktor setzt " <<  m << " " << em;
  }
  else cout << "Falscher Aufruf des Poisson-Konstruktors. Aufruf mit fehlerrechnung( mittelwert , \"poisson\") " << endl;
}
//##########

//########## Neudefinitionen Operatoren
friend fehlerrechnung operator+(fehlerrechnung w1, fehlerrechnung w2){
 double m1 = w1.m;
 double em1= w1.em;
 double m2 = w2.m;
 double em2= w2.em;	
	return fehlerrechnung(m1 + m2, sqrt( em1*em1 + em2*em2 ) );
}
friend fehlerrechnung operator-(fehlerrechnung w1, fehlerrechnung w2){
 double m1 = w1.m;
 double em1= w1.em;
 double m2 = w2.m;
 double em2= w2.em;	
	return fehlerrechnung(m1 - m2, sqrt( em1*em1 + em2*em2 ) );
}

friend fehlerrechnung operator*(fehlerrechnung w1, fehlerrechnung w2){
 double m1 = w1.m;
 double em1= w1.em;
 double m2 = w2.m;
 double em2= w2.em;	
	return fehlerrechnung(m1*m2, sqrt( m2*m2*em1*em1 + m1*m1*em2*em2 ) );
}
friend fehlerrechnung operator*(double w1, fehlerrechnung w2){
 double m2 = w2.m;
 double em2= w2.em;	
	return fehlerrechnung(w1*m2, w1*em2 );
}
friend fehlerrechnung operator/(fehlerrechnung w1, fehlerrechnung w2){
 double m1 = w1.m;
 double em1= w1.em;
 double m2 = w2.m;
 double em2= w2.em;	
	if (m2 != 0) return fehlerrechnung( m1/m2,sqrt( (em1*em1)/(m2*m2) + m1*m1*em2*em2/(m2*m2*m2*m2) ));
	else {
	cout<<"der zweite Mittelwert ist null, dadurch kann ich nicht teilen, du Nase"<<endl;
	return 666;
	}
}
friend fehlerrechnung operator/(fehlerrechnung w1, double w2){
 double m1 = w1.m;
 double em1= w1.em;
	if (w2 != 0) return fehlerrechnung( m1/w2 , em1/w2  );
	else {
	cout<<"Du versuchst durch Null zu teilen. Vergeblich!"<<endl;
	return 666;
	}
}

friend bool operator==(fehlerrechnung w1, fehlerrechnung w2){
return (w1.m==w2.m && w1.em==w2.em);
}
friend bool operator!=(fehlerrechnung w1, fehlerrechnung w2){
return (w1.m!=w2.m || w1.em==w2.em);
}
//################

//########## spezielle Funktionen
// Exponential
friend fehlerrechnung fexp(fehlerrechnung w1){
return fehlerrechnung(exp(w1.m), w1.em*exp(w1.m));
}
// Wurzel
friend fehlerrechnung fsqrt(fehlerrechnung w1){
return fehlerrechnung(sqrt(w1.m), 0.5/sqrt(w1.m)*w1.em);
}
// Potenzen
template <class T>
friend fehlerrechnung fpow(fehlerrechnung w1, T expo){
return fehlerrechnung(pow(w1.m,expo), expo*pow(w1.m,expo-1)*w1.em);
}
// natuerlicher Logarithmus
friend fehlerrechnung flog(fehlerrechnung w1) {
return fehlerrechnung(log( w1.m ), w1.em/w1.m );
}
// sin
friend fehlerrechnung fsin(fehlerrechnung w1) {
return fehlerrechnung(sin( w1.m ), cos(w1.m)*w1.em );
}
// cos
friend fehlerrechnung fcos(fehlerrechnung w1) {
return fehlerrechnung(cos( w1.m ), sin(w1.m)*w1.em );
}
// sinh
friend fehlerrechnung fsinh(fehlerrechnung w1) {
return fehlerrechnung(sinh( w1.m ), cosh(w1.m)*w1.em );
}
// cosh
friend fehlerrechnung fcosh(fehlerrechnung w1) {
return fehlerrechnung(cosh( w1.m ), sinh(w1.m)*w1.em );
}


};



