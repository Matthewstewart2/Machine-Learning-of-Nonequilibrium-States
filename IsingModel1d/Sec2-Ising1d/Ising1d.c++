/*
Ising1d with periodic boundaries
prints on screen
time, magnetisation/spin
outputs one file ev 10 MC sweeps with spin values
 */

#include<iostream>
#include<stdlib.h>
#include<cmath>
#include<sstream>
#include<fstream>
#include<string>
#include<iomanip>
#include<vector>
using namespace std;

/* CONSTANTS */
const int N=100; //size
const int Tmax=10000;

/* ARRAYS */
int spin[N]; //spin
int spinnew;

/* PARAMETERS */
double j=1.0; //coupling constant

/**/
double energy;
double newenergy;

/* MAIN */
int main(int argc, char* argv[]){
srand(time(NULL));

double temp=atof(argv[1]);

/*Start loop over time*/
for(int t=0; t<=Tmax; t++){

rand();

if(t==0){
//DEFINE lattice at random up/down
for(int n=0;n<N; n++) spin[n]=(rand()%2*2-1); //(-1,1)

//compute initial energy
energy=0;
for(int n=0;n<N; n++){
    int left=n-1;
    if(left==-1)left=N-1;
    //summing only s_left*s_n to avoid double-counting
    energy+=-j*(spin[left]*spin[n]);
}
cout <<"inital energy per spin " << energy*1.0/N <<endl;
//    for(int n=0;n<N; n++){
//    write << n << " " << spin[n]<< endl;
//}
}
    
/*Monte Carlo sweep*/
for(int n=0;n<N; n++){

//pick a site at random
int i=int(rand()*1.0/RAND_MAX*N);
//w periodic boundaries
int left=i-1; if(left==-1)left+=N;
int right=i+1; if(right==N)right-=N;
    
//ATTEMPT a spin flip
spinnew=-spin[i];

//compute delta energy (new - old)
double deltaenergy=-j*(spinnew*spin[left]+spinnew*spin[right]-spin[i]*spin[left]-spin[i]*spin[right]);
    
    //cout << "Flipped " << i << " from " << spin[left] << spin[i] << spin[right] << " to " << spin[left] << spinnew << spin[right] << " -> " << deltaenergy <<endl;
    
//accept move w Metropolis criterion
//i.e. w probability =exp(-DeltaEnergy/kBT)
double randf=rand()*1.0/RAND_MAX;
    if(randf<exp(-deltaenergy/temp)){
        spin[i]=spinnew;
        energy+=deltaenergy;
        //cout << "accepted!" << randf << " < " << exp(-deltaenergy/temp)<<endl;
    }
 
   // cin.get();
}


if(t%10==0){
stringstream file;
file<<"out."<<t<<".dat";
ofstream write;
write.open(file.str().c_str());

double magnetisation=0;
for(int n=0;n<N; n++){
        write << n << " " << spin[n]<< endl;
        magnetisation+=spin[n];
}
cout << t << " " << magnetisation*1.0/N <<endl;
}
    
}//close loop over time
    

return 0;
}

