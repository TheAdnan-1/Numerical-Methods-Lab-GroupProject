#include <bits/stdc++.h>
using namespace std;

double f(double x,double c[],int deg){
    double res=0;
    for(int i=0;i<=deg;i++) res+=c[i]*pow(x,deg-i);
    return res;
}

int main(){
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int T;
    fin>>T;
    fout<<"Total Test Cases: "<<T<<"\n\n";

    for(int t=1;t<=T;t++){
        int deg;
        fin>>deg;
        double coef[20];
        for(int i=0;i<=deg;i++) fin>>coef[i];

        double a=-10,b=10,tol=0.001;
        int n=50;
        double fa=f(a,coef,deg),fb=f(b,coef,deg);

        fout<<"Function: f(x) = ";
        for(int i=0;i<=deg;i++){
            if(i>0 && coef[i]>=0) fout<<"+";
            fout<<coef[i];
            if(deg-i>0) fout<<"x^"<<(deg-i)<<" ";
        }
        fout<<"\n";
        fout<<"Degree: "<<deg<<"\n";
        fout<<fixed<<setprecision(6);
        fout<<"Error Tolerance: "<<tol<<"\n";
        fout<<"Search Interval: ["<<a<<", "<<b<<"]\n";
        fout<<"Root Approximation:\n";

        if(fa*fb>0){
            fout<<"  Invalid interval (no sign change)\n";
        }else{
            double c;
            for(int i=0;i<n;i++){
                c=(a*fb-b*fa)/(fb-fa);
                double fc=f(c,coef,deg);
                fout<<"  Iter "<<i+1<<": x = "<<c<<"\n";
                if(fabs(fc)<tol) break;
                if(fa*fc<0){
                    b=c;
                    fb=fc;
                }else{
                    a=c;
                    fa=fc;
                }
            }
            fout<<"  Final root ≈ "<<c<<"\n";
        }
    }
    fin.close();
    fout.close();
    return 0;
}
