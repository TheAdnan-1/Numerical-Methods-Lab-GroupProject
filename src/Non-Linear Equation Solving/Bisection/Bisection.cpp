#include <bits/stdc++.h>
using namespace std;

int degree;
double coeff[20];

double f(double x){
    double result=0;
    for(int i=0;i<=degree;i++) result+=coeff[i]*pow(x,degree-i);
    return result;
}

int main(){
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    double a,b,mid;
    int n;
    fin>>degree;
    for(int i=0;i<=degree;i++) fin>>coeff[i];
    fin>>a>>b>>n;
    if(f(a)*f(b)>=0){
        fout<<"Invalid interval\n";
        return 0;
    }
    for(int i=0;i<n;i++){
        mid=(a+b)/2;
        if(f(mid)==0) break;
        else if(f(a)*f(mid)<0) b=mid;
        else a=mid;
    }
    fout<<fixed<<setprecision(6)<<mid<<'\n';

    fin.close();
    fout.close();
    return 0;
}
