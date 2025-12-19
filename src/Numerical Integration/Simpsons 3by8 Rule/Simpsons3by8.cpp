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

    int t;
    fin>>t;
    while(t--){
        double a,b,h,sum;
        int n;
        fin>>degree;
        for(int i=0;i<=degree;i++) fin>>coeff[i];
        fin>>a>>b>>n;
        if(n%3!=0){
            fout<<"Invalid n (must be multiple of 3)\n";
            continue;
        }
        h=(b-a)/n;
        sum=f(a)+f(b);
        for(int i=1;i<n;i++){
            double x=a+i*h;
            if(i%3==0) sum+=2*f(x);
            else sum+=3*f(x);
        }
        double result=(3*h/8)*sum;
        fout<<fixed<<setprecision(6);
        fout<<"Integral: "<<result<<'\n';
    }
    fin.close();
    fout.close();
    return 0;
}
