#include <bits/stdc++.h>
using namespace std;

double f(const vector<double>& coeffs,double x){
    double res=0,xn=1;
    for(double c:coeffs){
        res+=c*xn;
        xn*=x;
    }
    return res;
}

double df(const vector<double>& coeffs,double x){
    double res=0,xn=1;
    for(int i=1;i<coeffs.size();i++){
        res+=i*coeffs[i]*xn;
        xn*=x;
    }
    return res;
}

double newtonRaphson(const vector<double>& coeffs,double x0,double tol=1e-6,int maxIter=1000){
    for(int i=0;i<maxIter;i++){
        double dfx=df(coeffs,x0);
        if(fabs(dfx)<1e-12) break;
        double x1=x0-f(coeffs,x0)/dfx;
        if(fabs(x1-x0)<tol) return x1;
        x0=x1;
    }
    return x0;
}

int main(){
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int t;
    fin>>t;
    fout<<fixed<<setprecision(6);

    while(t--){
        int degree;
        fin>>degree;
        vector<double> coeffs(degree+1);
        for(int i=0;i<=degree;i++) fin>>coeffs[i];
        double xmin,xmax;
        fin>>xmin>>xmax;

        double step=0.5;
        vector<double> roots;

        for(double x=xmin;x<xmax;x+=step){
            if(f(coeffs,x)*f(coeffs,x+step)<=0){
                double root=newtonRaphson(coeffs,x);
                if(none_of(roots.begin(),roots.end(),
                    [&](double r){return fabs(r-root)<1e-4;})){
                    roots.push_back(root);
                }
            }
        }

        if(roots.empty()){
            fout<<"No real roots found\n";
        }else{
            for(double r:roots) fout<<"Root: "<<r<<"\n";
        }
        fout<<"----\n";
    }
    fin.close();
    fout.close();
    return 0;
}
