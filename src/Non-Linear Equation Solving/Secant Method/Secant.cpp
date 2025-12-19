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
        double x0,x1,x2,tol;
        int maxIter,iter=0;
        fin>>degree;
        for(int i=0;i<=degree;i++) fin>>coeff[i];
        fin>>x0>>x1>>tol>>maxIter;

        bool found=false;
        while(iter<maxIter){
            if(f(x1)-f(x0)==0) break;
            x2=x1-(f(x1)*(x1-x0))/(f(x1)-f(x0));
            if(fabs(x2-x1)<tol){
                fout<<fixed<<setprecision(6);
                fout<<"Root: "<<x2<<'\n';
                found=true;
                break;
            }
            x0=x1;
            x1=x2;
            iter++;
        }
        if(!found) fout<<"Not Converged\n";
    }
    fin.close();
    fout.close();
    return 0;
}
