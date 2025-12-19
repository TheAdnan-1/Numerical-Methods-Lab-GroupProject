#include <bits/stdc++.h>
using namespace std;
#define arr vector<double>
int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    int t;
    fin>>t;
    while(t--){
    int model,n,i;
    fin>>model>>n;
    arr x(n),y(n);
    for(i=0;i<n;i++)fin>>x[i]>>y[i];
    fout<<fixed<<setprecision(3);
    if(model==1){
        // y=ae^(bx)
        arr Y(n);
        for(i=0;i<n;i++) Y[i] = log(y[i]);
        double sx=0, sy=0, sxx=0, sxy=0,b,a;
        for(i=0;i<n;i++){
            sx += x[i];
            sy += Y[i];
            sxx += (x[i]*x[i]);
            sxy += (x[i]*Y[i]);
        }
        b=(n*sxy-sx*sy)/(n*sxx-sx*sx);
        a=exp((sy-b*sx)/n);
        fout<<"a = "<<a<<'\n';
        fout<<"b = "<<b<<'\n';
        fout<<"Equation: y = "<<a<<" e^("<<b<<" x)\n";
    }
    else if(model==2){
        // y=ax^b
        arr X(n),Y(n);
        for(i=0;i<n;i++){
            X[i]=log(x[i]);
            Y[i]=log(y[i]);
        }
        double sx=0, sy=0, sxx=0, sxy=0,b,a;
        for(i=0;i<n;i++){
            sx += X[i];
            sy += Y[i];
            sxx += X[i]*X[i];
            sxy += X[i]*Y[i];
        }
        b=(n*sxy-sx*sy)/(n*sxx-sx*sx);
        a=exp((sy-b*sx)/n);
        fout<<"a = "<<a<<'\n';
        fout<<"b = "<<b<<'\n';
        fout<<"Equation: y = "<<a<<" x^"<<b<<'\n';
    }
    else if(model == 3){
        // y = a + b e^(x/z)
        double z;
        fin>>z;
        arr X(n);
        for(i=0;i<n;i++) X[i]=exp(x[i]/z);
        double sx=0, sy=0, sxx=0, sxy=0,a,b;
        for(i=0;i<n;i++){
            sx += X[i];
            sy += y[i];
            sxx += X[i]*X[i];
            sxy += X[i]*y[i];
        }
        b=(n*sxy-sx*sy)/(n*sxx-sx*sx);
        a=(sy-b*sx)/n;
        fout<<"a = "<<a<<'\n';
        fout<<"b = "<<b<<'\n';
        fout<<"Equation: y = "<<a<<" + "<<b<<" e^(x/"<<z<<")\n";
    }
    else fout<<"Invalid Option\n";
    }
    fin.close();
    fout.close();
    return 0;
}
