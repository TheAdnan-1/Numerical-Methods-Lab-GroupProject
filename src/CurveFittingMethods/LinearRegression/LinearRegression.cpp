#include <bits/stdc++.h>
using namespace std;
#define arr vector<double>
int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    int t;
    fin>>t;
    while(t--){
    int n,i;
    fin>>n;
    arr x(n),y(n);
    for(i=0;i<n;i++) fin>>x[i]>>y[i];
    double sumx=0,sumy=0,sumx2=0,sumxy=0,b,a;
    for(i=0;i<n;i++) {
        sumx += x[i];
        sumy += y[i];
        sumx2 += (x[i]*x[i]);
        sumxy += (x[i]*y[i]);
    }
    b=(n*sumxy-sumx*sumy)/(n*sumx2-sumx*sumx);
    a=(sumy-b*sumx)/n;
    fout<<fixed<<setprecision(3);
    fout<<"Linear Regression Line:\n";
    fout<<"y = a + b x\n";
    fout<<"a = "<<a<<'\n';
    fout<<"b = "<<b<< '\n';
    }
    fin.close();
    fout.close();
    return 0;
}
