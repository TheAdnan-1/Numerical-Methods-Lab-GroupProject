#include <bits/stdc++.h>
using namespace std;
#define arr vector<double>
int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    int t;
    fin>>t;
    while(t--){
    int n, m,i,row,col,k,j;
    double sum;
    fin>>n>>m;
    arr x(n), y(n);
    for(i=0;i<n;i++) fin>>x[i]>>y[i];
    vector<arr>A(m+1,arr(m+2));
    for(row=0;row<m+1;row++){
        for(col=0;col<m+1;col++){
            sum=0;
            for(i=0;i<n;i++) sum += pow(x[i],row+col);
            A[row][col]=sum;
        }
        sum=0;
        for(i=0;i<n;i++)
            sum += y[i]*pow(x[i], row);
        A[row][m+1] = sum;
    }
    for(i=0;i<m+1;i++){
        for(k=i+1;k<m+1;k++){
            double factor = A[k][i]/A[i][i];
            for(j=i; j <= m+1; j++)
                A[k][j] -= factor * A[i][j];
        }
    }
    arr a(m+1);
    for(i=m;i>=0;i--) {
        a[i] = A[i][m+1];
        for(j=i+1;j<m+1;j++) a[i] -= A[i][j] * a[j];
        a[i] /= A[i][i];
    }
    fout<<fixed<<setprecision(3);
    fout<<"Polynomial Regression (degree "<<m<<"):\n";
    for(i=0;i<=m;i++)fout<<'a'<<i<<" = "<<a[i]<<'\n';
    }
    fin.close();
    fout.close();
    return 0;
}
