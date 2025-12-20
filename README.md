# NUMERICAL METHODS LABORATORY GROUP PROJECT

A Numerical Methods group project built with C++.  

## Project Introduction
This repository contains implementations of some fundamental numerical methods written in C++. Every method supports multiple test cases and allows input to be read from files with results displayed in both output files and the console. The methods are logically grouped based on their problem domains.



## Project Structure

```
README.md
├── Solution of Non-Linear Equations
│   ├── Bisection
│   ├── False Position
│   ├── Newton Raphson
│   └── Secant
│
├── Solution of Linear Equations
│   ├── Gauss Elimination
│   ├── Gauss Jordan Elimination
│   ├── LU Decomposition
│   └── Matrix Inversion
│
├── Differential Equation Solving
│   └── Runge Kutta 4th Order
│
├── Interpolation Methods
│   ├── Newton's Forward Interpolation
│   ├── Newton's Backward Interpolation
│   └── Newton's Divided Difference Interpolation
│
├── Numerical Differentiation
│   ├── Differentiation by Forward Interpolation
│   └── Differentiation by Backward Interpolation
│
├── Curve Fitting and Regression
│   ├── Linear Regression
│   ├── Polynomial Regression
│   └── Transcendental Regression
│
└── Numerical Integration
    ├── Simpson's 1/3 Rule
    └── Simpson's 3/8 Rule
```

Each method folder contains:
- `Method.txt` - Briefly explained theory and features
- `Method.cpp` - C++ implementation of the method
- `input.txt` - Sample test case inputs
- `output.txt` - Expected outputs and results

---

## Table of Contents

- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
  - [False Position Method](#false-position-method)
  - [Newton Raphson Method](#newton-raphson-method)
  - [Secant Method](#secant-method)
- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
  - [LU Decomposition Method](#lu-decomposition-method)
  - [Matrix Inversion Method](#matrix-inversion-method)
- [Differential Equation Solving](#differential-equation-solving)
  - [Runge-Kutta 4th Order Method](#runge-kutta-4th-order-method)
- [Interpolation Methods](#interpolation-methods)
  - [Newton Forward Interpolation Method](#newton-forward-interpolation-method)
  - [Newton Backward Interpolation Method](#newton-backward-interpolation-method)
  - [Newton Divided Difference Interpolation Method](#newton-divided-difference-interpolation-method)
- [Numerical Differentiation](#numerical-differentiation)
  - [Differentiation by Forward Interpolation Method](#differentiation-by-forward-interpolation-method)
  - [Differentiation by Backward Interpolation Method](#differentiation-by-backward-interpolation-method)
- [Curve Fitting & Regression](#curve-fitting--regression)
  - [Linear Regression Method](#linear-regression-method)
  - [Polynomial Regression Method](#polynomial-regression-method)
  - [Transcendental Regression Method](#transcendental-regression-method)
- [Numerical Integration](#numerical-integration)
  - [Simpson's 1/3 Rule Method](#simpsons-13-rule-method)
  - [Simpson's 3/8 Rule Method](#simpsons-38-rule-method)
- [Purpose of This Project](#purpose-of-this-project)
- [Authors](#authors)

---

## Solution of Non-Linear Equations

### Bisection Method

**Theory**
- The Bisection Method is a numerical technique used to find a real root of a nonlinear equation f(x) = 0.
- It requires an initial interval [a, b] such that the function is continuous and f(a)·f(b) < 0, ensuring the existence of a root.
- The interval is repeatedly divided into two equal halves by computing the midpoint c = (a + b) / 2.
- The sign of f(c) determines the new interval: if f(a)·f(c) < 0, the root lies in [a, c]; otherwise, it lies in [c, b].
- This process continues until the interval length becomes sufficiently small or the function value at the midpoint is close to zero.
- The method always converges for a valid initial interval but has linear (slow) convergence.
- The maximum error after n iterations is (b − a) / 2ⁿ.


**Code**
```cpp
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
```
**Input**

```bash
3
2
1 0 -4
3
1 0 -1 -2
2
1 -3 2
```

**Output**

```bash
2
1.52197
1
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### False Position Method

**Theory**

-The False Position method is a numerical technique used to find a real root of a nonlinear equation f(x) = 0.
-It requires an initial interval [a, b] where the function is continuous and f(a)·f(b) < 0.
-The root is approximated by a straight-line interpolation between the points (a, f(a)) and (b, f(b)).
-The interval is updated by replacing the endpoint having the same sign as f(c), and the process is repeated until convergence.

**Code**

```cpp
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
```

**Input**

```bash
3
4
4 -9 0 0 4
3
1 -6 11 -6
2
1 0 -4
```

**Output**

```bash
Total Test Cases: 3

Function: f(x) = 4x^4 -9x^3 +0x^2 +0x^1 +4
Degree: 4
Error Tolerance: 0.001000
Search Interval: [-10.000000, 10.000000]
Root Approximation:
  Iter 1: x = 0.444444
  Iter 2: x = 0.694444
  Iter 3: x = 0.826087
  Iter 4: x = 0.887097
  Iter 5: x = 0.915254
  ...
  Final root ≈ 0.906250

Function: f(x) = 1x^3 -6x^2 +11x^1 -6
Degree: 3
Error Tolerance: 0.001000
Search Interval: [-10.000000, 10.000000]
Root Approximation:
  Iter 1: x = 1.000000
  Iter 2: x = 2.000000
  Iter 3: x = 3.000000
  ...
  Final root ≈ 0.999554

Function: f(x) = 1x^2 +0x^1 -4
Degree: 2
Error Tolerance: 0.001000
Search Interval: [-10.000000, 10.000000]
Root Approximation:
  Iter 1: x = 0.400000
  Iter 2: x = 1.333333
  Iter 3: x = 1.600000
  Iter 4: x = 1.777778
  Iter 5: x = 1.846154
  ...
  Final root ≈ 2.000186
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Newton Raphson Method

**Theory**

-The Newton–Raphson method is a numerical technique used to find a root of a nonlinear equation f(x) = 0.
-It starts with an initial approximation and uses the derivative of the function to improve accuracy.
-The iterative formula is xₙ₊₁ = xₙ − f(xₙ) / f′(xₙ).
-The method converges rapidly (quadratic convergence) when the initial guess is close to the root.

**Code**

```cpp
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
```

**Input**

```bash
2
3
-6 11 -6 1
0 4
2
-4 0 1
-3 3

```

**Output**

```bash
Root: 1.000000
Root: 2.000000
Root: 3.000000
----
Root: -2.000000
Root: 2.000000
----
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Secant Method

**Theory**

-The Secant method is a numerical technique used to find a root of a nonlinear equation f(x) = 0.
-It starts with two initial approximations and does not require the derivative of the function.
-The iterative formula uses a secant line through the points (xₙ₋₁, f(xₙ₋₁)) and (xₙ, f(xₙ)).
-The method converges faster than the bisection method but slower than Newton–Raphson.

**Code**

```cpp
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
```

**Input**

```bash
2
3
1 0 -1 -2
1 2
0.0001
100
2
1 -3 2
0 3
0.0001
100
```

**Output**

```bash
Root: 1.52138
Root: 1
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

## Solution of Linear Equations

### Gauss Elimination Method

**Theory**
- Gaussian elimination is a method to solve systems of linear equations.  
- It works by converting the augmented matrix [A|B] into an upper-triangular form using row operations.  
- Partial pivoting improves stability: rows are swapped so the largest pivot in each column is placed on the diagonal. This avoids dividing by very small numbers and reduces numerical errors.  
- If a pivot column becomes all zeros but the augmented part is non-zero, the system has no solution.  
- If both are zero, the system may have infinitely many solutions.  
- Once the matrix is triangular with non-zero diagonal entries, the system Ux = C is solved using back substitution, starting from the last row upward.  


**Algorithm Steps**

1. Apply partial pivoting in each column to bring the largest pivot to the diagonal.  
2. Perform forward elimination to form an upper-triangular matrix.  
3. Check for consistency - zero pivots or zero rows indicate special cases.  
4. If the system is consistent and has a unique solution, solve using back substitution.  


**Notes**

- A tolerance of `1e-10` is used to treat very small pivots or coefficients as zero.  
- This prevents instability caused by floating-point precision limits.  


**Code**

```cpp
#include <bits/stdc++.h>
using namespace std;

const double eps = 1e-9;

void printMatrix(ofstream &fout, vector<vector<double>> &mat, int n, int m) {
    fout << "Matrix:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            fout << setw(10) << fixed << setprecision(3) << mat[i][j] << " ";
        }
        fout << endl;
    }
    fout << endl;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    fout << "GAUSSIAN ELIMINATION\n\n";

    int testCases;
    fin >> testCases;
    fout << "Number of test cases: " << testCases << "\n\n";

    for (int tc = 1; tc <= testCases; tc++) {
        int n;
        fin >> n;

        vector<vector<double>> mat(n, vector<double>(n + 1));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= n; j++) {
                fin >> mat[i][j];
            }
        }

        fout << "\nTest Case " << tc << " (" << n << "x" << n << ")\n";
        printMatrix(fout, mat, n, n + 1);

        int row = 0;
        for (int col = 0; col < n && row < n; col++) {
            int pivot = row;
            for (int i = row + 1; i < n; i++) {
                if (fabs(mat[i][col]) > fabs(mat[pivot][col]))
                    pivot = i;
            }

            if (fabs(mat[pivot][col]) < eps)
                continue;

            swap(mat[row], mat[pivot]);

            for (int i = row + 1; i < n; i++) {
                double factor = mat[i][col] / mat[row][col];
                for (int j = col; j <= n; j++) {
                    mat[i][j] -= factor * mat[row][j];
                }
            }
            row++;
        }

        int rankA = 0, rankAug = 0;

        for (int i = 0; i < n; i++) {
            bool nonZeroA = false;
            bool nonZeroAug = false;

            for (int j = 0; j < n; j++) {
                if (fabs(mat[i][j]) > eps)
                    nonZeroA = true;
            }
            if (fabs(mat[i][n]) > eps)
                nonZeroAug = true;

            if (nonZeroA)
                rankA++;
            if (nonZeroA || nonZeroAug)
                rankAug++;
        }

        if (rankA < rankAug) {
            fout << "No Solution!\n";
        }
        else if (rankA < n) {
            fout << "Infinite Solutions!\n";
        }
        else {
            vector<double> solution(n);

            for (int i = n - 1; i >= 0; i--) {
                solution[i] = mat[i][n];
                for (int j = i + 1; j < n; j++) {
                    solution[i] -= mat[i][j] * solution[j];
                }
                solution[i] /= mat[i][i];
            }

            fout << "Unique Solution!\n";
            for (int i = 0; i < n; i++) {
                fout << "x" << i + 1 << " = " 
                     << fixed << setprecision(3) 
                     << solution[i] << endl;
            }
        }
    }

    fin.close();
    fout.close();
    return 0;
}
```

**Input**

```bash
3
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
2
1 2 3
2 4 5
3
2 4 6 8
1 2 3 4
3 6 9 12
```

**Output**

```bash
GAUSSIAN ELIMINATION

Test Case 1 (3x3)
Matrix:
     2.000      1.000     -1.000      8.000
    -3.000     -1.000      2.000    -11.000
    -2.000      1.000      2.000     -3.000

Unique Solution!
x1 = 2.000
x2 = 3.000
x3 = -1.000

Test Case 2 (2x2)
Matrix:
     1.000      2.000      3.000
     2.000      4.000      5.000

No Solution!

Test Case 3 (3x3)
Matrix:
     2.000      4.000      6.000      8.000
     1.000      2.000      3.000      4.000
     3.000      6.000      9.000     12.000

Infinite Solutions!
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Gauss Jordan Elimination Method

**Theory**

- Gauss-Jordan applies row operations to convert [A|B] directly to RREF, where each pivot is 1 and the pivot columns have zeros elsewhere. 

- Partial pivoting improves numerical stability and avoids dividing by tiny pivots.

- Rank comparison:  if rank(A) < rank([A|B]) → inconsistent; if rank(A) < n but rank(A) == rank([A|B]) → infinite solutions; if rank(A) = n → unique solution.

### Algorithm steps

1. **Partial pivoting**: swap the current row with the row having the largest absolute pivot in the column. 

2. **Scale the pivot row** so the pivot becomes 1.

3. **Eliminate the pivot column** in all other rows to reach RREF.

4. **Check ranks** to classify:  inconsistent, infinite solutions, or unique solution (read directly from RREF).

### What it does: 

- Solves linear systems using Gauss-Jordan elimination with partial pivoting. 

- Reduces the augmented matrix to Reduced Row Echelon Form (RREF).

- Detects three outcomes per test case:  Unique Solution, No Solution (inconsistent), Infinite Solutions (dependent).

- Handles multiple test cases in one run, reading from Input. txt and writing to Output.txt while printing to console. 

**Code**

```cpp
#include <bits/stdc++.h>
using namespace std;
const double eps = 1e-9;

void printMatrix(ofstream &fout, vector<vector<double>> &mat, int n, int m) {
    fout << "Matrix:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            fout << setw(10) << fixed << setprecision(3) << mat[i][j] << " ";
        }
        fout << endl;
    }
    fout << endl;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    fout << "GAUSS JORDAN ELIMINATION\n\n";

    int testCases;
    fin >> testCases;
    fout << "Number of test cases: " << testCases << "\n\n";

    for (int tc = 1; tc <= testCases; tc++) {

        int n;
        fin >> n;

        vector<vector<double>> mat(n, vector<double>(n + 1));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= n; j++) {
                fin >> mat[i][j];
            }
        }

        fout << "\nTest Case " << tc << " (" << n << "x" << n << ")\n";
        printMatrix(fout, mat, n, n + 1);

        int row = 0;
        for (int col = 0; col < n && row < n; col++) {

            int pivot = row;
            for (int i = row; i < n; i++) {
                if (fabs(mat[i][col]) > fabs(mat[pivot][col]))
                    pivot = i;
            }

            if (fabs(mat[pivot][col]) < eps)
                continue;

            swap(mat[row], mat[pivot]);

            double pivotValue = mat[row][col];
            for (int j = 0; j <= n; j++) {
                mat[row][j] /= pivotValue;
            }

            for (int i = 0; i < n; i++) {
                if (i == row) continue;
                double factor = mat[i][col];
                for (int j = 0; j <= n; j++) {
                    mat[i][j] -= factor * mat[row][j];
                }
            }

            row++;
        }

        int rankA = 0, rankAug = 0;

        for (int i = 0; i < n; i++) {
            bool nonZeroA = false;
            bool nonZeroAug = false;

            for (int j = 0; j < n; j++) {
                if (fabs(mat[i][j]) > eps)
                    nonZeroA = true;
            }
            if (fabs(mat[i][n]) > eps)
                nonZeroAug = true;

            if (nonZeroA)
                rankA++;
            if (nonZeroA || nonZeroAug)
                rankAug++;
        }

        if (rankA < rankAug) {
            fout << "No Solution!\n";
        }
        else if (rankA < n) {
            fout << "Infinite Solutions!\n";
        }
        else {
            fout << "Unique Solution!\n";
            for (int i = 0; i < n; i++) {
                fout << "x" << i + 1 << " = "
                     << fixed << setprecision(3)
                     << mat[i][n] << endl;
            }
        }
    }

    fin.close();
    fout.close();
    return 0;
}
```

**Input**

```bash
3
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
2
1 2 3
2 4 5
3
2 4 6 8
1 2 3 4
3 6 9 12
```

**Output**

```bash
GAUSS JORDAN ELIMINATION

Test Case 1 (3x3)
Matrix:
     2.000      1.000     -1.000      8.000
    -3.000     -1.000      2.000    -11.000
    -2.000      1.000      2.000     -3.000

Unique Solution!
x1 = 2.000
x2 = 3.000
x3 = -1.000

Test Case 2 (2x2)
Matrix:
     1.000      2.000      3.000
     2.000      4.000      5.000

No Solution!

Test Case 3 (3x3)
Matrix:
     2.000      4.000      6.000      8.000
     1.000      2.000      3.000      4.000
     3.000      6.000      9.000     12.000

Infinite Solutions!
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### LU Decomposition Method

**Theory**
- LU decomposition breaks down matrix A into two triangular matrices: L (lower) and U (upper), where A = LU and L has unit diagonal elements.

- Computational complexity: decomposition requires O(n^3) operations; solving triangular systems takes O(n^2) steps.

- Singular matrix detection: when any diagonal element of U equals zero, the determinant is zero.

### Algorithm steps:

1. **Apply Doolittle factorization**: construct L with ones on diagonal and U with computed pivot values from matrix A.

2. **Forward pass**: solve Ly = b using forward substitution to find intermediate vector y.

3. **Singularity check**: examine diagonal entries of U; if |U[i][i]| falls below epsilon threshold, mark as singular.

4. **Classify singular cases**:  zero pivot with non-zero corresponding y-value indicates inconsistency; zero pivot with zero y-value suggests dependent equations.

5. **Backward pass**: when matrix is non-singular, apply back substitution to solve Ux = y and output the solution; display L and U matrices in all cases.

### What it does:

- Decomposes coefficient matrix into lower and upper triangular factors using Doolittle's algorithm for solving linear equations.

- Processes multiple systems from Input.txt file, generates solutions in Output.txt, and displays results on screen.

- Outputs both L and U decomposition matrices; determines system status:  single solution, inconsistent system, or dependent system with multiple solutions.


**Code**

```cpp
#include <bits/stdc++.h>
using namespace std;
const double eps = 1e-9;

void printMatrix(ofstream &fout, vector<vector<double>> &a, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            fout << setw(10) << fixed << setprecision(3) << a[i][j] << " ";
        }
        fout << endl;
    }
    fout << endl;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    fout << "LU DECOMPOSITION (Doolittle)\n\n";

    int testCases;
    fin >> testCases;

    for (int tc = 1; tc <= testCases; tc++) {
        int n;
        fin >> n;

        vector<vector<double>> A(n, vector<double>(n));
        vector<double> B(n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                fin >> A[i][j];
            fin >> B[i];
        }

        fout << "Test Case " << tc << " (" << n << "x" << n << ")\n";
        fout << "Matrix:\n";

        vector<vector<double>> aug(n, vector<double>(n + 1));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                aug[i][j] = A[i][j];
            aug[i][n] = B[i];
        }
        printMatrix(fout, aug, n, n + 1);

        vector<vector<double>> L(n, vector<double>(n, 0.0));
        vector<vector<double>> U(n, vector<double>(n, 0.0));

        bool singular = false;

        for (int i = 0; i < n; i++) {
            L[i][i] = 1.0;

            for (int j = i; j < n; j++) {
                double sum = 0;
                for (int k = 0; k < i; k++)
                    sum += L[i][k] * U[k][j];
                U[i][j] = A[i][j] - sum;
            }

            if (fabs(U[i][i]) < eps) {
                singular = true;
                break;
            }

            for (int j = i + 1; j < n; j++) {
                double sum = 0;
                for (int k = 0; k < i; k++)
                    sum += L[j][k] * U[k][i];
                L[j][i] = (A[j][i] - sum) / U[i][i];
            }
        }

        fout << "L matrix:\n";
        printMatrix(fout, L, n, n);

        fout << "U matrix:\n";
        printMatrix(fout, U, n, n);

        if (singular) {
            int rankA = 0, rankAug = 0;

            for (int i = 0; i < n; i++) {
                bool nzA = false, nzAug = false;

                for (int j = 0; j < n; j++)
                    if (fabs(aug[i][j]) > eps) nzA = true;
                if (fabs(aug[i][n]) > eps) nzAug = true;

                if (nzA) rankA++;
                if (nzA || nzAug) rankAug++;
            }

            fout << "Matrix A is singular (det A = 0)\n";
            if (rankA < rankAug)
                fout << "-> System has \"NO\" solution.\n\n";
            else
                fout << "-> System has infinite solutions.\n\n";
            continue;
        }

        vector<double> y(n);
        for (int i = 0; i < n; i++) {
            y[i] = B[i];
            for (int j = 0; j < i; j++)
                y[i] -= L[i][j] * y[j];
        }

        vector<double> x(n);
        for (int i = n - 1; i >= 0; i--) {
            x[i] = y[i];
            for (int j = i + 1; j < n; j++)
                x[i] -= U[i][j] * x[j];
            x[i] /= U[i][i];
        }

        fout << "UNIQUE Solution!\n";
        for (int i = 0; i < n; i++) {
            fout << "x" << i + 1 << " = "
                 << fixed << setprecision(3)
                 << x[i] << endl;
        }
        fout << endl;
    }

    fin.close();
    fout.close();
    return 0;
}
```

**Input**

```bash
4
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
2
1 2 3
2 4 5
3
2 4 6 8
1 2 3 4
3 6 9 12
2
1 1 2
1 1 3
```

**Output**

```bash
LU DECOMPOSITION (Doolittle)

Test Case 1 (3x3)
Matrix:
     2.000      1.000     -1.000      8.000
    -3.000     -1.000      2.000    -11.000
    -2.000      1.000      2.000     -3.000

L matrix:
     1.000      0.000      0.000
    -1.500      1.000      0.000
    -1.000      4.000      1.000

U matrix:
     2.000      1.000     -1.000
     0.000      0.500      0.500
     0.000      0.000     -1.000

UNIQUE Solution!
x1 = 2.000
x2 = 3.000
x3 = -1.000


Test Case 2 (2x2)
Matrix:
     1.000      2.000      3.000
     2.000      4.000      5.000

L matrix:
     1.000      0.000
     2.000      1.000

U matrix:
     1.000      2.000
     0.000      0.000

Matrix A is singular (det A = 0)
-> System has "NO" solution.


Test Case 3 (3x3)
Matrix:
     2.000      4.000      6.000      8.000
     1.000      2.000      3.000      4.000
     3.000      6.000      9.000     12.000

L matrix:
     1.000      0.000      0.000
     0.500      1.000      0.000
     1.500      0.000      1.000

U matrix:
     2.000      4.000      6.000
     0.000      0.000      0.000
     0.000      0.000      0.000

Matrix A is singular (det A = 0)
-> System has infinite solutions.


Test Case 4 (2x2)
Matrix:
     1.000      1.000      2.000
     1.000      1.000      3.000

L matrix:
     1.000      0.000
     1.000      1.000

U matrix:
     1.000      1.000
     0.000      0.000

Matrix A is singular (det A = 0)
-> System has "NO" solution.
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Matrix Inversion Method

**Theory**
- Matrix inversion uses the formula A^-1 = (1/det(A)) * adj(A), where adj(A) is the adjoint (transpose of cofactor matrix).

- For each element, compute the cofactor C_ij = (-1)^(i+j) * M_ij, where M_ij is the minor (determinant of submatrix).

- Computational cost is O(n^4) for cofactor expansion; more efficient than Gauss-Jordan for small matrices but impractical for large ones.

### Algorithm steps:

1. **Calculate determinant of A** using cofactor expansion or recursive method.

2. **Check if det(A) ≠ 0**:  if determinant is zero or |det(A)| < EPS, matrix is singular; report failure.

3. **Compute cofactor matrix**:  for each element A[i][j], calculate cofactor C[i][j] = (-1)^(i+j) * det(M[i][j]).

4. **Find adjoint matrix**: transpose the cofactor matrix to get adj(A) = C^T.

5. **Compute inverse**:  calculate A^-1 = (1/det(A)) * adj(A) and solve x = A^-1 * b.

### What it does:

- Solves linear systems AX = b by computing A^-1 using determinant and adjoint matrix method.

- Handles multiple test cases, reading from Input.txt (matrix A and vector b per test case), writing to Output.txt, and printing to console.

- Reports the inverse matrix A^-1 and the solution vector x = A^-1 * b. 

- Detects singular matrices (det = 0) and reports failure.

**Code**

```cpp
#include <bits/stdc++.h>
using namespace std;
const double eps = 1e-9;

void getCofactor(vector<vector<double>> &A, vector<vector<double>> &temp, int p, int q, int n) {
    int i = 0, j = 0;
    for(int row = 0; row < n; row++) {
        for(int col = 0; col < n; col++) {
            if(row != p && col != q) {
                temp[i][j++] = A[row][col];
                if(j == n-1) { j = 0; i++; }
            }
        }
    }
}

double determinant(vector<vector<double>> &A, int n) {
    if(n == 1) return A[0][0];
    double det = 0;
    int sign = 1;
    vector<vector<double>> temp(n, vector<double>(n));
    for(int f = 0; f < n; f++) {
        getCofactor(A, temp, 0, f, n);
        det += sign * A[0][f] * determinant(temp, n-1);
        sign = -sign;
    }
    return det;
}

void adjoint(vector<vector<double>> &A, vector<vector<double>> &adj) {
    int n = A.size();
    if(n == 1) { adj[0][0] = 1; return; }
    int sign;
    vector<vector<double>> temp(n, vector<double>(n));
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            getCofactor(A, temp, i, j, n);
            sign = ((i+j)%2==0)? 1 : -1;
            adj[j][i] = sign * determinant(temp, n-1); 
        }
    }
}

vector<double> multiply(vector<vector<double>> &mat, vector<double> &v) {
    int n = mat.size();
    vector<double> res(n,0);
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            res[i]+=mat[i][j]*v[j];
    return res;
}

void printMatrix(ofstream &fout, vector<vector<double>> &mat) {
    int n = mat.size(), m = mat[0].size();
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++)
            fout << setw(10) << fixed << setprecision(3) << mat[i][j] << " ";
        fout << endl;
    }
    fout << endl;
}

void printVector(ofstream &fout, vector<double> &v) {
    for(double x:v)
        fout << setw(10) << fixed << setprecision(3) << x << endl;
    fout << endl;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    fout << "MATRIX INVERSION METHOD\n\n";
    int testCases; 
    fin >> testCases;

    for(int tc=1;tc<=testCases;tc++){
        int n;
        fin >> n;
        vector<vector<double>> A(n, vector<double>(n));
        vector<double> B(n);

        fout << "Test Case " << tc << " (" << n << "x" << n << ")\n\n";

        fout << "Coefficient Matrix A:\n";
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)
                fin >> A[i][j];
        printMatrix(fout, A);

        fout << "Vector b:\n";
        for(int i=0;i<n;i++) fin >> B[i];
        printVector(fout, B);

        double det = determinant(A,n);
        fout << "Determinant det(A): " << fixed << setprecision(3) << det << "\n\n";

        if(fabs(det)<eps){
            fout << "Matrix is singular (det=0)\n";
            fout << "-> System has no unique solution.\n\n";
            continue;
        }

        vector<vector<double>> adj(n, vector<double>(n));
        adjoint(A, adj);

        fout << "Adjoint Matrix adj(A):\n";
        printMatrix(fout, adj);

        vector<vector<double>> inv(n, vector<double>(n));
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)
                inv[i][j] = adj[i][j]/det;

        fout << "Inverse Matrix A^-1:\n";
        printMatrix(fout, inv);

        vector<double> x = multiply(inv,B);
        fout << "Solution (x = A^-1 * b):\n";
        for(int i=0;i<n;i++)
            fout << "x" << i+1 << " = " << fixed << setprecision(3) << x[i] << endl;
        fout << endl;
    }

    fin.close();
    fout.close();
    return 0;
}
```

**Input**

```bash
4
3
2 1 -1
-3 -1 2
-2 1 2
8
-11
-3
2
1 2
2 4
3
5
3
2 4 6
1 2 3
3 6 9
8
4
12
2
1 1
1 1
2
3
```
**Output**

```bash
MATRIX INVERSION METHOD 

-----------------------

Test Case 1 (3x3)

Coefficient Matrix A:
     2.000      1.000     -1.000
    -3.000     -1.000      2.000
    -2.000      1.000      2.000

Determinant det(A): 1.000

Adjoint Matrix adj(A):
     4.000      3.000     -1.000
    -2.000     -2.000      1.000
     5.000      4.000     -1.000

Inverse Matrix A⁻¹:
     4.000      3.000     -1.000
    -2.000     -2.000      1.000
     5.000      4.000     -1.000

Unique Solution exists.

Solution (x = A⁻¹ · b):
x1 =  2.000
x2 =  3.000
x3 = -1.000

--------------------------------------------------

Test Case 2 (2x2)

Coefficient Matrix A:
     1.000      2.000
     2.000      4.000

Determinant det(A): 0.000

Matrix is singular (det(A) = 0)
→ Inverse does not exist
→ System has NO unique solution

--------------------------------------------------

Test Case 3 (3x3)

Coefficient Matrix A:
     2.000      4.000      6.000
     1.000      2.000      3.000
     3.000      6.000      9.000

Determinant det(A): 0.000

Matrix is singular (det(A) = 0)
→ Inverse does not exist
→ System has INFINITE solutions

--------------------------------------------------

Test Case 4 (2x2)

Coefficient Matrix A:
     1.000      1.000
     1.000      1.000

Determinant det(A): 0.000

Matrix is singular (det(A) = 0)
→ Inverse does not exist
→ System has NO solution

--------------------------------------------------
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

## Differential Equation Solving

### Runge-Kutta 4th Order Method

**Theory**
This application numerically approximates solutions to ordinary differential equations utilizing the fourth-order Runge-Kutta technique. Processes multiple scenarios with varying starting values. 

### FUNCTION USED: 

```
dy/dx = f(x,y) = xy + y
```

### RUNGE-KUTTA 4TH ORDER FORMULA:

```
y(n+1) = y(n) + (h/6)[k1 + 2k2 + 2k3 + k4]
```

where: 

```
k1 = h*f(x(n), y(n))
k2 = h*f(x(n) + h/2, y(n) + k1/2)
k3 = h*f(x(n) + h/2, y(n) + k2/2)
k4 = h*f(x(n) + h, y(n) + k3)
```

### FEATURES:

- Batch processing capability for multiple test cases

- Constant step increment:  h = 0.001

- Dynamic iteration adjustment according to x interval

- Fourth-order accuracy guaranteed

- Controlled decimal precision in results

**Code**

```cpp
#include <bits/stdc++.h>
using namespace std;

double f(double x, double y) {
    return x*y + y;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    fout << "RUNGE KUTTA 4th ORDER METHOD\n\n";
    fout << "Function: f(x, y) = x*y + y\n\n";

    int T;
    fin >> T;
    fout << "Total Test Cases: " << T << "\n\n";

    for(int tc = 1; tc <= T; tc++) {
        double x0, y0, xn, h;
        fin >> x0 >> y0 >> xn >> h;

        int steps = (int)((xn - x0) / h);
        double x = x0;
        double y = y0;

        fout << "----------------------------------\n";
        fout << "TEST CASE #" << tc << "\n\n";
        fout << fixed << setprecision(3);
        fout << "Initial x0: " << x0 << ", y0: " << y0 << "\n";
        fout << "Final x: " << xn << "\n";
        fout << "Step size (h): " << h << "\n";
        fout << "Number of steps: " << steps << "\n\n";

        for(int i = 0; i < steps; i++) {
            double k1 = h * f(x, y);
            double k2 = h * f(x + h/2.0, y + k1/2.0);
            double k3 = h * f(x + h/2.0, y + k2/2.0);
            double k4 = h * f(x + h, y + k3);

            y = y + (k1 + 2*k2 + 2*k3 + k4)/6.0;
            x = x + h;
        }

        fout << "Result:\n";
        fout << "y(" << xn << ") = " << fixed << setprecision(3) << y << "\n\n";
    }

    fin.close();
    fout.close();
    return 0;
}
```

**Input**

```bash
3
0 1 1 0.001
0 2 2 0.001
1 1 3 0.001
```

**Output**

```bash
RUNGE KUTTA 4th ORDER METHOD

Function: f(x, y) = x*y + y

Total Test Cases: 3

----------------------------------
TEST CASE #1

Initial x0: 0.000, y0: 1.000
Final x: 1.000
Step size (h): 0.001
Number of steps: 1000

Result:
y(1.000) = 4.482

----------------------------------
TEST CASE #2

Initial x0: 0.000, y0: 2.000
Final x: 2.000
Step size (h): 0.001
Number of steps: 2000

Result:
y(2.000) = 109.196

----------------------------------
TEST CASE #3

Initial x0: 1.000, y0: 1.000
Final x: 3.000
Step size (h): 0.001
Number of steps: 2000

Result:
y(3.000) = 403.429
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

## Interpolation Methods

### Newton Forward Interpolation Method

**Theory**

#### Newton’s Forward Interpolation Method
Newton’s Forward Interpolation Method is a numerical technique used to estimate the value of a
function at a given point when the values of the function are known at equally spaced points. This
method is particularly suitable when the interpolation point lies near the beginning of the data set.
#### Basic Concept
Assume that the values of a function f(x) are known at n equally spaced points x0, x1, x2, ..., xn.
The spacing between consecutive x-values is constant and is denoted by h. Using these values, a
polynomial is constructed that passes through all the given points. This polynomial is then used to
estimate intermediate values of the function.
#### Forward Differences
Forward differences are used to build the interpolation polynomial. The first forward difference is
defined as:
∆y(i) = y(i+1) − y(i)
Higher-order forward differences are defined recursively as:
∆²y(i) = ∆y(i+1) − ∆y(i)
#### Newton’s Forward Interpolation Formula
The Newton forward interpolation polynomial is given by:
f(x) = y0 + u∆y0 + [u(u − 1) / 2!]∆²y0 + [u(u − 1)(u − 2) / 3!]∆³y0 + ...
where the parameter u is defined as:
u = (x − x0) / h
#### Advantages
- Simple and systematic method
- Efficient for equally spaced data
- Easy to implement in computer
programs
#### Applications
Newton’s Forward Interpolation Method is widely used in engineering, scientific computations,
numerical analysis, and data approximation problems where intermediate values need to be
estimated accurately


**Code**

```cpp
#include <bits/stdc++.h>
using namespace std;
#define arr vector<double>
double fun(const arr& coeff, double x) {
    double result = 0;
    int n = coeff.size(),i;
    for (i=0; i < n; i++) result += coeff[i] * pow(x, n - i - 1);
    return result;
}
int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    int t,cnt=1;
    fin>>t;
    while(t--){
    int degree,i,n,j;
    fin >> degree;
    arr coeff(degree + 1);
    for (i = 0; i <= degree; i++) fin >> coeff[i];
    fin >> n;
    arr x(n), y(n);
    for (i = 0; i < n; i++) {
        fin >> x[i];
        y[i] = fun(coeff, x[i]);
    }
    double value;
    fin >> value;
    vector<arr> diff(n, arr(n));
    for (i = 0; i < n; i++) diff[i][0] = y[i];
    for (j = 1; j < n; j++) {
        for (i = 0; i < n - j; i++) diff[i][j] = diff[i + 1][j - 1] - diff[i][j - 1];
    }
    double h = x[1] - x[0];
    double u = (value - x[0]) / h;
    double result = diff[0][0];
    double u2 = 1;
    double fact = 1;
    for ( i = 1; i < n; i++) {
        u2 *= (u - (i - 1));
        fact *= i;
        result += (u2 * diff[0][i]) / fact;
    }
    fout<<"Test Case: "<<cnt<<'\n';
    fout<<"Number of Data Points: "<<n<<'\n';
    fout<<"Data Points: (x,y):\n";
    for(i=0;i<n;i++) fout<<'('<<x[i]<<','<<y[i]<<")\n";
    fout<<"Forward Difference Table:\n";
    for(i=0;i<n;i++){
        fout<<"y["<<i+1<<"] = ";
        for(auto  &pk:diff[i]) fout<<pk<<' ';
        fout<<'\n';
    }
    fout << "Interpolated value at x = " << value << " is " << result << endl;
    cnt++;
    }
    fin.close();
    fout.close();
    return 0;
}
```

**Input**

```bash
1
3
1 -6 11 -6
4
0 1 2 3
1.5
```

**Output**

```bash
Test Case: 1
Number of Data Points: 4
Data Points: (x,y):
(0,-6)
(1,0)
(2,0)
(3,0)
Forward Difference Table:
y[1] = -6 6 -6 6 
y[2] = 0 0 0 0 
y[3] = 0 0 0 0 
y[4] = 0 0 0 0 
Interpolated value at x = 1.5 is 0.375

```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Newton Backward Interpolation Method

**Theory**

#### Newton’s Backward Interpolation Method
Newton’s Backward Interpolation Method is a numerical technique used to estimate the value of a
function when the known data points are equally spaced and the interpolation point lies near the
end of the data set. It is particularly useful when forward interpolation becomes less accurate.
#### Basic Concept
Let the values of a function f(x) be known at equally spaced points x0, x1, x2, ..., xn with a constant
spacing h. Newton’s backward interpolation constructs a polynomial using backward differences,
which are calculated from the end of the data table.
#### Backward Differences
Backward differences are defined as follows:
∇y(i) = y(i) − y(i−1)
Higher-order backward differences are defined recursively as:
∇²y(i) = ∇y(i) − ∇y(i−1)
#### Newton’s Backward Interpolation Formula
The Newton backward interpolation polynomial is given by:
f(x) = y(n) + u∇y(n) + [u(u + 1) / 2!]∇²y(n) + [u(u + 1)(u + 2) / 3!]∇³y(n) + ...
where the parameter u is defined as:
u = (x − x(n)) / h
#### Advantages
- Suitable when interpolation point is near the end of the data set
- Accurate for equally spaced data
- Easy to implement computationally
#### Applications
Newton’s Backward Interpolation Method is widely used in numerical analysis, engineering
problems, scientific computations, and data approximation tasks where estimation near the end of
tabulated data is required

**Code**

```cpp
#include <bits/stdc++.h>
using namespace std;
#define arr vector<double>
double fun(const arr& coeff, double x) {
    double result = 0;
    int n = coeff.size(),i;
    for (i = 0; i < n; i++) {
        result += coeff[i] * pow(x, n - i - 1);
    }
    return result;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    int t,cnt=1;
    fin>>t;
    while(t--){
    int degree,n,i,j;
    fin >> degree;
    arr coeff(degree + 1);
    for (i = 0; i <= degree; i++) fin >> coeff[i];
    fin >> n;
    arr x(n), y(n);
    for (i = 0; i < n; i++) {
        fin >> x[i];
        y[i] = fun(coeff, x[i]);
    }
    double value;
    fin >> value;
    vector<arr> diff(n, arr(n));
    for ( i = 0; i < n; i++) diff[i][0] = y[i];
    for ( j = 1; j < n; j++) {
        for (i = n - 1; i >= j; i--) {
            diff[i][j] = diff[i][j - 1] - diff[i - 1][j - 1];
        }
    }
    double h = x[1] - x[0];
    double u = (value - x[n - 1]) / h;
    double result = diff[n - 1][0];
    double u_term = 1;
    double fact = 1;
    for ( i = 1; i < n; i++) {
        u_term *= (u + (i - 1));
        fact *= i;
        result += (u_term * diff[n - 1][i]) / fact;
    }
    fout<<"Test Case: "<<cnt<<'\n';
    fout<<"Number of Data Points: "<<n<<'\n';
    fout<<"Data Points: (x,y):\n";
    for(i=0;i<n;i++) fout<<'('<<x[i]<<','<<y[i]<<")\n";
    fout<<"Backward Difference Table:\n";
    for(i=0;i<n;i++){
        fout<<"y["<<i+1<<"] = ";
        for(auto  &pk:diff[i]) fout<<pk<<' ';
        fout<<'\n';
    }
    fout << "Interpolated value at x = " << value << " is " << result << endl;
    cnt++;
    }
    fin.close();
    fout.close();
    return 0;
}

```

**Input**

```bash
1
3
1 -6 11 -6
4
0 1 2 3
2.5

```

**Output**

```bash
Test Case: 1
Number of Data Points: 4
Data Points: (x,y):
(0,-6)
(1,0)
(2,0)
(3,0)
Backward Difference Table:
y[1] = -6 0 0 0 
y[2] = 0 6 0 0 
y[3] = 0 0 -6 0 
y[4] = 0 0 0 6 
Interpolated value at x = 2.5 is -0.375

```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Newton Divided Difference Interpolation Method

**Theory**

#### Newton’s Divided Difference Interpolation Method
Newton’s Divided Difference Interpolation Method is a numerical technique used to construct an
interpolation polynomial when the given data points are not necessarily equally spaced. This
method generalizes Newton’s interpolation approach and is suitable for arbitrary data distributions.
#### Basic Concept
Suppose the values of a function f(x) are known at n distinct points x0, x1, x2, ..., xn. Newton’s
divided difference method constructs a polynomial that passes through all these points using
divided differences instead of finite differences.
#### Divided Differences
The first divided difference is defined as:
f[xi, xi+1] = ( f(xi+1) − f(xi) ) / ( xi+1 − xi )
Higher-order divided differences are defined recursively as:
f[xi, xi+1, ..., xi+k] = ( f[xi+1, ..., xi+k] − f[xi, ..., xi+k−1] ) / ( xi+k − xi )
#### Newton’s Divided Difference Formula
The interpolation polynomial is given by:
P(x) = f[x0] + (x − x0)f[x0, x1] + (x − x0)(x − x1)f[x0, x1, x2] + ...
#### Advantages
- Works for unequally spaced data
- Easy to extend when new data points are added 
- Computationally efficient
#### Applications
Newton’s Divided Difference Interpolation Method is widely used in numerical analysis, data fitting,
engineering computations, and scientific modeling where data points are irregularly spaced

**Code**

```cpp
#include <bits/stdc++.h>
using namespace std;
#define arr vector<double>
double fun(const arr& coeff, double x) {
    double result = 0;
    int n = coeff.size(),i;
    for ( i = 0; i < n; i++) result += coeff[i] * pow(x, n - i - 1);
    return result;
}
int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    int t,cnt=1;
    fin>>t;
    while(t--){
    int degree,n,i,j;
    fin >> degree;
    arr coeff(degree + 1);
    for ( i = 0; i <= degree; i++) fin >> coeff[i];
    fin >> n;
    arr x(n), y(n);
    for ( i = 0; i < n; i++){
        fin >> x[i];
        y[i] = fun(coeff, x[i]);
    }
    double value;
    fin >> value;
    vector<arr> div(n, arr(n));
    for (i = 0; i < n; i++) div[i][0] = y[i];
    for (j = 1; j < n; j++) {
        for (i = 0; i < n - j; i++) {
            div[i][j] = (div[i + 1][j - 1] - div[i][j - 1]) /
                        (x[i + j] - x[i]);
        }
    }
    double result = div[0][0];
    double term = 1.0;
    for (i = 1; i < n; i++) {
        term *= (value - x[i - 1]);
        result += term * div[0][i];
    }
    fout<<"Test Case: "<<cnt<<'\n';
    fout<<"Number of Data Points: "<<n<<'\n';
    fout<<"Data Points: (x,y):\n";
    for(i=0;i<n;i++) fout<<'('<<x[i]<<','<<y[i]<<")\n";
    fout<<"Difference Table:\n";
    for(i=0;i<n;i++){
        fout<<"y["<<i+1<<"] = ";
        for(auto  &pk:div[i]) fout<<pk<<' ';
        fout<<'\n';
    }
    fout << "Interpolated value at x = " << value << " is " << result << endl;
    cnt++;
    }
    fin.close();
    fout.close();
    return 0;
}

```

**Input**

```bash
1
3
1 -6 11 -6
4
0 1.2 2.1 3.5
2.5

```

**Output**

```bash
Test Case: 1
Number of Data Points: 4
Data Points: (x,y):
(0,-6)
(1.2,0.288)
(2.1,-0.099)
(3.5,1.875)
Difference Table:
y[1] = -6 5.24 -2.7 1 
y[2] = 0.288 -0.43 0.8 0 
y[3] = -0.099 1.41 0 0 
y[4] = 1.875 0 0 0 
Interpolated value at x = 2.5 is -0.375

```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

## Numerical Differentiation

### Differentiation by Forward Interpolation Method

**Theory**
Newton's Forward Interpolation serves as a computational approach for approximating function values along with their derivatives through a collection of discrete coordinate points. When a function f(x)f(x)f(x) is evaluated at uniformly spaced coordinates x0,x1,... ,xn, the derivative values at point X,X,X can be estimated using Forward Difference Table.

**• First Derivative f'(x)**

```
Δ1yi = yi+1-yi
Δ2yi = Δ1yi+1-Δ1yi
Δ3yi = Δ2yi+1-Δ2yi
```

### Where

```
yi = f(xi)
```

### Derivative Formulas: 

```
f'(X) ≈ (y0 + (2u-1)Δ2y0/2!  + (3u2-6u+2)Δ3y0/3!  +--)/h
```

**• Second Derivative**

```
f''(X) ≈ (Δ2y0+(u -1)Δ3y0--) /h^2
```

### Where: 

```
u=(X-x0)/h
h=x(i+1)-xi
```

### Error Calculation: 

Numerical derivatives are validated against exact analytical derivatives f'(X) and f''(X):

```
Error = |Analytical-Numerical|/Analytical×100
```

## Algorithm / Steps

1. Extract total number of test scenarios T.

2. Process each test scenario:

• Extract n,a,b,X from input source.

• Determine step increment h = (b -a)/n

• Create uniformly distributed coordinates xi = a + i ∗ h and evaluate yi = f(xi).

• Build forward difference table structure.

• Evaluate first derivative f'(X)) applying Newton's approach.

• Evaluate second derivative f''(X) applying Newton's approach.

• Validate against exact derivatives to determine percentage deviation.

Display inputs, difference table, computed derivatives, and error metrics to terminal and output document.

### Features

• Handles multiple test scenarios

• Computes first and second derivatives through numerical methods

• Constructs complete forward difference table

• Determines percentage deviation from analytical derivatives

• Writes results to terminal and file simultaneously


**Code**

```cpp
#include <bits/stdc++.h>
using namespace std;

double trueDerivative(vector<double> &coeffs, double x) {
    double result = 0;
    int n = coeffs.size();
    for (int i = 1; i < n; i++) {
        result += i * coeffs[i] * pow(x, i - 1);
    }
    return result;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int testCases;
    fin >> testCases;

    for (int tc = 1; tc <= testCases; tc++) {
        int n, deg;
        fin >> deg >> n;
        vector<double> coeffs(deg + 1);
        for (int i = 0; i <= deg; i++)
            fin >> coeffs[i];

        vector<double> x(n), y(n);
        for (int i = 0; i < n; i++)
            fin >> x[i];
        for (int i = 0; i < n; i++) {
            double val = 0;
            for (int j = 0; j <= deg; j++)
                val += coeffs[j] * pow(x[i], j);
            y[i] = val;
        }

        double h = x[1] - x[0];
        double diffPoint;
        fin >> diffPoint;

        vector<vector<double>> F(n, vector<double>(n));
        for (int i = 0; i < n; i++)
            F[i][0] = y[i];

        for (int j = 1; j < n; j++) {
            for (int i = 0; i < n - j; i++)
                F[i][j] = F[i+1][j-1] - F[i][j-1];
        }

        double u = (diffPoint - x[0]) / h;
        double derivative = F[0][1] 
                          + ((2*u - 1)/2.0) * F[0][2]
                          + ((3*u*u - 6*u + 2)/6.0) * F[0][3]
                          + ((4*u*u*u - 18*u*u + 22*u - 6)/24.0) * F[0][4]; // up to 4th order
                          
        double trueDeriv = trueDerivative(coeffs, diffPoint);
        double error = fabs((trueDeriv - derivative) / trueDeriv) * 100;

        fout << "Test Case " << tc << "\n";
        fout << "Polynomial: f(x) = ";
        for(int i=0;i<=deg;i++){
            fout << fixed << setprecision(4) << coeffs[i];
            if(i>0) fout << "x^" << i;
            if(i<deg) fout << " + ";
        }
        fout << "\nNumber of points: " << n << "\n";
        fout << "x-values: ";
        for(double xi:x) fout << fixed << setprecision(6) << xi << " ";
        fout << "\n";
        fout << "y-values: ";
        for(double yi:y) fout << fixed << setprecision(6) << yi << " ";
        fout << "\n";
        fout << "Step size (h): " << fixed << setprecision(6) << h << "\n";
        fout << "Differentiation point: " << fixed << setprecision(6) << diffPoint << "\n";

        fout << "Forward Difference Table:\n";
        for(int i=0;i<n;i++){
            fout << "Row " << i << ": ";
            for(int j=0;j<n-i;j++)
                fout << fixed << setprecision(6) << F[i][j] << " ";
            fout << "\n";
        }

        fout << "Approximate derivative : " << fixed << setprecision(6) << derivative << "\n";
        fout << "True derivative: " << fixed << setprecision(6) << trueDeriv << "\n";
        fout << "Percentage error: " << fixed << setprecision(6) << error << " %\n\n";
    }

    fin.close();
    fout.close();
    return 0;
}
```

**Input**

```bash
2
3 6
1 2 -1 0.5
0 1 2 3 4 5
4.5
2 7
0 0 1
0 0.5 1 1.5 2 2.5 3
2.0
```

**Output**

```bash
Test Case 1
Polynomial: f(x) = 1.0000 + 2.0000x^1 + -1.0000x^2 + 0.5000x^3
Number of points: 6
x-values: 0.000000 1.000000 2.000000 3.000000 4.000000 5.000000 
y-values: 1.000000 2.500000 5.000000 11.500000 25.000000 48.500000 
Step size (h): 1.000000
Differentiation point: 4.500000
Forward Difference Table:
Row 0: 1.000000 1.500000 1.000000 0.500000 0.500000 0.000000 
Row 1: 2.500000 2.500000 1.500000 1.000000 0.500000 
Row 2: 5.000000 6.500000 4.000000 1.500000 
Row 3: 11.500000 13.500000 7.000000 
Row 4: 25.000000 23.500000 
Row 5: 48.500000 
Approximate derivative : 23.375000
True derivative: 23.375000
Percentage error: 0.000000 %

Test Case 2
Polynomial: f(x) = 0.0000 + 0.0000x^1 + 1.0000x^2
Number of points: 7
x-values: 0.000000 0.500000 1.000000 1.500000 2.000000 2.500000 3.000000 
y-values: 0.000000 0.250000 1.000000 2.250000 4.000000 6.250000 9.000000 
Step size (h): 0.500000
Differentiation point: 2.000000
Forward Difference Table:
Row 0: 0.000000 0.250000 0.500000 0.500000 0.500000 0.500000 0.500000 
Row 1: 0.250000 0.750000 0.500000 0.000000 0.000000 0.000000 
Row 2: 1.000000 1.250000 0.500000 0.000000 0.000000 
Row 3: 2.250000 1.750000 0.500000 0.000000 
Row 4: 4.000000 2.250000 0.500000 
Row 5: 6.250000 2.750000 
Row 6: 9.000000 
Approximate derivative : 4.000000
True derivative: 4.000000
Percentage error: 0.000000 %
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Differentiation by Backward Interpolation Method

**Theory**

This approach employs Newton's backward difference technique to evaluate derivatives at positions close to the final entries of a dataset. It generates backward differences and utilizes Newton's backward differentiation expressions.

### BACKWARD DIFFERENCE DERIVATIVE FORMULAS: 

```
f'(xₙ) = [∇f(xₙ) + (2s+1)∇²f(xₙ)/2!  + (3s²+6s+2)∇³f(xₙ)/6 +... ] / h
```

### where

```
s = (X - xₙ)/h
```

Applied for polynomial derivative computation using backward differences.


**Code**

```cpp
#include <bits/stdc++.h>
using namespace std;

double trueDerivative(vector<double> &coeffs, double x) {
    double result = 0;
    int n = coeffs.size();
    for (int i = 1; i < n; i++)
        result += i * coeffs[i] * pow(x, i - 1);
    return result;
}

double factorial(int n) {
    double f = 1;
    for(int i=2;i<=n;i++) f *= i;
    return f;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int testCases;
    fin >> testCases;

    for(int tc=1; tc<=testCases; tc++){
        int deg, n;
        fin >> deg >> n;

        vector<double> coeffs(deg+1);
        for(int i=0;i<=deg;i++) fin >> coeffs[i];

        vector<double> x(n), y(n);
        for(int i=0;i<n;i++) fin >> x[i];

        for(int i=0;i<n;i++){
            double val = 0;
            for(int j=0;j<=deg;j++)
                val += coeffs[j]*pow(x[i], j);
            y[i] = val;
        }

        double h = x[1]-x[0];
        double diffPoint;
        fin >> diffPoint;

        vector<vector<double>> B(n, vector<double>(n));
        for(int i=0;i<n;i++) B[i][0] = y[i];
        for(int j=1;j<n;j++)
            for(int i=j;i<n;i++)
                B[i][j] = B[i][j-1] - B[i-1][j-1];

        int m = n-1;
        double u = (diffPoint - x[m])/h;

        double derivative = B[m][1];
        derivative += ((2*u+1)/2.0)*B[m][2];
        derivative += ((3*u*u + 6*u + 2)/6.0)*B[m][3];
        derivative += ((4*u*u*u + 18*u*u + 22*u + 6)/24.0)*B[m][4];

        double trueDeriv = trueDerivative(coeffs, diffPoint);
        double error = fabs((trueDeriv - derivative)/trueDeriv)*100;

        fout << "Test Case " << tc << "\n";
        fout << "Polynomial: f(x) = ";
        for(int i=0;i<=deg;i++){
            fout << fixed << setprecision(4) << coeffs[i];
            if(i>0) fout << "x^" << i;
            if(i<deg) fout << " + ";
        }
        fout << "\nNumber of points: " << n << "\n";

        fout << "x-values: ";
        for(double xi:x) fout << fixed << setprecision(6) << xi << " ";
        fout << "\n";

        fout << "y-values: ";
        for(double yi:y) fout << fixed << setprecision(6) << yi << " ";
        fout << "\n";

        fout << "Step size (h): " << fixed << setprecision(6) << h << "\n";
        fout << "Differentiation point: " << fixed << setprecision(6) << diffPoint << "\n";

        fout << "Backward Difference Table:\n";
        for(int i=0;i<n;i++){
            fout << "Row " << i << ": ";
            for(int j=0;j<=i;j++)
                fout << fixed << setprecision(6) << B[i][j] << " ";
            fout << "\n";
        }

        fout << "Approximate derivative : " << fixed << setprecision(6) << derivative << "\n";
        fout << "True derivative: " << fixed << setprecision(6) << trueDeriv << "\n";
        fout << "Percentage error: " << fixed << setprecision(6) << error << " %\n\n";
    }

    fin.close();
    fout.close();
    return 0;
}
```

**Input**

```bash
2
3 6
1 2 -1 0.5
0 1 2 3 4 5
4.5
2 7
0 0 1
0 0.5 1 1.5 2 2.5 3
2.0
```

**Output**

```bash
Test Case 1
Polynomial: f(x) = 1.0000 + 2.0000x^1 + -1.0000x^2 + 0.5000x^3
Number of points: 6
x-values: 0.000000 1.000000 2.000000 3.000000 4.000000 5.000000 
y-values: 1.000000 2.500000 5.000000 11.500000 25.000000 48.500000 
Step size (h): 1.000000
Differentiation point: 4.500000
Backward Difference Table:
Row 0: 1.000000 
Row 1: 2.500000 1.500000 
Row 2: 5.000000 2.500000 1.000000 
Row 3: 11.500000 6.500000 4.000000 3.000000 
Row 4: 25.000000 13.500000 7.000000 3.000000 0.000000 
Row 5: 48.500000 23.500000 10.000000 3.000000 0.000000 0.000000 
Approximate derivative : 23.375000
True derivative: 23.375000
Percentage error: 0.000000 %

Test Case 2
Polynomial: f(x) = 0.0000 + 0.0000x^1 + 1.0000x^2
Number of points: 7
x-values: 0.000000 0.500000 1.000000 1.500000 2.000000 2.500000 3.000000 
y-values: 0.000000 0.250000 1.000000 2.250000 4.000000 6.250000 9.000000 
Step size (h): 0.500000
Differentiation point: 2.000000
Backward Difference Table:
Row 0: 0.000000 
Row 1: 0.250000 0.250000 
Row 2: 1.000000 0.750000 0.500000 
Row 3: 2.250000 1.250000 0.500000 0.000000 
Row 4: 4.000000 1.750000 0.500000 0.000000 0.000000 
Row 5: 6.250000 2.250000 0.500000 0.000000 0.000000 0.000000 
Row 6: 9.000000 2.750000 0.500000 0.000000 0.000000 0.000000 0.000000 
Approximate derivative : 4.000000
True derivative: 4.000000
Percentage error: 0.000000 %
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

## Curve Fitting & Regression

### Linear Regression Method

**Theory**
Linear regression is a numerical technique used to establish a linear relationship between two
variables based on observed data. It is one of the simplest and most widely used curve fitting
methods in numerical analysis.
In linear regression, one variable is considered independent and is denoted by x, while the other
variable is dependent and is denoted by y. The objective is to find a straight line that best
represents the relationship between x and y.
The mathematical model of linear regression is given by y = a + bx, where a is the intercept and b is
the slope of the line. These constants are chosen such that the line fits the given data points as
closely as possible.
The method of least squares is used to determine the values of a and b. According to this method,
the sum of the squares of the deviations between the observed values and the estimated values is
minimized.
By applying the least squares principle, a set of normal equations is obtained. Solving these
equations yields the best-fit values of the regression coefficients a and b.
Linear regression is commonly used in engineering and scientific applications to model trends,
make predictions, and analyze experimental data.
#### CONCEPT

For n data points (xi, yi) the regression line: y = a + bx is determined by minimizing:
∑(yi − (a + bxi))²

This leads to two closed-form formulas:

Slope (b):
b = [ n∑xy − (∑x)(∑y) ] / [ n∑x² − (∑x)² ]

Intercept (a):
a = (∑y − b∑x) / n

These values define the best-fit straight line.

ALGORITHM (Least Squares Method):

• Read number of data points n  
• Read arrays x[n] and y[n]  
• Compute required sums:
  - Σx
  - Σy
  - Σxy
  - Σx²

**Code**

```cpp
#include <bits/stdc++.h>
using namespace std;
#define arr vector<double>
int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    int t,cnt=1;
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
    fout<<"Test Case: "<<cnt<<'\n';
    fout<<"x - values: ";
    for(i=0;i<n;i++) fout<<x[i]<<' ';
    fout<<'\n';
    fout<<"y - values: ";
    for(i=0;i<n;i++) fout<<y[i]<<' ';
    fout<<'\n';
    fout<<"Linear Regression Line:\n";
    fout<<"y = a + b x\n";
    fout<<"a = "<<a<<'\n';
    fout<<"b = "<<b<< '\n';
    cnt++;
    }
    fin.close();
    fout.close();
    return 0;
}

```

**Input**

```bash
1
5
1 2
2 4
3 5
4 4
5 6
```

**Output**

```bash
Test Case: 1
x - values: 1.000 2.000 3.000 4.000 5.000 
y - values: 2.000 4.000 5.000 4.000 6.000 
Linear Regression Line:
y = a + b x
a = 1.800
b = 0.800

```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Polynomial Regression Method

**Theory**

Polynomial regression is a curve fitting technique used when the relationship between the variables
is nonlinear. It extends linear regression by fitting a polynomial of higher degree to the given data
points.
In polynomial regression, the dependent variable y is expressed as a polynomial function of the
independent variable x. The general form of the polynomial regression model is y = a0 + a1x +
a2x^2 + ... + amx^m.
The degree of the polynomial m is chosen by the user based on the nature of the data. A higher
degree polynomial can provide a closer fit but may also lead to overfitting.
The coefficients of the polynomial are determined using the method of least squares. Normal
equations are formed by minimizing the sum of squared errors between the observed and predicted
values.
These normal equations form a system of linear equations which can be solved using numerical
methods such as Gaussian elimination.
Polynomial regression is useful in cases where experimental data exhibits curvature and cannot be
accurately modeled using a straight line.
Extends linear regression to fit polynomial curves of degree m:
- y = a₀ + a₁x + a₂x² + ... + aₘxᵐ
  
The coefficients are found by solving the normal equations: 
- ∑xy = a₀∑x + a₁∑x² + a₂∑x³ + ... + aₘ∑xᵐ⁺¹ ... ∑xᵐy = a₀∑xᵐ + a₁∑xᵐ⁺¹ + ... + aₘ∑x²ᵐ

**Code**

```cpp
#include <bits/stdc++.h>
using namespace std;
#define arr vector<double>
int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    int t,cnt=1;
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
    fout<<"Test Case: "<<cnt<<'\n';
    fout<<"x - values: ";
    for(i=0;i<n;i++) fout<<x[i]<<' ';
    fout<<'\n';
    fout<<"y - values: ";
    for(i=0;i<n;i++) fout<<y[i]<<' ';
    fout<<'\n';
    fout<<"Polynomial Regression (degree "<<m<<"):\n";
    for(i=0;i<=m;i++)fout<<'a'<<i<<" = "<<a[i]<<'\n';
    cnt++;
    }
    fin.close();
    fout.close();
    return 0;
}

```

**Input**

```bash
2
5 2
0 1
1 6
2 17
3 34
4 57
6 3
-2 -4
-1 0
0 2
1 2
2 8
3 26

```

**Output**

```bash
Test Case: 1
x - values: 0.000 1.000 2.000 3.000 4.000 
y - values: 1.000 6.000 17.000 34.000 57.000 
Polynomial Regression (degree 2):
a0 = 1.000
a1 = 2.000
a2 = 3.000
Test Case: 2
x - values: -2.000 -1.000 0.000 1.000 2.000 3.000 
y - values: -4.000 0.000 2.000 2.000 8.000 26.000 
Polynomial Regression (degree 3):
a0 = 1.175
a1 = -0.307
a2 = 0.230
a3 = 0.870

```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Transcendental Regression Method

**Theory**

Transcendental regression is used when the relationship between variables involves exponential or
power functions rather than simple polynomials. Such models frequently arise in physical,
biological, and engineering systems.
One common transcendental model is y = ae^(bx), which represents exponential growth or decay.
This model is transformed into a linear form by taking the natural logarithm of both sides.
Another important model is y = ax^b, known as the power model. This model is also linearized
using logarithmic transformation, allowing the coefficients to be estimated using least squares.
A third transcendental model is y = a + be^(x/z), where z is a known constant provided by the user.
In this case, e^(x/z) is treated as an independent variable.
After suitable transformation or substitution, the method of least squares is applied to estimate the
unknown parameters of the model.
Transcendental regression provides an effective way to fit complex nonlinear relationships using
linear regression techniques.

**Code**

```cpp
#include <bits/stdc++.h>
using namespace std;
#define arr vector<double>
int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    int t,cnt=1;
    fin>>t;
    while(t--){
    fout<<"Test Case: "<<cnt<<'\n';
    cnt++;
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
        fout<<"Model: y=ae^(bx)\n";
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
        fout<<"Model: y=ax^b\n";
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
        fout<<"Model: y = a + b e^(x/z)\n";
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

```

**Input**

```bash
3
1 5
1 2 3 4 5
2.7 7.4 20.1 54.6 148.4
2 4
2 3 4 5
8 27 64 125
3 6
0 4 8 12 16 20
6.0 13.4 29.6 65.1 143.5 316.0
4

```

**Output**

```bash
Test Case: 1
Model: y=ae^(bx)
a = 3.311
b = 0.072
Equation: y = 3.311 e^(0.072 x)
Test Case: 2
Model: y=ax^b
a = 1.522
b = 1.100
Equation: y = 1.522 x^1.100
Test Case: 3
Model: y = a + b e^(x/z)
a = 22.900
b = 0.000
Equation: y = 22.900 + 0.000 e^(x/4.000)

```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

## Numerical Integration

### Simpson's 1/3 Rule Method

**Theory**

-Simpson’s 1/3 rule is a numerical integration method used to approximate definite integrals.
-The interval is divided into an even number of equal subintervals, and the integrand is approximated by parabolic arcs.
-It uses function values at equally spaced points with weights 1, 4, and 2 alternately.
-The method provides higher accuracy than the trapezoidal rule for smooth functions.

**Code**

```cpp
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
        if(n%2!=0){
            fout<<"Invalid n (must be even)\n";
            continue;
        }
        h=(b-a)/n;
        sum=f(a)+f(b);
        for(int i=1;i<n;i++){
            double x=a+i*h;
            if(i%2==0) sum+=2*f(x);
            else sum+=4*f(x);
        }
        double result=(h/3)*sum;
        fout<<fixed<<setprecision(6);
        fout<<"Integral: "<<result<<'\n';
    }
    fin.close();
    fout.close();
    return 0;
}
```

**Input**

```bash
2
2
1 0 0
0 1
10
3
1 0 -1 -2
1 2
10
```

**Output**

```bash
Integral: 0.333333
Integral: 2.25
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Simpson's 3/8 Rule Method

**Theory**

-Simpson’s 3/8 rule is a numerical integration technique used to evaluate definite integrals.
-The interval is divided into subintervals in multiples of three, and the function is approximated using cubic polynomials.
-It uses function values with weights 1, 3, 3, and 1 over each group of three subintervals.
-The method is slightly less accurate than Simpson’s 1/3 rule but useful when the number of subintervals is a multiple of three.

**Code**

```cpp
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
```

**Input**

```bash
3
2
1 0 0
0 1
6
3
1 0 -1 -2
1 2
12
1
1 1
0 3
9
```

**Output**

```bash
Integral: 0.333333
Integral: 2.25
Integral: 7.5
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

## Purpose of This Project

This project was developed as a requirement of the Numerical Methods Laboratory course(CSE-2208) to:

- Develop C++ implementations of classical numerical methods
- Explore the practical significance of numerical analysis
- Apply file handling and modular programming concepts
- Establish a reusable numerical methods library for solving practical problems

---

## Authors

![](https://img.shields.io/badge/1-Mahabubul%20Alam%20|%20Roll:%202207001-ECC5C0?style=for-the-badge)

![](https://img.shields.io/badge/2-Md%20Sulaiman%20|%20Roll:%202207011-FFCC00?style=for-the-badge)

![](https://img.shields.io/badge/3-Md%20Redoanul%20Karim%20Rusho%20|%20Roll:%202207014-CD7F32?style=for-the-badge)

---

