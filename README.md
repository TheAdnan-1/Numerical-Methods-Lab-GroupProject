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

```bash
cat Bisection/Bisection.txt
```

**Code**

```bash
cat Bisection/Bisection.cpp
```

**Input**

```bash
cat Bisection/input.txt
```

**Output**

```bash
cat Bisection/output.txt
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### False Position Method

**Theory**

```bash
cat "False Position/False Position.txt"
```

**Code**

```bash
cat "False Position/FalsePosition.cpp"
```

**Input**

```bash
cat "False Position/input.txt"
```

**Output**

```bash
cat "False Position/output.txt"
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Newton Raphson Method

**Theory**

```bash
cat "Newton Raphson/Newton Raphson.txt"
```

**Code**

```bash
cat "Newton Raphson/NewtonRaphson. cpp"
```

**Input**

```bash
cat "Newton Raphson/input.txt"
```

**Output**

```bash
cat "Newton Raphson/output.txt"
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Secant Method

**Theory**

```bash
cat Secant/Secant.txt
```

**Code**

```bash
cat Secant/Secant.cpp
```

**Input**

```bash
cat Secant/input.txt
```

**Output**

```bash
cat Secant/output.txt
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

## Solution of Linear Equations

### Gauss Elimination Method

**Theory**

```bash
cat "Gauss Elimination/Gauss Elimination.txt"
```

**Code**

```bash
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

```bash
cat "Gauss Jordan Elimination/Gauss Jordan Elimination.txt"
```

**Code**

```bash
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

```bash
cat "LU Decomposition/LU Decomposition. txt"
```

**Code**

```bash
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

```bash
cat "Matrix Inversion/Matrix Inversion.txt"
```

**Code**

```bash
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

```bash
cat "Runge Kutta 4th Order/Runge Kutta 4th Order.txt"
```

**Code**

```bash
cat "Runge Kutta 4th Order/RungeKutta. cpp"
```

**Input**

```bash
cat "Runge Kutta 4th Order/input.txt"
```

**Output**

```bash
cat "Runge Kutta 4th Order/output.txt"
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

## Interpolation Methods

### Newton Forward Interpolation Method

**Theory**

```bash
cat "Newton Forward Interpolation/Newton Forward Interpolation.txt"
```

**Code**

```bash
cat "Newton Forward Interpolation/NewtonForward.cpp"
```

**Input**

```bash
cat "Newton Forward Interpolation/input.txt"
```

**Output**

```bash
cat "Newton Forward Interpolation/output.txt"
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Newton Backward Interpolation Method

**Theory**

```bash
cat "Newton Backward Interpolation/Newton Backward Interpolation.txt"
```

**Code**

```bash
cat "Newton Backward Interpolation/NewtonBackward.cpp"
```

**Input**

```bash
cat "Newton Backward Interpolation/input.txt"
```

**Output**

```bash
cat "Newton Backward Interpolation/output.txt"
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Newton Divided Difference Interpolation Method

**Theory**

```bash
cat "Newton Divided Difference Interpolation/Newton Divided Difference Interpolation.txt"
```

**Code**

```bash
cat "Newton Divided Difference Interpolation/NewtonDividedDifference.cpp"
```

**Input**

```bash
cat "Newton Divided Difference Interpolation/input.txt"
```

**Output**

```bash
cat "Newton Divided Difference Interpolation/output.txt"
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

## Numerical Differentiation

### Differentiation by Forward Interpolation Method

**Theory**

```bash
cat "Differentiation by Forward Interpolation/Differentiation by Forward Interpolation.txt"
```

**Code**

```bash
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

```bash
cat "Differentiation by Backward Interpolation/Differentiation by Backward Interpolation. txt"
```

**Code**

```bash
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

```bash
cat "Linear Regression/Linear Regression.txt"
```

**Code**

```bash
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
```

**Input**

```bash
cat "Linear Regression/input.txt"
```

**Output**

```bash
cat "Linear Regression/output.txt"
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Polynomial Regression Method

**Theory**

```bash
cat "Polynomial Regression/Polynomial Regression.txt"
```

**Code**

```bash
cat "Polynomial Regression/PolynomialRegression.cpp"
```

**Input**

```bash
cat "Polynomial Regression/input.txt"
```

**Output**

```bash
cat "Polynomial Regression/output.txt"
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Transcendental Regression Method

**Theory**

```bash
cat "Transcendental Regression/Transcendental Regression.txt"
```

**Code**

```bash
cat "Transcendental Regression/TranscendentalRegression.cpp"
```

**Input**

```bash
cat "Transcendental Regression/input.txt"
```

**Output**

```bash
cat "Transcendental Regression/output.txt"
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

## Numerical Integration

### Simpson's 1/3 Rule Method

**Theory**

```bash
cat "Simpson 1_3/Simpson 1_3.txt"
```

**Code**

```bash
cat "Simpson 1_3/Simpson13.cpp"
```

**Input**

```bash
cat "Simpson 1_3/input.txt"
```

**Output**

```bash
cat "Simpson 1_3/output.txt"
```

[⬆ Back to top](#numerical-methods-laboratory-group-project)

---

### Simpson's 3/8 Rule Method

**Theory**

```bash
cat "Simpson 3_8/Simpson 3_8.txt"
```

**Code**

```bash
cat "Simpson 3_8/Simpson38.cpp"
```

**Input**

```bash
cat "Simpson 3_8/input.txt"
```

**Output**

```bash
cat "Simpson 3_8/output.txt"
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

![](https://img.shields.io/badge/3-Md%20Redoyanul%20Karim%20Rusho%20|%20Roll:%202207014-CD7F32?style=for-the-badge)

---

