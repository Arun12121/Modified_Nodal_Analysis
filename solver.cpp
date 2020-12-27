#include <string>
#include <vector>
#include <emscripten/bind.h>

using namespace std;

vector<string> split(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}

/*
*R- Resistor
*I - Current Source
*V - Voltage Source
*G - VCCS
*E - VCVS
*F - CCCS
*H - CCVS
*/

class Circuit {
    private:
    double **cond;
    double *curr;
    int n, rn, vSrcCurrent;

    public:
    Circuit(int N, int RN) 
    {
        int i;
        /*A = calloc(n, sizeof(int *));
        R = calloc(rn, sizeof(int *));
        B = calloc(n, sizeof(int));
        for(i = 0; i<n; i++)
            A[i] = calloc(n, sizeof(int));
        for(i = 0; i<rn; i++)
            R[i] = calloc(rn, sizeof(int));*/
        n = N;
        rn = RN;
        vSrcCurrent = rn;
        cond = new double*[n];
        curr = new double[n];

        for(i = 0; i<n; i++) {
            cond[i] = new double[n];
            curr[i] = 0;
        }
        for(i = 0; i<n; i++)
            for(int j = 0; j < n; j++)
                cond[i][j] = 0;
    }

    void addResistor(int k, int l, double R) 
    {
        if(k != 0)
            cond[k-1][k-1] += 1/R;
        if(l != 0)
            cond[l-1][l-1] += 1/R;
        if(l != 0 && k != 0) {
            cond[l-1][k-1] -= 1/R;
            cond[k-1][l-1] -= 1/R;
        }
    }

    void addCurrentSource(int k, int l, double I)
    {
        if(k != 0)
            curr[k-1] += I;
        if(l != 0)
            curr[l-1] -= I;
    }

    void addVoltageSource(int k, int l, double V)
    {
        if(k != 0) {
            cond[k-1][vSrcCurrent] += 1;
            cond[vSrcCurrent][k-1] += 1;
        }
        if(l != 0) {
            cond[l-1][vSrcCurrent] -= 1;
            cond[vSrcCurrent][l-1] -= 1;
        }
        curr[vSrcCurrent] += V;
        vSrcCurrent++;
    }

    void addVCCS(int k, int l, int m, int n, double G)
    {
        if(m != 0)
        {
            if(l != 0) cond[l-1][m-1] -= G;
            if(k != 0) cond[k-1][m-1] += G;
        }
        if(n != 0)
        {
            if(l != 0) cond[l-1][n-1] += G;
            if(k != 0) cond[k-1][n-1] -= G;
        }
    }

    void addVCVS(int k, int l, int m, int n, double A)
    {
        if(k != 0)
        {
            cond[k-1][vSrcCurrent] += 1;
            cond[vSrcCurrent][k-1] += 1;
        }
        if(l != 0)
        {
            cond[l-1][vSrcCurrent] -= 1;
            cond[vSrcCurrent][l-1] -= 1;
        }
        if(m != 0) cond[vSrcCurrent][m-1] -= A;
        if(n != 0) cond[vSrcCurrent][n-1] += A;
    }

    void addCCCS(int k, int l, int m, int n, double A)
    {
        if(m != 0)
        {
            cond[m-1][vSrcCurrent] += 1;
            cond[vSrcCurrent][m-1] += 1;
        }
        if(n != 0)
        {
            cond[n-1][vSrcCurrent] -= 1;
            cond[vSrcCurrent][n-1] -= 1;
        }
        if(k != 0) cond[k-1][vSrcCurrent] += A;
        if(l != 0) cond[l-1][vSrcCurrent] -= A;
    }

    void addCCVS(int k, int l, int m, int n, double A)
    {
        if(m != 0)
        {
            cond[m-1][vSrcCurrent] += 1;
            cond[vSrcCurrent][m-1] += 1;
        }
        if(n != 0)
        {
            cond[n-1][vSrcCurrent] -= 1;
            cond[vSrcCurrent][n-1] -= 1;
        }
        vSrcCurrent++;
        if(k != 0)
        {
            cond[k-1][vSrcCurrent] += 1;
            cond[vSrcCurrent][k-1] += 1;
        }
        if(l != 0)
        {
            cond[l-1][vSrcCurrent] -= 1;
            cond[vSrcCurrent][l-1] -= 1;
        }
        cond[vSrcCurrent][vSrcCurrent-1] -= A;
        vSrcCurrent++;
    }

    void compute(double X[])
    {
        //LU Decomposition
        double L[n][n], U[n][n], Y[n];
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                L[i][j] = U[i][j] = 0;
            }
            X[i] = Y[i] = 0;
        }
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                if(i==j)
                    L[i][j]=1;
                else if(i<j)
                    L[i][j]=0;
                U[i][j]=cond[i][j];
            }
        }
        for(int j=0;j<n;j++)
        {
            for(int i=j+1;i<n;i++)
            {
                // Find Coeff
                double r = U[i][j]/U[j][j];
                for(int p=0;p<n;p++)
                    U[i][p] = U[i][p] - r*U[j][p];
                L[i][j] = r;
            }
        }
        
        //finding solution using Forward and Backward method
        for(int i=0; i<n; i++)
        {
            Y[i]=curr[i];
            for(int j=0; j<i; j++)
            {
                Y[i]-=L[i][j]*Y[j];
            }
        }
        
        for(int i=n-1; i>=0; i--)
        {
            X[i]= Y[i];
            for(int j=i+1; j<n; j++)
            {
                X[i]-=U[i][j]*X[j];
            }
            X[i]/=U[i][i];
        }
    }


    ~Circuit() 
    {

        for(int i = 0; i<n; i++)
        {
            free(cond[i]);
            cond[i] = NULL;
        }
        free(cond); cond = NULL; 
        free(curr); curr = NULL;
    }
};

extern "C" {
    char* MNA_solver(int ncomp, char* data);
}

char* MNA_solver(int ncomp, char* data)
{
    auto compsline = split(data, "/");
    int i, n = 0, rn = 0;
    vector<vector<string>> comps;
    for(i = 0; i < ncomp; i++)
    {
        comps.push_back(split(compsline[i], " "));
        if(compsline[i][0]=='V' || compsline[i][0] == 'E')
            n++;
        else if(compsline[i][0]=='H')
            n += 2;
        rn = max(rn, max(atoi(comps[i][1].c_str()), atoi(comps[i][2].c_str())));
    }
    n += rn;
    Circuit ckt(n, rn);
    for(i = 0; i < ncomp; i++)
    {
        if(comps[i][0][0] == 'R')
            ckt.addResistor(stoi(comps[i][1]), stoi(comps[i][2]), stod(comps[i][3]));
        else if(comps[i][0][0] == 'I')
            ckt.addCurrentSource(stoi(comps[i][1]), stoi(comps[i][2]), stod(comps[i][3]));
        else if(comps[i][0][0] == 'V')
            ckt.addVoltageSource(stoi(comps[i][1]), stoi(comps[i][2]), stod(comps[i][3]));
        else if(comps[i][0][0] == 'G')
            ckt.addVCCS(stoi(comps[i][1]), stoi(comps[i][2]), stoi(comps[i][3]), stoi(comps[i][4]), stod(comps[i][5]));
        else if(comps[i][0][0] == 'E')
            ckt.addVCVS(stoi(comps[i][1]), stoi(comps[i][2]), stoi(comps[i][3]), stoi(comps[i][4]), stod(comps[i][5]));
        else if(comps[i][0][0] == 'F')
            ckt.addCCCS(stoi(comps[i][1]), stoi(comps[i][2]), stoi(comps[i][3]), stoi(comps[i][4]), stod(comps[i][5]));
        else if(comps[i][0][0] == 'H')
            ckt.addCCVS(stoi(comps[i][1]), stoi(comps[i][2]), stoi(comps[i][3]), stoi(comps[i][4]), stod(comps[i][5]));
    }
    double result[n];
    ckt.compute(result);
    string result_str = "";
    for(i = 0; i < n; i++)
        result_str += to_string(result[i]) + " ";
    return &result_str[0];
}