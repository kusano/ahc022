#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

int main()
{
    int L, N, S;
    cin>>L>>N>>S;
    cerr<<L<<" "<<N<<" "<<S<<endl;
    vector<int> X(N), Y(N);
    for (int i=0; i<N; i++)
        cin>>Y[i]>>X[i];

    vector<int> I;
    for (int i=0; i<N/2; i++)
    {
        I.push_back(i);
        I.push_back(N-i-1);
    }
    if (N%2!=0)
        I.push_back(N/2);

    vector<vector<int>> P(L, vector<int>(L, 500));
    vector<vector<bool>> PF(L, vector<bool>(L));
    for (int i=0; i<N; i++)
    {
        P[Y[I[i]]][X[I[i]]] = 1000*i/(N-1);
        PF[Y[I[i]]][X[I[i]]] = true;
    }

    for (int i=0; i<100; i++)
    {
        vector<vector<int>> T = P;
        for (int y=0; y<L; y++)
            for (int x=0; x<L; x++)
                if (PF[y][x])
                    P[y][x] = T[y][x];
                else
                {
                    int s = (T[(y+L-1)%L][x] + T[(y+1)%L][x] + T[y][(x+L-1)%L] + T[y][(x+1)%L]+2)/4;
                    P[y][x] = s;
                }
    }

    for (int y=0; y<L; y++)
    {
        for (int x=0; x<L; x++)
            cout<<(x==0?"":" ")<<P[y][x];
        cout<<endl;
    }

    vector<int> E(N);
    for (int i=0; i<N; i++)
    {
        vector<int> T;
        for (int j=0; j<99; j++)
        {
            cout<<i<<" "<<0<<" "<<0<<endl;
            int t;
            cin>>t;
            T.push_back(t);
        }
        sort(T.begin(), T.end());
        int t = T[49];

        int m = 9999;
        for (int j=0; j<N; j++)
        {
            int d = abs(P[Y[j]][X[j]]-t);
            if (d<m)
            {
                m = d;
                E[i] = j;
            }
        }
    }

    cout<<-1<<" "<<-1<<" "<<-1<<endl;
    for (int e: E)
        cout<<e<<endl;
}
