#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
using namespace std;

int xor64() {
    static uint64_t x = 88172645463345263ULL;
    x ^= x<<13;
    x ^= x>> 7;
    x ^= x<<17;
    return int(x&0x7fffffff);
}

int main()
{
    int L, N, S;
    cin>>L>>N>>S;
    cerr<<L<<" "<<N<<" "<<S<<endl;
    vector<int> X(N), Y(N);
    for (int i=0; i<N; i++)
        cin>>Y[i]>>X[i];

    vector<vector<int>> PP(L, vector<int>(L, -1));
    vector<int> DX, DY;
    vector<vector<int>> EI(N);
    while (true)
    {
        int dx, dy;
        while (true)
        {
            dx = xor64()%21-10;
            dy = xor64()%21-10;
            if (abs(dx)+abs(dy)<=10) {
                bool ok = true;
                for (int i=0; i<(int)DX.size(); i++)
                    if (dx==DX[i] && dy==DY[i])
                        ok = false;
                if (ok)
                    break;
            }
        }
        DX.push_back(dx);
        DY.push_back(dy);

        for (int i=0; i<N; i++)
        {
            int tx = (X[i]+dx+L)%L;
            int ty = (Y[i]+dy+L)%L;
            int b = -1;
            if (PP[ty][tx]==-1)
            {
                b = xor64()%2;
                PP[ty][tx] = b;
            }
            else
                b = PP[ty][tx];
            EI[i].push_back(b);
        }
        if (set<vector<int>>(EI.begin(), EI.end()).size()==N)
            break;
    }

    vector<vector<int>> P(L, vector<int>(L, 500));
    for (int y=0; y<L; y++)
        for (int x=0; x<L; x++)
        {
            if (PP[y][x]==-1)
                P[y][x] = 500;
            else
                P[y][x] = PP[y][x]*1000;
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
        vector<int> B;
        for (int j=0; j<(int)DX.size(); j++)
        {
            vector<int> T;
            for (int k=0; k<5; k++)
            {
                cout<<i<<" "<<DY[j]<<" "<<DX[j]<<endl;
                int t;
                cin>>t;
                T.push_back(t);
            }
            sort(T.begin(), T.end());
            B.push_back(T[2]>500?1:0);
        }
        int e = 0;
        for (int j=0; j<N; j++)
            if (B==EI[j])
                e = j;
        E[i] = e;
    }

    cout<<-1<<" "<<-1<<" "<<-1<<endl;
    for (int e: E)
        cout<<e<<endl;
}
