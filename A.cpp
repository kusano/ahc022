#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
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

        map<vector<int>, vector<int>> M;
        for (int i=0; i<N; i++)
            M[EI[i]] = vector<int>(2);

        for (int i=0; i<N; i++)
        {
            int tx = (X[i]+dx+L)%L;
            int ty = (Y[i]+dy+L)%L;
            if (PP[ty][tx]!=-1)
                M[EI[i]][PP[ty][tx]]++;
        }

        for (int i=0; i<N; i++)
        {
            int tx = (X[i]+dx+L)%L;
            int ty = (Y[i]+dy+L)%L;
            if (PP[ty][tx]==-1)
            {
                int b;
                if (M[EI[i]][0]<M[EI[i]][1])
                    b = 0;
                else
                    b = 1;
                PP[ty][tx] = b;
                M[EI[i]][b]++;
            }
            EI[i].push_back(PP[ty][tx]);
        }
        if (set<vector<int>>(EI.begin(), EI.end()).size()==N)
            break;
    }
    cerr<<"BW: "<<EI[0].size()<<endl;

    vector<vector<int>> P(L, vector<int>(L, 500));
    for (int y=0; y<L; y++)
        for (int x=0; x<L; x++)
        {
            if (PP[y][x]==-1)
                P[y][x] = 500;
            else
                if (PP[y][x]==0)
                    P[y][x] = max(0, 500-2*S);
                else
                    P[y][x] = min(1000, 500+2*S);
        }

    for (int i=0; i<100; i++)
    {
        vector<vector<int>> T(L, vector<int>(L));
        T.swap(P);

        for (int y=0; y<L; y++)
            for (int x=0; x<L; x++)
                if (PP[y][x]==-1)
                    P[y][x] = (T[(y+L-1)%L][x]+T[(y+1)%L][x]+T[y][(x+L-1)%L]+T[y][(x+1)%L]+2)/4;
                else
                    P[y][x] = T[y][x];
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
        int n = 10000/N/(int)DX.size();

        vector<int> B;
        for (int j=0; j<(int)DX.size(); j++)
        {
            vector<int> T;
            for (int k=0; k<n; k++)
            {
                cout<<i<<" "<<DY[j]<<" "<<DX[j]<<endl;
                int t;
                cin>>t;
                T.push_back(t);
            }
            sort(T.begin(), T.end());
            int t;
            if (n%2==0)
                t = (T[n/2-1]+T[n/2]);
            else
                t = T[n/2];
            B.push_back(t>500?1:0);
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
