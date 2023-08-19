#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <cmath>
using namespace std;

int xor64() {
    static uint64_t x = 88172645463345263ULL;
    x ^= x<<13;
    x ^= x>> 7;
    x ^= x<<17;
    return int(x&0x7fffffff);
}

int predict(int S, vector<int> C, vector<int> R)
{
    double pi = acos(-1.0);

    vector<double> T(2001);
    for (int i=-1000; i<=1000; i++)
        T[i+1000] = 1./(sqrt(2*pi)*S)*exp(-(double)i*i/(2*S*S));

    double rest = 1.0;
    for (double t: T)
        rest -= t;

    vector<double> P(C.size());
    for (int r: R)
    {
        for (int i=0; i<(int)C.size(); i++)
        {
            int c = C[i];
            double p = 0.;
            if (r==0)
            {
                p += rest/2.;
                for (int j=0; j<=r-c+1000; j++)
                    p += T[j];
            }
            else if (r==1000)
            {
                for (int j=r-c+1000; j<=2001; j++)
                    p += T[j];
                p += rest/2.;
            }
            else
                p = T[r-c+1000];

            if (p==0.)
                P[i] += -1e10;
            else
                P[i] += log(p);
        }
    }

    int r = 0;
    for (int i=1; i<(int)P.size(); i++)
        if (P[i]>P[r])
            r = i;
    return r;
}

int main()
{
    //cout<<predict(4, {0, 250, 500, 750, 1000}, {251, 250, 240})<<endl;
    //return 0;

    int L, N, S;
    cin>>L>>N>>S;
    cerr<<L<<" "<<N<<" "<<S<<endl;
    vector<int> X(N), Y(N);
    for (int i=0; i<N; i++)
        cin>>Y[i]>>X[i];

    vector<vector<int>> PP(L, vector<int>(L, -1));
    vector<int> DX, DY;
    vector<vector<int>> EI(N);
    //  階調の分割数
    int W = 5;
    if (S>=324)
        W = 3;
    if (S>=529)
        W = 2;
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
            M[EI[i]] = vector<int>(W);

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
                int b = 0;
                for (int j=0; j<W; j++)
                    if (M[EI[i]][j]<M[EI[i]][b])
                        b = j;
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
            int w = 1000/(W-1);
            if (PP[y][x]==-1)
                P[y][x] = 500;
            else
                P[y][x] = PP[y][x]*w;
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

            vector<int> cand;
            int w = 1000/(W-1);
            for (int k=0; k<W; k++)
                cand.push_back(k*w);
            int b = predict(S, cand, T);
            B.push_back(b);
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
