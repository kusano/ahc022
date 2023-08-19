#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <cmath>
#include <utility>
using namespace std;

int xor64() {
    static uint64_t x = 88172645463345263ULL;
    x ^= x<<13;
    x ^= x>> 7;
    x ^= x<<17;
    return int(x&0x7fffffff);
}

class Comm
{
public:
    virtual void get_param(int *L, int *N, int *S, vector<int> *X, vector<int> *Y) = 0;
    virtual void place(vector<vector<int>> P) = 0;
    virtual int measure(int i, int x, int y) = 0;
    virtual long long answer(vector<int> E) = 0;
};

// 提出用。
class CommJudge: public Comm
{
public:
    void get_param(int *L, int *N, int *S, vector<int> *X, vector<int> *Y)
    {
        cin>>*L>>*N>>*S;
        X->resize(*N);
        Y->resize(*N);
        for (int i=0; i<*N; i++)
            cin>>(*Y)[i]>>(*X)[i];
    }

    void place(vector<vector<int>> P)
    {
        int L = (int)P.size();
        for (int y=0; y<L; y++)
        {
            for (int x=0; x<L; x++)
                cout<<(x==0?"":" ")<<P[y][x];
            cout<<endl;
        }
    }

    int measure(int i, int x, int y)
    {
        cout<<i<<" "<<y<<" "<<x<<endl;
        int m;
        cin>>m;
        return m;
    }

    long long answer(vector<int> E)
    {
        int N = (int)E.size();
        cout<<-1<<" "<<-1<<" "<<-1<<endl;
        for (int e: E)
            cout<<e<<endl;

        return 0;
    }
};

// ファイルから読み込み。
class CommFile: public Comm
{
    int L, N, S;
    vector<int> X, Y;
    vector<int> A;
    vector<int> f;
    vector<vector<int>> P;
    long long place_cost = 0;
    long long measure_cost = 0;
    int measure_count = 0;

public:
    CommFile()
    {
        cin>>L>>N>>S;
        X.resize(N);
        Y.resize(N);
        for (int i=0; i<N; i++)
            cin>>Y[i]>>X[i];
        A.resize(N);
        for (int &a: A)
            cin>>a;
        f.resize(10000);
        for (int &x: f)
            cin>>x;
    }

    void get_param(int *L, int *N, int *S, vector<int> *X, vector<int> *Y)
    {
        *L = this->L;
        *N = this->N;
        *S = this->S;
        *X = this->X;
        *Y = this->Y;
    }

    void place(vector<vector<int>> P)
    {
        this->P = P;

        for (int y=0; y<L; y++)
        {
            for (int x=0; x<L; x++)
                cout<<(x==0?"":" ")<<P[y][x];
            cout<<endl;
        }

        place_cost = 0;
        for (int y=0; y<L; y++)
            for (int x=0; x<L; x++)
            {
                int p = P[y][x];
                int px = P[y][(x+1)%L];
                int py = P[(y+1)%L][x];
                place_cost += (p-px)*(p-px)+(p-py)*(p-py);
            }
    }

    int measure(int i, int x, int y)
    {
        cout<<i<<" "<<y<<" "<<x<<endl;
        int m = P[((Y[A[i]]+y)%L+L)%L][((X[A[i]]+x)%L+L)%L] + f[measure_count++];
        m = max(0, min(1000, m));

        measure_cost += 100*(10+abs(x)+abs(y));

        return m;
    }

    long long answer(vector<int> E)
    {
        cout<<-1<<" "<<-1<<" "<<-1<<endl;
        for (int e: E)
            cout<<e<<endl;

        int W = 0;
        for (int i=0; i<N; i++)
            if (A[i]!=E[i])
                W++;

        long long score = (long long)ceil(1e14*pow(.8, W)/(place_cost+measure_cost+1e5));

        fprintf(
            stderr,
            "%2d %3d %3d %12lld %3d %12lld %12lld %5d\n",
            L,
            N,
            S,
            score,
            W,
            place_cost,
            measure_cost,
            measure_count);

        return score;
    }
};

// パラメタを与えて生成。
class CommParam: public Comm
{
    int L, N, S;
    vector<int> X, Y;
    vector<int> A;
    vector<int> f;
    vector<vector<int>> P;
    long long place_cost = 0;
    long long measure_cost = 0;
    int measure_count = 0;

    int xor64() {
        static uint64_t x = 88172645463345263ULL;
        x ^= x<<13;
        x ^= x>> 7;
        x ^= x<<17;
        return int(x&0x7fffffff);
    }

public:
    CommParam(int L, int N, int S)
    {
        if (L==-1)
            L = xor64()%41+10;
        if (N==-1)
            N = xor64()%41+60;
        if (S==-1)
        {
            int s = xor64()%30+1;
            S = s*s;
        }
        this->L = L;
        this->N = N;
        this->S = S;

        set<pair<int, int>> XY;

        for (int i=0; i<N; i++)
        {
            while (true)
            {
                int x = xor64()%L;
                int y = xor64()%L;
                if (XY.count({y, x})==0)
                {
                    XY.insert({y, x});
                    break;
                }
            }
        }
        for (auto xy: XY)
        {
            X.push_back(xy.second);
            Y.push_back(xy.first);
        }

        for (int i=0; i<N; i++)
            A.push_back(i);
        for (int i=0; i<N; i++)
            swap(A[i], A[xor64()%(N-i)+i]);

        double pi = acos(-1);
        for (int i=0; i<10000; i++)
        {
            double x = (double)(xor64()+1)/0x80000001;
            double y = (double)(xor64()+1)/0x80000001;
            f.push_back((int)round(sqrt(-2*log(x))*cos(2*pi*y)*S));
        }
    }

    void get_param(int *L, int *N, int *S, vector<int> *X, vector<int> *Y)
    {
        *L = this->L;
        *N = this->N;
        *S = this->S;
        *X = this->X;
        *Y = this->Y;
    }

    void place(vector<vector<int>> P)
    {
        this->P = P;

        place_cost = 0;
        for (int y=0; y<L; y++)
            for (int x=0; x<L; x++)
            {
                int p = P[y][x];
                int px = P[y][(x+1)%L];
                int py = P[(y+1)%L][x];
                place_cost += (p-px)*(p-px)+(p-py)*(p-py);
            }
    }

    int measure(int i, int x, int y)
    {
        int m = P[((Y[A[i]]+y)%L+L)%L][((X[A[i]]+x)%L+L)%L] + f[measure_count++];
        m = max(0, min(1000, m));

        measure_cost += 100*(10+abs(x)+abs(y));

        return m;
    }

    long long answer(vector<int> E)
    {
        int W = 0;
        for (int i=0; i<N; i++)
            if (A[i]!=E[i])
                W++;

        long long score = (long long)ceil(1e14*pow(.8, W)/(place_cost+measure_cost+1e5));

        fprintf(
            stderr,
            "%2d %3d %3d %12lld %3d %12lld %12lld %5d\n",
            L,
            N,
            S,
            score,
            W,
            place_cost,
            measure_cost,
            measure_count);

        return score;
    }
};

// R は C のうちどれから発生したものかを予測する。
int predict(int S, vector<int> C, vector<int> R)
{
    static int initS = -1;
    static vector<double> T;
    static double rest;

    if (initS!=S)
    {
        double pi = acos(-1.0);

        T.resize(2001);
        for (int i=-1000; i<=1000; i++)
            T[i+1000] = 1./(sqrt(2*pi)*S)*exp(-(double)i*i/(2*S*S));

        rest = 1.0;
        for (double t: T)
            rest -= t;
    }

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
                for (int j=r-c+1000; j<=2000; j++)
                    p += T[j];
                p += rest/2.;
            }
            else
                p = T[r-c+1000];

            if (p<1e-100)
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

void solve(Comm *comm)
{
    int L, N, S;
    vector<int> X, Y;
    comm->get_param(&L, &N, &S, &X, &Y);

    vector<vector<int>> PP(L, vector<int>(L, -1));
    vector<int> DX, DY;
    vector<vector<int>> EI(N);
    //  階調の分割数
    int W = 10;
    if (S>=100)
        W = 5;
    if (S>=144)
        W = 4;
    if (S>=196)
        W = 3;
    if (S>=256)
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
        if (set<vector<int>>(EI.begin(), EI.end()).size()==N || DX.size()>=16)
            break;
    }
    //cerr<<"BW: "<<EI[0].size()<<endl;

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

    comm->place(P);

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
                int t = comm->measure(i, DX[j], DY[j]);
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

    comm->answer(E);
}

int main(int argc, char **argv)
{
    /*
    for (int i=0; i<30; i++)
    {
        for (int j=0; j<32; j++)
        {
            CommParam comm(-1, -1, (i+1)*(i+1));
            solve(&comm);
        }
    }
    return 0;
    //*/

    Comm *comm;

    if (argc==1)
        comm = new CommJudge;
    else
        comm = new CommFile;

    solve(comm);
}
