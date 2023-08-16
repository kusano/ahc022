#include <iostream>
#include <vector>
using namespace std;

int main()
{
    int L, N, S;
    cin>>L>>N>>S;
    vector<int> X(N), Y(N);
    for (int i=0; i<N; i++)
        cin>>Y[i]>>X[i];

    vector<vector<int>> P(L, vector<int>(L));
    for (int i=0; i<N; i++)
        P[Y[i]][X[i]] = i*10;

    for (int y=0; y<L; y++)
    {
        for (int x=0; x<L; x++)
            cout<<(x==0?"":" ")<<P[y][x];
        cout<<endl;
    }

    vector<int> E(N);
    for (int i=0; i<N; i++)
    {
        cout<<i<<" "<<0<<" "<<0<<endl;
        int t;
        cin>>t;
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
