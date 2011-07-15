#define _CRT_SECURE_NO_WARNINGS
#include <algorithm>
#include <numeric>
#include <string>
#include <cstring>
#include <set>
#include <map>
#include <vector>
#include <queue>
#include <iostream>
#include <iterator>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>

using namespace std;

#define REP(i,n) for (int i=0,_n=(n); i < _n; i++)
#define REPD(i,n) for (int i=(n)-1; i >= 0; i--)
#define FOR(i,a,b) for (int _b=(b), i=(a); i <= _b; i++)
#define FORD(i,a,b) for(int i=(a),_b=(b);i>=_b;i--)
#define ALL(c) (c).begin(), (c).end()
#define SORT(c) sort(ALL(c))

#define CLEAR(x) memset(x,0,sizeof x);
#define CLEARA(x) memset(&x,0,sizeof x);
#define FILL(x,v) memset(x,v,sizeof x);
#define FILLA(x,v) memset(&x,v,sizeof x);

#define VAR(a,b) __typeof(b) a=(b)
#define FOREACH(it,c) for(VAR(it,(c).begin());it!=(c).end();++it)

#define REVERSE(c) reverse(ALL(c))
#define UNIQUE(c) SORT(c),(c).resize(unique(ALL(c))-(c).begin())
#define INF 0x7fffffff
#define X first
#define Y second
#define pb push_back
#define SZ(c) (int)(c).size()
#define MP make_pair
#define eps 1.0e-11
const double pi = acos(-1.0);

typedef pair<int, int> PII;
typedef vector<PII> VPII;
typedef vector<int> VI;
typedef vector<VI> VVI;
typedef long long LL;

const int kGeneMaxLen = 1024;
const int kHairpinDistMin = 4;
const int kHairpinDistMax = 20;
const int kHelixMin = 4;
const int kGeneMinLen = 10;
const double kPairedMin = 0.6;
const double kFatMin = 0.35;
const double kFatMax = 0.65;

const int kShiftLeft = 1000000000;
const int kShiftRight = 1000000001;

int n;
char* s;
const char* nucle = "ATGC";
int a[kGeneMaxLen][kGeneMaxLen][2];
int prv[kGeneMaxLen][kGeneMaxLen][2];
int fhlx[kGeneMaxLen][kGeneMaxLen];
int nhlx[kGeneMaxLen][kGeneMaxLen];
int phlx[kGeneMaxLen][kGeneMaxLen];
int* fwd;
int* back;
int* fatSums;

#define A(x) a[(x)%kGeneMaxLen]
#define P(x) prv[(x)%kGeneMaxLen]
#define FHLX(x) fhlx[(x)%kGeneMaxLen]
#define NHLX(x) nhlx[(x)%kGeneMaxLen]
#define PHLX(x) phlx[(x)%kGeneMaxLen]

void readData()
{
	scanf("%d",&n);
	s = new char[n+1];
	scanf("%s",s);
	REP(i,n)
	{
		s[i] = toupper(s[i]);
		if (s[i] == 'N')
			s[i] = -1;
		else
			s[i] = strchr(nucle,s[i])-nucle;
	}
}

void prepareCodes()
{
	fwd = new int[n];
	back = new int[n];
	REP(i,n)
	{
		fwd[i] = 0;
		REP(j,kHelixMin)
		{
			if (i+j >= n || s[i+j] == -1)
			{
				fwd[i] = -1;
				break;
			}
			fwd[i] = fwd[i]*4 + s[i+j];
		}
		back[i] = 0;
		REP(j,kHelixMin)
		{
			if (i-j < 0 || s[i-j] == -1)
			{
				back[i] = -1;
				break;
			}
			back[i] = back[i]*4 + (s[i-j]^1);
		}
	}
	fatSums = new int[n+1];
	fatSums[0]=0;
	REP(i,n)
	{
		fatSums[i+1] = (0 <= s[i] && s[i] < 2) + fatSums[i];
	}
}

bool isFat(int s, int e)
{
	double fat = (fatSums[e+1] - fatSums[s])/(double)(e-s+1);
	return kFatMin > fat || fat > kFatMax;
}

struct Helix
{
	int start,end,len;
	Helix(int s, int e, int l): start(s),end(e),len(l) {}
};

typedef vector<Helix> Gene;

void collect(int end, int len, vector<Helix> &res)
{
	if (len <= 0 || A(end)[len][0] == 0)
		return;
	int p = P(end)[len][0];
	if (p < 0)
	{
		p=-p;
		res.pb(Helix(end-len+p,end-p+1,p));
		collect(end-p,len-2*p,res);
	}
	else
	{
		collect(end-p,len-p,res);
		collect(end,p,res);
	}
}

vector<Gene> genes;

void solve()
{
	REP(end,n)
	{
		CLEAR(A(end));
		FOR(len,kHelixMin*2 + kHairpinDistMin, min(kGeneMaxLen-1, end+1))
		{
			int start = end-len+1;
			if (fwd[start] != -1 && fwd[start] == back[end]) // helix
			{
				A(end)[len][1] = A(end-kHelixMin)[len-2*kHelixMin][0] + kHelixMin;
				P(end)[len][1] = -kHelixMin;
				if (A(end-1)[len-2][1]+1 >= A(end)[len][1])
				{
					A(end)[len][1] = A(end-1)[len-2][1]+1;
					P(end)[len][1] = P(end-1)[len-2][1]-1;
				}
				A(end)[len][0] = A(end)[len][1];
				P(end)[len][0] = P(end)[len][1];
			}
			// TODO: optimize
			for (int cut = 1; cut < len; ++cut)
			{
				int t = A(end)[cut][0] + A(end-cut)[len-cut][0];
				if (t > A(end)[len][0])
				{
					A(end)[len][0] = t;
					P(end)[len][0] = cut;
				}
			}
			if (len >= kGeneMinLen && A(end)[len][1]*2.0/len >= kPairedMin)
				if (!isFat(end-len+1, end))
				{
					// Candidate gene
					genes.pb(Gene());
					collect(end,len,genes.back());
					printf("Gene end %d len %d, score %.6lf\n",end,len,A(end)[len][1]*2.0/len);
				}
		}
	}
}

void write(Gene g)
{
	printf("%d\n",SZ(g));
	REP(i,SZ(g))
		printf("%d %d %d\n",g[i].start,g[i].end,g[i].len);
}

int main()
{
	freopen("ref_chr7_00.gbk","r",stdin);
	readData();
	n=2000;
	prepareCodes();
	solve();
	printf("ans = %d\n",A(n-1)[n][0]);
	//REP(i,SZ(genes)) write(genes[i]);
	return 0;
}