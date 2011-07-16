#define _CRT_SECURE_NO_WARNINGS
#include <algorithm>
#include <cassert>
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
#include <ctime>

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

class Timer
{
	clock_t time_start;
public:
	Timer()
	{
		time_start = clock();
	}
	double elapsed()
	{
		return (clock()-time_start)/(double)CLOCKS_PER_SEC;
	}
	void debug(const string& s)
	{
		double e = elapsed();
		int min = (int)(e/60.0);
		fprintf(stderr,"%d:%08.5lf %s\n",min,e-min*60,s.data());
	}
};

Timer timer;

const int kGeneMaxLen = 1024;
const int kHairpinDistMin = 4;
const int kHairpinDistMax = 20;
const int kHelixMin = 4;
const int kGeneMinLen = 100;
const int kPairedMin100 = 60;
const int kFatMin100 = 35;
const int kFatMax100 = 65;
const double kPairedMin = kPairedMin100 * 0.01;
const double kFatMin = kFatMin100 * 0.01;
const double kFatMax = kFatMax100 * 0.01;

const int kHeartBeat = 100000;
const int kMaxLen = 16<<20;

int n;
char s[kMaxLen+1];
const char* nucle = "ATGC";
int a[kGeneMaxLen][kGeneMaxLen][2];
int prv[kGeneMaxLen][kGeneMaxLen][2];
int fhlx[kGeneMaxLen];
int nhlx[kGeneMaxLen][kGeneMaxLen];
int fwd[kMaxLen];
int back[kMaxLen];
int fatSums[kMaxLen];
int offset;

#define A(x) a[(x)%kGeneMaxLen]
#define P(x) prv[(x)%kGeneMaxLen]
#define FHLX(x) fhlx[(x)%kGeneMaxLen]
#define NHLX(x) nhlx[(x)%kGeneMaxLen]

#define LOGTO(t,x) { stringstream message; message x; t.debug(message.str()); }
#define LOG(x) LOGTO(timer,x)

void readGbk()
{
	scanf("%d",&n);
	assert(n < kMaxLen);
	scanf("%s",s);
	offset = 0;
	LOG( << "Read gbk file, length " << n);
}

bool readSplitted()
{
	if (scanf("%d",&offset) != 1)
		return false;
	assert(offset >= 0);
	scanf("%s",s);
	n = strlen(s);
	LOG( << "Read piece, offset " << offset << " length " << n);
	return true;
}

void prepareCodes()
{
	REP(i,n)
	{
		s[i] = toupper(s[i]);
		if (s[i] == 'N')
			s[i] = -1;
		else
			s[i] = strchr(nucle,s[i])-nucle;
	}
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
	fatSums[0]=0;
	REP(i,n)
	{
		fatSums[i+1] = (0 <= s[i] && s[i] < 2) + fatSums[i];
	}
}

inline bool isFat(int s, int e)
{
	int sum = fatSums[e+1] - fatSums[s];
	int len = e-s+1;
	return kFatMin100*len > sum*100 || sum*100 > kFatMax100*len;
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
		res.pb(Helix(end-len+p+offset,end-p+1+offset,p));
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
		if (end % kHeartBeat == 0)
			LOG( << "Heartbeat " << end << "/" << n << "\t" << SZ(genes) << " c-genes so far");
		FHLX(end) = -1;
		CLEAR(A(end));
		FOR(len,kHelixMin*2 + kHairpinDistMin, min(kGeneMaxLen-1, end+1))
		{
			int start = end-len+1;
			for (int w = FHLX(end); w != -1; w = NHLX(end)[w])
				if (A(end-w)[len-w][0] > 0 && P(end-w)[len-w][0] != (len-w-1))
				{
					int t = A(end)[w][0] + A(end-w)[len-w][0];
					if (t > A(end)[len][0])
					{
						A(end)[len][0] = t;
						P(end)[len][0] = w;
					}
				}
			if (fwd[start] != -1 && fwd[start] == back[end]) // helix
			{
				//printf("?? %d %d\n",start,end);
				A(end)[len][1] = A(end-kHelixMin)[len-2*kHelixMin][0];
				if (A(end)[len][1] > 0 || len-2*kHelixMin <= kHairpinDistMax)
				{
					A(end)[len][1] += kHelixMin;
					P(end)[len][1] = -kHelixMin;
				}
				if (A(end-1)[len-2][1] > 0 && A(end-1)[len-2][1]+1 >= A(end)[len][1])
				{
					A(end)[len][1] = A(end-1)[len-2][1]+1;
					P(end)[len][1] = P(end-1)[len-2][1]-1;
				}
				if (A(end)[len][1] > 0)
				{
					NHLX(end)[len] = FHLX(end);
					FHLX(end) = len;
					if (A(end)[len][1] > A(end)[len][0])
					{
						A(end)[len][0] = A(end)[len][1];
						P(end)[len][0] = P(end)[len][1];
					}
				}
			}
			if (len >= kGeneMinLen && A(end)[len][0]*2.0/len > kPairedMin)
				if (!isFat(end-len+1, end))
				{
					// Candidate gene
					genes.pb(Gene());
					collect(end,len,genes.back());
					//fprintf(stderr,"Gene end %d len %d, score %.6lf\n",end,len,A(end)[len][0]*2.0/len);
				}
			if (A(end)[len-1][0] > A(end)[len][0])
			{
				A(end)[len][0] = A(end)[len-1][0];
				P(end)[len][0] = len-1;
			}
			if (A(end-1)[len-1][0] > A(end)[len][0])
			{
				A(end)[len][0] = A(end-1)[len-1][0];
				P(end)[len][0] = 1;
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

string serialize(Gene g)
{
	stringstream str;
	REP(i,SZ(g))
	{
		if (i) str << '#';
		str << g[i].start << " " << g[i].end << " " << g[i].len;
	}
	return str.str();
}

void writeAsMapper()
{
	REP(i,SZ(genes))
		printf("%s\n",serialize(genes[i]));
}

void writeAsAnswer()
{
	printf("%d\n",SZ(genes));
	REP(i,SZ(genes))
		write(genes[i]);
}

void processGbk()
{
	readGbk();
	//n=100000;
	prepareCodes();
	solve();
	//printf("ans = %d\n",A(n-1)[n][0]);
}

void processSplitted()
{
	while (readSplitted())
	{
		prepareCodes();
		solve();
	}
}

int main()
{
	//freopen("ref_chr7_00.gbk","r",stdin);
	//freopen("split.txt","r",stdin);
	//freopen("results.out","w",stdout);

	//processGbk();
	processSplitted();
	LOG( << "total c-genes " << SZ(genes));
	writeAsMapper();

	timer.debug("Done!");
	return 0;
}