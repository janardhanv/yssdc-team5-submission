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

const int kSideBlocksOverlap = 1024;

const int kHeartBeat = 100000;
const int kMaxLen = 16<<20;

int n;
char s[kMaxLen+1];
const char* nucle = "ATGC";
int a[kGeneMaxLen][kGeneMaxLen][2];
int prv[kGeneMaxLen][kGeneMaxLen][2];
int fhlx[kGeneMaxLen];
int nhlx[kGeneMaxLen][kGeneMaxLen];
struct State* st[kGeneMaxLen][kGeneMaxLen];
int fwd[kMaxLen];
int back[kMaxLen];
int fatSums[kMaxLen];
int offset;

#define A(x) a[(x)%kGeneMaxLen]
#define P(x) prv[(x)%kGeneMaxLen]
#define ST(x) st[(x)%kGeneMaxLen]
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
	int getJStart() { return start+len-1+offset; }
	int getJEnd() { return end-len+1+offset; }
};

int creations = 0;

struct State
{
	Helix helix;
	State *left,*right;
	int value;
	bool sel;
	State(Helix h, int v, State *l=0, State *r=0): helix(h),value(v),left(l),right(r),sel(false) {++creations;}
	void collectHelixes(vector<Helix>& res) const
	{
		if (helix.len > 0)
			res.pb(helix);
		if (left) left->collectHelixes(res);
		if (right) right->collectHelixes(res);
	}
};

State* retrievePath(int end, int len)
{
	if (len <= 0 || A(end)[len][0] == 0)
		return NULL;
	State* &res = ST(end)[len];
	if (res != NULL)
		return res;
	int p = P(end)[len][0];
	if (p < 0)
	{
		p=-p;
		res = new State(Helix(end-len+1,end,p),A(end)[len][0],retrievePath(end-p,len-2*p));
	}
	else
	{
		State* s1 = retrievePath(end-p,len-p);
		State* s2 = retrievePath(end,p);
		if (s2 == 0)
			res = s1;
		else if (s1 == 0)
			res = s2;
		else
			res = new State(Helix(end-len+1,end,0),A(end)[len][0],s1,s2);
	}
	return res;
}

vector<State*> genes;

void solve()
{
	REP(end,n)
	{
		if (end % kHeartBeat == 0)
			LOG( << "Heartbeat " << end << "/" << n << "\t" << SZ(genes) << " c-genes so far");
		FHLX(end) = -1;
		CLEAR(A(end));
		CLEAR(ST(end));
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
			if (len >= kGeneMinLen && A(end)[len][0]*2*100 > kPairedMin100*len)
				if (!isFat(end-len+1, end))
				{
					// Candidate gene
					genes.pb(retrievePath(end,len));
					//fprintf(stderr,"Gene end %d len %d, score %.6lf\n",end,len,A(end)[len][0]*2.0/len);
				}
			if (A(end)[len-1][0] > A(end)[len][0])
			{
				A(end)[len][0] = A(end)[len-1][0];
				P(end)[len][0] = len-1;
				ST(end)[len] = NULL;
			}
			if (A(end-1)[len-1][0] > A(end)[len][0])
			{
				A(end)[len][0] = A(end-1)[len-1][0];
				P(end)[len][0] = 1;
				ST(end)[len] = NULL;
			}
		}
	}
}

void write(const State* g)
{
	vector<Helix> h;
	g->collectHelixes(h);
	printf("%d\n",SZ(h));
	REP(i,SZ(h))
		printf("%d %d %d\n",h[i].getJStart(),h[i].getJEnd(),h[i].len);
}

string serialize(const State* g)
{
	vector<Helix> h;
	g->collectHelixes(h);
	stringstream str;
	REP(i,SZ(h))
	{
		if (i) str << '#';
		str << h[i].getJStart() << " " << h[i].getJEnd() << " " << h[i].len;
	}
	return str.str();
}

void writeAsMapper()
{
	REP(i,SZ(genes))
		if (genes[i]->sel)
			printf("cgene\t%s\n",serialize(genes[i]).c_str());
}

void writeAsAnswer()
{
	int k = 0;
	REP(i,SZ(genes))
		k += genes[i]->sel;
	printf("%d\n",SZ(genes));
	REP(i,SZ(genes))
		if (genes[i]->sel)
			write(genes[i]);
}

void processGbk()
{
	readGbk();
	//n=1000000;
	prepareCodes();
	solve();
	//printf("ans = %d\n",A(n-1)[n][0]);
	LOG( << "total c-genes " << SZ(genes));
	//LOG( << creations << " creations");
	selectGenes();
	writeAsMapper();
}

void processSplitted()
{
	while (readSplitted())
	{
		prepareCodes();
		solve();
		LOG( << "total c-genes " << SZ(genes));
		LOG( << "========================");
		selectGenes();
		writeAsMapper();
	}
}

bool cmpGeneEnd(const State* x, const State* y)
{
	return x->helix.end < y->helix.end;
}

void selectGenes()
{
	static int sum[kMaxLen];
	static State* prev[kMaxLen];
	sort(ALL(genes),cmpGeneEnd);
	int p = 0;
	sum[0] = 0;
	prev[0] = NULL;
	FOR(i,1,n-1)
	{
		if (i)
			sum[i] = sum[i-1];
		for (; p < SZ(genes) && genes[p]->helix.end == i; p++)
		{
			int t = genes[p]->value;
			if (genes[p]->helix.start > 0)
				t += sum[genes[p]->helix.start-1];
			if (t > sum[i])
			{
				sum[i] = t;
				prev[i] = genes[p];
			}
		}
	}
	for (int i = n-1; i >= 0; i = prev[i] ? prev[i]->helix.start-1 : i-1)
	{
		if (prev[i])
			prev[i]->sel = true;
	}
	REP(i,SZ(genes))
		if (genes[i]->helix.start < kSideBlocksOverlap || genes[i]->helix.end >= n-kSideBlocksOverlap)
			genes[i]->sel = true;
}

int main()
{
	//freopen("ref_chr7_00.gbk","r",stdin);
	//freopen("split.txt","r",stdin);
	//freopen("results.out","w",stdout);

	//processGbk();
	processSplitted();

	timer.debug("Done!");
	return 0;
}
