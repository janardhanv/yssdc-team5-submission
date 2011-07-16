#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <bitset>
#include <set>
#include <sstream>
#include <stdlib.h>
#include <map>
#include <queue>


using namespace std;


#define sz(x) ((int)x.size())
#define all(x) (x).begin(), (x).end()
#define pb(x) push_back(x)
#define mp(x, y) make_pair(x, y)

typedef long long int64;

typedef vector <int> vi;
typedef vector < vi > vvi;


const int kBlockSize = 1e5;
const int kOverlapSize = 1e3;


int main () {
	//freopen("", "rt", stdin);
	//freopen("", "wt", stdout);

    string s;

    getline(cin, s);
    getline(cin, s);

    for (int pos = 0; pos < sz(s); pos += kBlockSize - kOverlapSize) {
        cout << pos << "\t" << s.substr(pos, kBlockSize) << endl;
    }



    return 0;
}

