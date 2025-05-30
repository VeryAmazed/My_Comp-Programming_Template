#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
// if you end up using long double, you need to set the floating point notation to fixed, and set the percision to be very high
typedef long double ld;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;
typedef vector<int> vi;

// contrsuct umaps like this, unordered_map<long long, int, custom_hash> safe_map;
// FIXED_RANDOM is static so it doesn not get redeclared between function calls
struct custom_hash {
    static uint64_t splitmix64(uint64_t x) {
        // http://xorshift.di.unimi.it/splitmix64.c
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return x ^ (x >> 31);
    }

    size_t operator()(uint64_t x) const {
		
        static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x + FIXED_RANDOM);
    }
};


const ll INF = 2001001001;
const ll INF2 = 2e18;
const ll MOD = 1000000007;

#define max3(a, b, c) max(a, max(b, c))
#define min3(a, b, c) min(a, min(b, c))
#define pb push_back 
#define pf push_front
#define f first
#define s second
#define mp make_pair
#define sz(x) (int)(x).size()
#define all(x) begin(x), end(x)

const bool testing = 0;
mt19937 rng(63);

// first four are north, west, east ,south
int dir1[] = {1, 0, -1, 0, 1, 1, -1, -1};
int dir2[] = {0, 1, 0, -1, 1, -1, 1, -1};

int main() {
	// apparently this does fast i/o
	cin.tie(0) , ios::sync_with_stdio(0);
	cin.exceptions(cin.failbit);
	// use this if you read in from a file
	/*
	freopen("in.txt", "r", stdin);
    freopen("out.txt", "w", stdout);
	*/
	
	stringstream ss;
	
	// Do it once. Do it right.
	// Read the problem statement carefully, if it's a physical copy, underline key points
	// Plan out the steps in words on a piece of paper before implementing
	// after RTE(obviously) but also WA, run valgrind!!!
	// Testing your solution on samples before coding is a great way to see if you read the problem correctly!!!
	// Also take notes about key elements in the problem statement while reading the problem!!!
	// If you're stuck, try small cases (not just a trick for math problems)
	// When debugging, especially when it's not your code, go through the code line by line, don't just scan
	// over sections, and if it's not your code, take notes on what each line does.
	
	//cout << fixed << setprecision(12);
	// if you use ld, use the above and don't use string stream
	
	// use instead of ceil(a, b) if a and b are positive
	// (a + b - 1) / b
	cout << ss.str();
	return 0;
}
