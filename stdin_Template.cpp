#include <vector>
#include <set>
#include <map>
#include <queue>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iterator>
#include <utility>
#include <deque>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <cstdio>
#include <cstring>
#include <stack>
#include <functional>
using namespace std;
typedef long long ll;
// if you end up using long double, you need to set the floating point notation to fixed, and set the percision to be very high
typedef long double ld;

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


#define INF 2001001001
#define INF2 2e18
#define MOD 1000000007

#define f0r(a, b) for (long long a = 0; a < b; a++)
#define f1r(a, b, c) for(long long a = b; a < c; a++)
#define max3(a, b, c) max(a, max(b, c))
#define min3(a, b, c) min(a, min(b, c))
#define pb push_back 
#define f first
#define s second
#define mp make_pair
#define pll pair<ll, ll>
#define pii pair<int, int>
#define tp make_tuple


int main() {
	// apparently this does fast i/o
	cin.tie(0) , ios::sync_with_stdio(0);
	
	// use this if you read in from a file
	/*
	freopen("in.txt", "r", stdin);
    freopen("out.txt", "w", stdout);
	*/
	
	stringstream ss;
	// Do it once. Do it right.
	// Plan out the steps in words on a piece of paper before implementing
	
	//cout << fixed << setprecision(12);
	// if you use ld, use the above and don't use string stream

	cout << ss.str();
	return 0;
}

///////////////////////////////////////////////////////////////////////////////

//template for compare functor
struct cmp {
	bool operator()(const Class& x, const Class& y) const { return x.a < y.a; }
};

// size of components are stored as negative values
// any node with a negative value is the representative node
struct DSU{
	vector<int> sets;
	void init(int n){ sets = vector<int>(n+1 ,-1);}
	// uses path compression
	int find(int x){
		return sets[x] < 0 ? x : sets[x] = find(sets[x]);
	}
	bool unite(int x, int y){ // union by size
		x = find(x);
		y = find(y);
		if(x == y) return 0;
		if(sets[x] > sets[y]) swap(x, y);
		sets[x] += sets[y];
		sets[y] = x;
		return 1;
	}
	
};

ll gcd(ll a, ll b) {
	if (b == 0) return a;
	return gcd(b, a % b);
}

// for mod exponentiation simply pass in a third parameter, m (mod), and then mod after every time you multiply
ll binpow(ll a, ll b) {
    ll res = 1;
    while (b > 0) {
        if (b & 1)
            res = res * a; // for mod exponentiation add a mod m here
        a = a * a; // for mod exponentiation add a mod m here
        b >>= 1;
    }
    return res;
}
// find mod inverse of a number. Uses fermat's little theorem so only works when mod is prime
ll mod_inv(ll a, ll mod){
	if(gcd(a, mod) == 1){
		return binpow(a , mod-2);
	}
	else{
		return -1; // a doesn't have a mod inverse in this mod space
	}
}

// Sieve of Eratosthenes
// prime numbers are set to 0 in sieve[], composite numbers are set to some multiple of a prime
ll sieve[n+1] = {0};
for (int x = 2; x <= n; x++) {
	if (sieve[x]) continue;
	for (int u = 2*x; u <= n; u += x) {
		sieve[u] = x;
	}
}

// Z-algorithm, generates an array z, where z[i] is the length of a substring substr, where substr is a prefix of the string str, and it starts at index i in str
// can be used be used to find in number of occurences of a string a, as a substring of another string b, simply input a + '$' + b, as the parameter, '$' is any character that will not be in string A or B
// doesn't need to take a string as input, because a string is just a vector of characters, which is just a vector is integers, you can also just input a vector of integers 
vector<int> z(string& s) {
	int n = s.size();
	vector<int> z(n);
	int x = 0, y = 0;
	for (int i = 1; i < n; i++) {
		z[i] = max(0,min(z[i-x],y-i+1));
		while (i+z[i] < n && s[z[i]] == s[i+z[i]]) {
			x = i; y = i+z[i]; z[i]++;
		}
	}
	return z;
}

// String Hashing Stuff
// All the precalculation for 2 base stringhash
// to get the hash of a string: 
// 		((poly1[end+1]-poly1[start]*pow1[end-start+1])%M1 + M1)%M1; 
ll M1 = 1e+7;
ll M2 = 1e+9;
const ll P1 = (uint64_t)chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count() % (uint64_t)M1;
const ll P2 = (uint64_t)chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count() % (uint64_t)M2;

vector<ll> pow1(str.size()+1, 0);
vector<ll> pow2(str.size()+1, 0);
vector<ll> poly1(str.size()+1, 0);
vector<ll> poly2(str.size()+1, 0);
pow1[0] =1;
pow2[0] =1;
for(int i =0; i < n; i++){
	pow1[i+1] = (pow1[i]*P1)%M1;
	pow2[i+1] = (pow2[i]*P2)%M2;
	poly1[i+1] = ((poly1[i]*P1)%M1 + str[i])%M1;
	poly2[i+1] = ((poly2[i]*P2)%M2 + str[i])%M2;
}
	
// monkey8's Lazy Prop Seg Tree Template
// Modified by Me
// Modified again to to allow for different updates and query types

// p is curent pos in array, always start at 1
// L, R is current bounds of the segment, always start at 0, size-1
// i, j is the range you are querying, both inclusive
// default version of this tempalte is for adding
struct Lazy_Seg_T{
	// modify the default value of the lazy array based on what is needed
	const ll LAZY_DEFAULT = 0;
	// modify the default value of lazy_set based on what is needed
	// make sure this is out of range of the set of posible values
	const ll SET_DEFAULT = 0;
	// modify the default value of the seg tree
	const ll INIT_DEFAULT = 0;
	vector<ll> st, lazy, lazy_set;
	int size;
	
	Lazy_Seg_T(int n, vector<ll>& in){
		size = n;
		while((size&(size-1))){
			size++;
		}
		
		st.assign(2*size, INIT_DEFAULT);
		lazy.assign(2*size, LAZY_DEFAULT);
		lazy_set.assign(2*size, SET_DEFAULT);
		for(int i = 0; i < in.size(); i++){
			int pos = i+size;
			st[pos] = in[i];
			for(pos/=2; pos >= 1; pos/=2){
				//change to query comb
				st[pos] = query_comb(st[left(pos)], st[right(pos)]);
			}
		}
	}
	
	Lazy_Seg_T(int n){
		size = n;
		while((size&(size-1))){
			size++;
		}
		st.assign(2*size, INIT_DEFAULT);
		lazy.assign(2*size, LAZY_DEFAULT);
		lazy_set.assign(2*size, SET_DEFAULT);
		
	}
	
	// change depending on how you want to combine for queries (i.e. min query, sum)
	// Template ver is for summing
	ll query_comb(ll a, ll b){
		return a + b;
	}
	// change depending on how you want to combine for updates (i.e. add to a range)
	// template is for summing
	ll update_comb(ll a, ll b){
		return a + b;
	}
	
	// helper
	int left(int p) {
		return (p << 1);
	}
	
	// helper
	int right(int p) {
		return (p << 1) + 1;
	}
	// Based on Sol to CSES Range Updates and Sums
	/* https://usaco.guide/plat/RURQ?lang=cpp */
	void push(int p, int L, int R) {
		
		if(lazy_set[p] != SET_DEFAULT){
			// change depending on what your query combine is
			// i.e. for min query, replace with st[p] = lazy_set[p] since now all the elements in this range are the same
			// default is for summing
			st[p] = lazy_set[p] * (R - L + 1);
			
			lazy[p] = LAZY_DEFAULT;
			if(L != R) {
				lazy[left(p)] = LAZY_DEFAULT;
				lazy[right(p)] = LAZY_DEFAULT;
				lazy_set[left(p)] = lazy_set[p];
				lazy_set[right(p)] = lazy_set[p];
				
			}
			lazy_set[p] = SET_DEFAULT;
		}
		else if(lazy[p] != LAZY_DEFAULT){
			// change depending on what your query combine is
			// i.e. for min query, replace with st[p] = st[p] + lazy[p], since all values in the range will grow by lazy[p]
			// default is for summing
			st[p] = update_comb(st[p], lazy[p] * (R - L + 1));
			
			if(L != R) {
				if(lazy_set[left(p)] != SET_DEFAULT){
					// change to update combine
					lazy_set[left(p)] = update_comb(lazy[p], lazy_set[left(p)]);
					lazy[left(p)] = LAZY_DEFAULT;
				}else{
					// change to update combine
					lazy[left(p)] = update_comb(lazy[left(p)], lazy[p]);
				}
				if(lazy_set[right(p)] != SET_DEFAULT){
					// change to update combine
					lazy_set[right(p)] = update_comb(lazy[p], lazy_set[right(p)]);
					lazy[right(p)] = LAZY_DEFAULT;
				}else{
					// change to update combine
					lazy[right(p)] = update_comb(lazy[right(p)], lazy[p]);
				}
				
			}
			lazy[p] = LAZY_DEFAULT;
		}
	}
	
	// for range/point updates (point updates i == j)
	// takes twice as long as normal point update seg tree, but same time complexity
	void update(int p, int L, int R, int i, int j, int val) {
		push(p, L, R);
		// completely outside the segment
		if(i > R || j < L) return;
		// fully inside the segment
		if(L >= i && R <= j) {
			if(lazy_set[p] != SET_DEFAULT){
				lazy[p] = LAZY_DEFAULT;
				// update combine
				lazy_set[p] = update_comb(lazy_set[p], val);
			}else{
				lazy[p] = val;
			}
			push(p, L, R);
			return;
		}
		// half in half out of segment
		update(left(p), L, (L + R)/2, i, j, val);
		update(right(p), (L + R)/2 + 1, R, i, j, val);
		// change to query combine
		st[p] = query_comb(st[left(p)], st[right(p)]);
	}
	
	// for range/point sets (point updates i==j)
	// takes twice as long as normal point update seg tree, but same time complexity
	void set(int p, int L, int R, int i, int j, int val) {
		push(p, L, R);
		// completely outside the segment
		if(i > R || j < L) return;
		// fully inside the segment
		if(L >= i && R <= j) {
			lazy_set[p] = val;
			lazy[p] = LAZY_DEFAULT;
			push(p, L, R);
			return;
		}
		// half in half out of segment
		set(left(p), L, (L + R)/2, i, j, val);
		set(right(p), (L + R)/2 + 1, R, i, j, val);
		// change to query combine
		st[p] = query_comb(st[left(p)], st[right(p)]);
	}
	
	ll query(int p, int L, int R, int i, int j) {
		push(p, L, R);
		// completely outside the segment
		// Default always loses in comb()
		if(i > R || j < L) return INIT_DEFAULT;
		// fully inside the segment
		if(L >= i && R <= j) return st[p];
		// change to query combine
		return query_comb(query(left(p), L, (L + R)/2, i, j), query(right(p), (L + R)/2 + 1, R, i, j));
	}
	
	// For debugging purposes
	void toString(){
		cout << endl;
		for(int i = 0; i < 2*size; i++){
			cout << st[i] << " ";
		}
		cout << endl;
		for(int i = 0; i < 2*size; i++){
			cout << lazy[i] << " ";
		}
		cout << endl;
		for(int i = 0; i < 2*size; i++){
			cout << lazy_set[i] << " ";
		}
		cout << endl;
	}
};
