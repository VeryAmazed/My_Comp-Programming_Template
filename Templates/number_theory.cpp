
// a and b must be positive because mod in c++ does not follow the mathematical definition
ll gcd(ll a, ll b) {
	if (b == 0) return a;
	return gcd(b, a % b);
}

// extended euclidean algorithm
// a and b must be positive because mod in c++ does not follow the mathematical definition
// solves ax + by = gcd(a,b)
ll gcd_ext(ll a, ll b, ll& x, ll& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    ll x1, y1;
    ll d = gcd_ext(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}

// finds a solution to a linear diophantine equation
// uses the above function
// a and b can be negative 
// solves ax0 + by0 = c and g = gcd(a, b)
// the solutions returned x0 and y0 are just one solution to the equation, 
// the general solutions, x and y can be expressed as x = x0 + k(b/g) and y = y0 - k(a/g), 
// we see that plugging these 2 expressions back into the euqation doesn't break the inequality
// for more info, refer to the cp-algos website
bool find_any_solution(ll a, ll b, ll c, ll &x0, ll &y0, ll &g) {
    g = gcd_ext(abs(a), abs(b), x0, y0);
    if (c % g) {
        return false;
    }

    x0 *= c / g;
    y0 *= c / g;
    if (a < 0) x0 = -x0;
    if (b < 0) y0 = -y0;
    return true;
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
// find mod inverse of a number. Uses eulaer's and fermat's little theorem so only works when 
// mod is prime, if m is prime x^(-1) = x^(m-2)
ll mod_inv(ll a, ll mod){
	if(gcd(a, mod) == 1){
		return binpow(a , mod-2);
	}
	else{
		return -1; // a doesn't have a mod inverse in this mod space
	}
}

// Sieve of Eratosthenes
// Primes are in the primes vector
// Also allows for the getting of all prime factors and their powers of some given number a[i] < mxn
// Basically we set primes to themselves and every composite number we store the smallest prime in it's 
// prime factorization then to get the prime factors just use the while loop below
// Basically inductively, for every number that is composed purely of primes less than the current prime, 
// we have a correct prime factorization
// Then for the current prime, the division either takes you to another number that is marked with the 
// current prime, or to a number for which following the path will produce a correct prime factorization
ll sieve[mxn+1] = {0};
vector<ll> primes;
for (ll i = 2; i <= mxn; i++) {
	if (sieve[i] != 0) continue;
	sieve[i] = i;
	primes.pb(i);
	for (ll j = 2*i; j <= mxn; j+=i) {
		if (sieve[j] != 0) continue;
		sieve[j] = i;
	}
}

while (a[i] > 1) {
	a[i] /= sieve[a[i]];
}

// From -is-this-fft-
// Fast integer sqrt that works like a binary search, i.e. O(logn)
ll int_sqrt (ll x) {
  ll ans = 0;
  for (ll k = (ll)1 << 31; k != 0; k /= 2) {
    if ((ans + k) * (ans + k) <= x) {
      ans += k;
    }
  }
  return ans;
}
