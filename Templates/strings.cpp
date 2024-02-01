// Z-algorithm, generates an array z, where z[i] is the length of a substring substr, where substr is a prefix of the string str, 
// and it starts at index i in str
// can be used be used to find in number of occurences of a string a, as a substring of another string b, simply input a + '$' + b, 
// as the parameter, '$' is any character that will not be in string A or B
// doesn't need to take a string as input, because a string is just a vector of characters, which is just a vector is integers, 
// you can also just input a vector of integers 
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
ll M1 = 1e9+7;
ll M2 = 1e9+9;
ll M3 = 1000000021;
const ll P1 = (uint64_t)chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now()
.time_since_epoch()).count() % (uint64_t)M1;
const ll P2 = (uint64_t)chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now()
.time_since_epoch()).count() % (uint64_t)M2;
const ll P3 = (uint64_t)chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now()
.time_since_epoch()).count() % (uint64_t)M3;

vector<ll> pow1(str.size()+1, 0);
vector<ll> pow2(str.size()+1, 0);
vector<ll> pow3(str.size()+1, 0);
vector<ll> poly1(str.size()+1, 0);
vector<ll> poly2(str.size()+1, 0);
vector<ll> poly3(str.size()+1, 0);
pow1[0] =1;
pow2[0] =1;
pow3[0] =1;
for(int i =0; i < n; i++){
	pow1[i+1] = (pow1[i]*P1)%M1;
	pow2[i+1] = (pow2[i]*P2)%M2;
	pow3[i+1] = (pow3[i]*P3)%M3;
	poly1[i+1] = ((poly1[i]*P1)%M1 + str[i])%M1;
	poly2[i+1] = ((poly2[i]*P2)%M2 + str[i])%M2;
	poly3[i+1] = ((poly3[i]*P3)%M3 + str[i])%M3;
}
