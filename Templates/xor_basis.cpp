// Get basis vectors
// Must note that the vectors you get can be a xor combination of your input vectors
// Vectors already in the basis pass off their smaller components to the newer vector
// i.e. if A can be a basis for 2^4 and 2^1 and A is already in the basis, then when B is introduced
// and can only be a basis for 2^4, when you take min(B, A^B) A^B is a basis for 2^1
vector<int> basis;
auto add_base = [&](int x){
	for (int i = 0; i < (int)basis.size(); i++) {  // reduce x using the current basis vectors
		x = min(x, x ^ basis[i]);
	}

	if (x != 0) { basis.push_back(x); }
};
for(auto a : vec){
	add_base(a);
}

// Transform a set of basis vectors into reduced canocical form 
// First loop takes your basis vectors and puts them into their position by index of most signficant bit 
// Second loop goes through the vectors in order of largest most sig bit to smallest least sig bit 
// and tries to reduce the vector with the larger most sig bit by xoring it with vectors with smaller least sig bits
vector<int> base(max_bit+1);
for(int x : basis){
    for(int b = max_bit; b >= 0; b--){i
        if(!(x & (1 << b))) continue;
        if(!base[b]){ base[b] = x; break;}
    }
}
for (int b = max_bit; b >= 0; --b) if (base[b]) {
    for (int c = b-1; c >= 0; --c)
        if (base[c] && (base[b] & (1<<c))) base[b] ^= base[c];
}
vector<int> bvec;
for (int b = 0; b <= max_bit; ++b) if (base[b]) bvec.push_back(base[b]);
