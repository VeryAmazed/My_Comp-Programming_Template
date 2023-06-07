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