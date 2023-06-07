// monkey8's Lazy Prop Seg Tree Template
// Modified by Me
// Modified again to to allow for different updates and query types

// p is curent pos in array, always start at 1
// L, R is current bounds of the segment, always start at 0, size-1
// i, j is the range you are querying, both inclusive
// default version of this tempalte is for adding

// Modify to be Templated, w/ 2 different types, one for vector and one for lazy, like below

template <class T, class U> 
struct Lazy_Seg_T{
	// modify the default value of the lazy array based on what is needed
	const U LAZY_DEFAULT = 0;
	
	// modify the default value of the seg tree
	const T INIT_DEFAULT = 0;
	
	vector<T> st;
	vector<U> lazy;
	int size;
	
	Lazy_Seg_T(int n, vector<T>& in){
		size = n;
		while((size&(size-1))){
			size++;
		}
		
		st.assign(2*size, INIT_DEFAULT);
		lazy.assign(2*size, LAZY_DEFAULT);
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
		
	}
	
	// change depending on how you want to combine for queries (i.e. min query, sum)
	// Template ver is for summing
	T query_comb(T a, T b){
		return a + b;
	}
	// change depending on how you want to combine for updates (i.e. add to a range)
	// template is for summing
	U update_comb(U a, U b){
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
	
	void push(int p, int L, int R) {
		
		if(lazy[p] != LAZY_DEFAULT){
			// change depending on what your query combine is
			// i.e. for min query, replace with st[p] = st[p] + lazy[p], since all values in the range will grow by lazy[p]
			// default is for summing
			st[p] += lazy[p] * (R - L + 1);
			
			if(L != R) {
				lazy[left(p)] = update_comb(lazy[left(p)], lazy[p]);
				
				lazy[right(p)] = update_comb(lazy[right(p)], lazy[p]);
			}
			lazy[p] = LAZY_DEFAULT;
		}
	}
	
	// for range/point updates (point updates i == j)
	// takes twice as long as normal point update seg tree, but same time complexity
	void update(int p, int L, int R, int i, int j, U val) {
		push(p, L, R);
		// completely outside the segment
		if(i > R || j < L) return;
		// fully inside the segment
		if(L >= i && R <= j) {
			
			lazy[p] = val;
			
			push(p, L, R);
			return;
		}
		// half in half out of segment
		update(left(p), L, (L + R)/2, i, j, val);
		update(right(p), (L + R)/2 + 1, R, i, j, val);
		// change to query combine
		st[p] = query_comb(st[left(p)], st[right(p)]);
	}
	
	T query(int p, int L, int R, int i, int j) {
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
	}
};
