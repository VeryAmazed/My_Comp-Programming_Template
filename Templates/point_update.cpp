// Taken from USACO.guide and then modified by me
/** A data structure that can answer point update & range minimum queries. */
template <class T> 
struct SegTree {
	/** The operation to use for combining two elements. (Must be associative) 
		Used for querying (when the parent combines its two children)
	 */
	 // change as necessary
	T query_comb(T a, T b) { 
		return a + b;
	}
	
	// operation for updating element at the specified index
	// change as necessary
	T update_comb(T a, T b) {
		return a + b;
	}
	
	const T DEFAULT = 0;  // Default value, change as necessary

	vector<T> segtree;
	int len;

	SegTree(int len) : len(len), segtree(len * 2, DEFAULT) {}
	
	// look at ASSERTS to see how you should be indexing things
	
	/** Sets the value at ind to val. */
	void set(int ind, T val) {
		// assert(0 <= ind && ind < len);
		ind += len;
		segtree[ind] = val;
		for (; ind > 1; ind /= 2) {
			segtree[ind >> 1] = query_comb(segtree[ind], segtree[ind ^ 1]);
		}
	}
	
	/** updates the value at ind to val. */
	void update(int ind, T val) {
		// assert(0 <= ind && ind < len);
		ind += len;
		segtree[ind] = update_comb(segtree[ind], val);
		for (; ind > 1; ind /= 2) {
			segtree[ind >> 1] = query_comb(segtree[ind], segtree[ind ^ 1]);
		}
	}

	/** queries the range [start, end) */
	T query(int start, int end) {
		// assert(0 <= start && start < len && 0 < end && end <= len);
		T sum = DEFAULT;
		for (start += len, end += len; start < end; start /= 2, end /= 2) {
			if ((start & 1) != 0) { sum = query_comb(sum, segtree[start++]); }
			if ((end & 1) != 0) { sum = query_comb(sum, segtree[--end]); }
		}
		return sum;
	}
};