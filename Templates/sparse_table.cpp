// Templated Sparse Table
// Only works for operations where double counting doesn't matter (Range Min/Max Queries)
template<typename T>
struct Sparse{
	// size of array to be queried
	int size;
	// value of the largest k s.t. 2^k <= size
	int k;
	// sparse table
	vector<vector<T>> st;
	// precomputed logs
	vector<int> lg;
	
	Sparse(vector<T>& in){
			size = in.size();
			
			// determine size of k
			k =0;
			int mult = 1;
			while(mult <= size){
				k++;
				mult *= 2;
			}
			
			// precompute sparse table
			st.assign(k+1, vector<T>());
			for(int i =0; i < size; i++){
				st[0].pb(in[i]);
			}
			for(int i = 1; i <= k; i++){
				for(int j =0; j + (1 << i) <= size; j++){
					st[i].pb(comp(st[i-1][j], st[i-1][j+ (1 << (i-1))]));
				}
			}
			
			// precompute all logs
			lg.assign(size+1, 0);
			lg[1] = 0;
			for(int i =2; i <= size; i++){
				lg[i] = lg[i/2] + 1;
			}
	}
	
	// comparison function
	// modify as needed
	T comp(T a, T b){
		return min(a, b);
	}
	
	// L and R are the 0-indexed indices we want to query between, both inclusive
	T query(int L, int R){
		int i = lg[R-L+1];
		return comp(st[i][L], st[i][R-(1<<i)+1]);
	}
};