// size of components are stored as negative values
// any node with a negative value is the representative node
/*
	if the values corresponding to the nodes do not fit in 1 through n
	just use an external unordered_map to map the actual values to spots
	inside the vector. Also note that the mapping of the values to indices 
	does not need to following the relative ordering of the values
*/
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

/*
DSU used for dynamic connectivity.
Doesn't have path compression because dynamic connectivity wouldn't work with it
*/
struct DSU{
	vector<int> sets;
	// the time at which a node becomes no longer a parent of itself
	vector<int> time_changed;
	void init(int n){ 
		sets = vector<int>(n+1 ,-1);
		time_changed = vector<int>(n+1);
	}
	// find root of node x at time t
	int find(int x, int time){
		return (sets[x] < 0  || time_changed[x] > time) ? x : find(sets[x], time);
	}
	bool unite(int x, int y, int time){ // union by size
		x = find(x, time);
		y = find(y, time);
		if(x == y) return 0;
		if(sets[x] > sets[y]) swap(x, y);
		sets[x] += sets[y];
		sets[y] = x;
		time_changed[y] = time;
		return 1;
	}
	
};

/*
Set 1 <= t <= m when you set connections. 
t = 0 handles the case for when a = b. Because it will be less than every value in time_changed. 
*/ 
int get_first_connection(DSU& dsu, int a, int b, int m){
	int l = 0, r = m+1;
	while(l < r){
		int mid = l + (r-l)/2;
		if(dsu.find(a, mid) == dsu.find(b, mid)){
			r = mid;
		}else{
			l = mid+1;
		}
	}
	return l;
}
