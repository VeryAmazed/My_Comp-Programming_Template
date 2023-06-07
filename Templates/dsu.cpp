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