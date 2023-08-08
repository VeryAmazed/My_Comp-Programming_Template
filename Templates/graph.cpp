// Everything in here probably works

// Djikstra's (from CPH)
vector<vector<pll>> vec; // the graph
vector<ll> dists;
vector<bool> visited;

void djik(int start){
	dists[start] = 0;
	priority_queue<pll, vector<pll>, greater<pll>> pq;
	pq.push(mp(0, start));
	while (!pq.empty()) {
		int a = pq.top().s; 
		pq.pop();
		if (visited[a]) continue;
		visited[a] = true;
		for (auto e : vec[a]) {
			int b = e.first, w = e.second;
			if (dists[a]+w < dists[b]) {
				dists[b] = dists[a]+w;
				pq.push(mp(dists[b],b));
			}
		}
	}
}

// Bellman-Ford (from CPH)
/*
on each of the iterations from i to n-1, dists will contain the shortest path 
from the starting node to any node that goes through i edges or less
*/
// To check for a negative cycle, go through the edges an n-th time and if any distance gets reduced, there is a negative cycle
vector<vector<pll>> vec; // the graph
vector<ll> dists;

dists[x] = 0;
for (int i = 1; i <= n-1; i++) {
	for(int a = 1; a <= n; a++){
		for(auto e : vec[a]){
			b = e.f;
			w = e.s;
			distance[b] = min(distance[b], distance[a]+w);
		}
	}
}

// Dense Djikstra's (from cp-algos)
/*
On every iteration from 1 to n, we find the node that has the smallest dist from start
and attempt to reduce the distance of any node it has an edge with.
Basically how the Djikstra's algorithm is visualized. You have a set of nodes which
have the minimum distance and every iteration you add in a new node into this set and when
you do you perform some relaxations
*/
vector<vector<pair<int, int>>> vec;
vector<ll> dists;
vector<int> paths;
vector<bool> visited;

void denseDjik(int start) {
    dists[start] = 0;
    for (int i = 0; i < n; i++) {
        int v = -1;
        for (int j = 0; j < n; j++) {
            if (!visited[j] && (v == -1 || dists[j] < dists[v]))
                v = j;
        }
		// Means none of the remaining nodes are reachable
        if (dists[v] == INF)
            break;

        visited[v] = true;
        for (auto e : vec[v]) {
            int to = e.first;
            int len = e.second;

            if (dists[v] + len < dists[to]) {
                dists[to] = dists[v] + len;
                paths[to] = v;
            }
        }
    }
}

