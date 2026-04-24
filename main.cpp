// Problem 081 - Vexoben 的橙汁加工厂
// Computes sum over all a<b of max flow between a and b.
// Undirected unit-capacity edges; degree<=3, n<=3000, m<=4500.
// Use Gomory–Hu tree construction with Dinic as maxflow oracle.

#include <bits/stdc++.h>
using namespace std;

struct Dinic {
    struct Edge { int to; int cap; int rev; };
    int N;
    vector<vector<Edge>> G;
    vector<int> level, it;
    Dinic(int n=0): N(n), G(n), level(n), it(n) {}
    void reset(int n){ N=n; G.assign(n,{}); level.assign(n,0); it.assign(n,0);}    
    inline void add_edge(int u,int v,int c){
        G[u].push_back(Edge{v,c,(int)G[v].size()});
        G[v].push_back(Edge{u,0,(int)G[u].size()-1});
    }
    bool bfs(int s,int t){
        fill(level.begin(), level.end(), -1);
        static int qbuf[5005]; // n<=3000
        int qb=0, qe=0;
        level[s]=0; qbuf[qe++]=s;
        while(qb<qe){
            int u=qbuf[qb++];
            for(const auto &e:G[u]) if(e.cap>0 && level[e.to]<0){ level[e.to]=level[u]+1; qbuf[qe++]=e.to; }
        }
        return level[t]>=0;
    }
    int dfs(int u,int t,int f){
        if(u==t) return f;
        for(int &i=it[u]; i<(int)G[u].size(); ++i){
            Edge &e=G[u][i];
            if(e.cap>0 && level[e.to]==level[u]+1){
                int ret=dfs(e.to,t,min(f,e.cap));
                if(ret>0){ e.cap-=ret; G[e.to][e.rev].cap+=ret; return ret; }
            }
        }
        return 0;
    }
    long long maxflow(int s,int t){
        long long flow=0;
        while(bfs(s,t)){
            fill(it.begin(), it.end(), 0);
            while(int f=dfs(s,t,INT_MAX)) flow+=f;
        }
        return flow;
    }
};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n,m; if(!(cin>>n>>m)) return 0; 
    vector<pair<int,int>> edges; edges.reserve(m);
    vector<int> deg(n,0);
    vector<vector<int>> adj(n);
    for(int i=0;i<m;++i){
        int a,b; cin>>a>>b; --a; --b;
        edges.emplace_back(a,b);
        deg[a]++; deg[b]++;
        adj[a].push_back(b); adj[b].push_back(a);
    }

    // Trivial cases
    if(n<=1){ cout<<0<<"\n"; return 0; }
    if(m==0){ cout<<0<<"\n"; return 0; }

    // Build Gomory–Hu tree using Dinic oracle on unit capacities.
    // Complexity: (n-1) maxflow computations. Each maxflow on sparse unit graph is reasonable.
    Dinic din(n);

    auto build_graph = [&](void){
        din.reset(n);
        for(auto &e: edges){
            int u=e.first, v=e.second;
            din.add_edge(u,v,1);
            din.add_edge(v,u,1);
        }
    };

    vector<int> parent(n,0);
    vector<long long> cutval(n,0);
    for(int i=1;i<n;++i) parent[i]=0; // root parent of others
    for(int s=1; s<n; ++s){
        int t=parent[s];
        build_graph();
        long long f = din.maxflow(s,t);
        cutval[s]=f;
        // BFS in residual graph to find reachable set from s
        vector<char> vis(n,0);
        queue<int> q; q.push(s); vis[s]=1;
        while(!q.empty()){
            int u=q.front(); q.pop();
            for(const auto &e: din.G[u]) if(e.cap>0 && !vis[e.to]){ vis[e.to]=1; q.push(e.to);}            
        }
        for(int v=s+1; v<n; ++v){ if(parent[v]==t && vis[v]) parent[v]=s; }
        if(vis[parent[t]]){ parent[s]=parent[t]; parent[t]=s; swap(cutval[s], cutval[t]); }
    }

    // Build GH tree adjacency list and list of edges
    vector<tuple<int,int,int>> tedges; tedges.reserve(n-1);
    vector<vector<pair<int,int>>> tree(n);
    for(int v=1; v<n; ++v){
        int u=parent[v]; int w=(int)cutval[v];
        tree[u].push_back({v,w});
        tree[v].push_back({u,w});
        tedges.emplace_back(u,v,w);
    }

    // Sum over all unordered pairs of min edge weight along their path in the tree.
    // This equals processing edges in non-decreasing weight and union components:
    // contribution += w * size[a] * size[b] when connecting components of sizes a,b.
    struct DSU{
        vector<int> p, sz;
        DSU(int n=0){init(n);}        
        void init(int n){ p.resize(n); sz.assign(n,1); iota(p.begin(), p.end(), 0);}        
        int find(int x){ return p[x]==x?x:p[x]=find(p[x]); }
        bool unite(int a,int b,long long &ans,int w){ a=find(a); b=find(b); if(a==b) return false; ans += 1LL*w*sz[a]*sz[b]; if(sz[a]<sz[b]) swap(a,b); p[b]=a; sz[a]+=sz[b]; return true; }
    } dsu(n);
    sort(tedges.begin(), tedges.end(), [](auto &A, auto &B){ return get<2>(A) > get<2>(B); });
    long long answer=0;
    for(auto &e: tedges){ int u,v,w; tie(u,v,w)=e; dsu.unite(u,v,answer,w); }

    cout << answer << "\n";
    return 0;
}
