// Problem 081 - Vexoben 的橙汁加工厂
// Computes sum over all a<b of max flow between a and b.
// Undirected unit-capacity edges; degree<=3, n<=3000, m<=4500.
// Use Gomory–Hu tree construction with Dinic as maxflow oracle.

#include <bits/stdc++.h>
using namespace std;

struct Dinic {
    int N;
    vector<int> head, to, cap, nxt;
    vector<int> level, it;
    int ec;
    Dinic(int n=0){ reset(n); }
    void reset(int n){ N=n; head.assign(n,-1); to.clear(); cap.clear(); nxt.clear(); level.assign(n,0); it.assign(n,0); ec=0; to.reserve(20000); cap.reserve(20000); nxt.reserve(20000);}    
    inline void add_dir(int u,int v,int c){ to.push_back(v); cap.push_back(c); nxt.push_back(head[u]); head[u]=ec++; }
    inline void add_edge(int u,int v,int c){ add_dir(u,v,c); add_dir(v,u,0); }
    bool bfs(int s,int t){
        fill(level.begin(), level.end(), -1);
        static int qbuf[6005]; // n<=3000
        int qb=0, qe=0;
        level[s]=0; qbuf[qe++]=s;
        while(qb<qe){
            int u=qbuf[qb++];
            for(int e=head[u]; e!=-1; e=nxt[e]) if(cap[e]>0){ int v=to[e]; if(level[v]<0){ level[v]=level[u]+1; qbuf[qe++]=v; if(v==t) return true; }}
        }
        return level[t]>=0;
    }
    int dfs(int u,int t,int f){
        if(u==t) return f;
        for(int &e=it[u]; e!=-1; e=nxt[e]) if(cap[e]>0){ int v=to[e]; if(level[v]==level[u]+1){ int ret=dfs(v,t,min(f,cap[e])); if(ret>0){ cap[e]-=ret; cap[e^1]+=ret; return ret; } } }
        return 0;
    }
    long long maxflow(int s,int t,int limit=INT_MAX){
        long long flow=0;
        while(flow<limit && bfs(s,t)){
            it = head;
            int add;
            while(flow<limit && (add=dfs(s,t,limit - (int)flow))>0) flow+=add;
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
        din.to.reserve(4*(int)edges.size()+5); din.cap.reserve(4*(int)edges.size()+5); din.nxt.reserve(4*(int)edges.size()+5);
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
        int lim = min(3, min(deg[s], deg[t]));
        long long f = din.maxflow(s,t,lim);
        cutval[s]=f;
        // BFS in residual graph to find reachable set from s
        vector<char> vis(n,0);
        queue<int> q; q.push(s); vis[s]=1;
        while(!q.empty()){
            int u=q.front(); q.pop();
            for(int e=din.head[u]; e!=-1; e=din.nxt[e]) if(din.cap[e]>0){ int v=din.to[e]; if(!vis[v]){ vis[v]=1; q.push(v);} }
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
