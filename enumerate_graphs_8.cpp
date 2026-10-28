// Self-contained exhaustive (isomorph-free) graph enumeration for n=8.
// All 2^21 seven-vertex graphs are enumerated; connected ones are canonically
// deduplicated. Attaching an eighth vertex to every nonempty neighbor subset
// yields every connected eight-vertex graph (remove a spanning tree leaf).
// Canonicalization explicitly considers every permutation that preserves the
// sorted degree sequence. Other permutations cannot attain a canonical form
// with sorted degree sequence.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <unordered_set>
#include <vector>
#include <cassert>
using namespace std;
using U=uint32_t;
struct Canon {
 unsigned char adj[8][8]{};
 int ord[8]{}, blockend[8]{}, groups{}, n{};
 U best;
 U pack(){
  U out=0;int k=0;
  for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++k)
   out|=U(adj[ord[i]][ord[j]])<<k;
  return out;
 }
 void permute(int g){
  if(g==groups){best=min(best,pack());return;}
  const int a=g?blockend[g-1]:0,b=blockend[g];
  do{permute(g+1);}while(next_permutation(ord+a,ord+b));
 }
 U operator()(U mask,int vertices){
  n=vertices;int deg[8]{};for(auto &row:adj)for(auto &v:row)v=0;
  int k=0;
  for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++k)
   if(mask>>k&1){adj[i][j]=adj[j][i]=1;deg[i]++;deg[j]++;}
  iota(ord,ord+n,0);
  sort(ord,ord+n,[&](int a,int b){return deg[a]==deg[b]?a<b:deg[a]<deg[b];});
  groups=0;
  for(int i=1;i<=n;i++)if(i==n||deg[ord[i]]!=deg[ord[i-1]])blockend[groups++]=i;
  best=UINT32_MAX;permute(0);return best;
 }
};
bool connected(U mask,int n){
 U neighbors[8]{};int k=0;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(mask>>k&1){neighbors[i]|=1u<<j;neighbors[j]|=1u<<i;}
 U seen=1;
 while(true){U grow=seen;for(int i=0;i<n;i++)if(seen>>i&1)grow|=neighbors[i];if(grow==seen)break;seen=grow;}
 return seen==((1u<<n)-1);
}
int main(){
 Canon ca;
 unordered_set<U> seven;seven.reserve(1200);
 for(U m=0;m<(1u<<21);++m){
  if(connected(m,7))seven.insert(ca(m,7));
 }
 cerr<<"7-vertex connected unlabeled classes: "<<seven.size()<<"\n";
 if(seven.size()!=853)return 1;
 unordered_set<U> eight;eight.reserve(12000);
 for(U m:seven){
  unsigned char adj[7][7]{};int k=0;
  for(int i=0;i<7;i++)for(int j=i+1;j<7;j++,k++)adj[i][j]=adj[j][i]=(m>>k&1);
  for(int add=1;add<128;add++){
   U mask=0;k=0;
   for(int i=0;i<8;i++)for(int j=i+1;j<8;j++,k++)
    if(j==7?(add>>i&1):adj[i][j])mask|=1u<<k;
   eight.insert(ca(mask,8));
  }
 }
 cerr<<"8-vertex connected unlabeled classes: "<<eight.size()<<"\n";
 if(eight.size()!=11117)return 2;
 vector<U> keys(eight.begin(),eight.end());sort(keys.begin(),keys.end());
 for(U v:keys)cout<<v<<"\n";
}
