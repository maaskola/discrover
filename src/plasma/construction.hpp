/*
 * This code is based on source code that accompanied:
 * Linear Work Suffix Array Construction
 * Juha Kärkkäinen
 * Peter Sanders
 * Stefan Burkhardt
 */

#include <vector>

template <typename X, typename Y>
inline bool leq(X a1, Y a2, X b1, Y b2) { // lexic. order for pairs
  return(a1 < b1 || (a1 == b1 && a2 <= b2)); 
}
template <typename X, typename Y, typename Z>
inline bool leq(X a1, Y a2, Z a3, X b1, Y b2, Z b3) {
  return(a1 < b1 || (a1 == b1 && leq(a2,a3, b2,b3))); 
};
// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
template <typename idx_t, typename Iter>
void radixPass(const std::vector<idx_t> &a, std::vector<idx_t> &b, Iter r, idx_t n, idx_t K) 
{ // count occurrences
  std::vector<idx_t> c(K + 1, 0);             // counter array
  for(idx_t i = 0; i < n; i++)
    c[r[a[i]]]++;                             // count occurrences
  for(idx_t k = 0, sum = 0; k <= K; k++) {
    idx_t t = c[k]; c[k] = sum; sum += t;     // exclusive prefix sums
  }
  for(idx_t i = 0; i < n;  i++)
    b[c[r[a[i]]]++] = a[i];                   // sort
};

// find the suffix array SA of s[0..n-1] in {1..K}^n
// require s[n]=s[n+1]=s[n+2]=0, n>=2
template <typename idx_t, typename Iter>
void suffixArray(Iter begin, Iter end, std::vector<idx_t> &SA, idx_t n, idx_t K) {
  idx_t n0=(n+2)/3, n1=(n+1)/3, n2=n/3, n02=n0+n2; 
  std::vector<idx_t> s12(n02 + 3);  s12[n02]= s12[n02+1]= s12[n02+2]=0; 
  std::vector<idx_t> SA12(n02 + 3); SA12[n02]=SA12[n02+1]=SA12[n02+2]=0;
  std::vector<idx_t> s0(n0);
  std::vector<idx_t> SA0(n0);

  // generate positions of mod 1 and mod  2 suffixes
  // the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
  for (idx_t i=0, j=0;  i < n+(n0-n1);  i++) if (i%3 != 0) s12[j++] = i;

  // least-significant-bit radix sort the mod 1 and mod 2 triples
  radixPass(s12 , SA12, begin+2, n02, K);
  radixPass(SA12, s12 , begin+1, n02, K);  
  radixPass(s12 , SA12, begin  , n02, K);

  // find lexicographic names of triples
  idx_t name = 0;
  typename Iter::value_type c0 = 0, c1 = 0, c2 = 0;
  bool first = true;
  for (idx_t i = 0;  i < n02;  i++) {
    if (first or begin[SA12[i]] != c0 || begin[SA12[i]+1] != c1 || begin[SA12[i]+2] != c2) { 
      name++;  c0 = begin[SA12[i]];  c1 = begin[SA12[i]+1];  c2 = begin[SA12[i]+2]; first = false;
    }
    if (SA12[i] % 3 == 1) { s12[SA12[i]/3]      = name; } // left half
    else                  { s12[SA12[i]/3 + n0] = name; } // right half
  }

  // recurse if names are not yet unique
  if (name < n02) {
    suffixArray(s12.begin(), s12.end(), SA12, n02, name);
    // suffixArray(s12, SA12, n02, name);
    // store unique names in s12 using the suffix array 
    for (idx_t i = 0;  i < n02;  i++) s12[SA12[i]] = i + 1;
  } else // generate the suffix array of s12 directly
    for (idx_t i = 0;  i < n02;  i++) SA12[s12[i] - 1] = i; 

  // stably sort the mod 0 suffixes from SA12 by their first character
  for (idx_t i=0, j=0;  i < n02;  i++) if (SA12[i] < n0) s0[j++] = 3*SA12[i];
  radixPass(s0, SA0, begin, n0, K);

  // merge sorted SA0 suffixes and sorted SA12 suffixes
  for (idx_t p=0,  t=n0-n1,  k=0;  k < n;  k++) {
#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
    idx_t i = GetI(); // pos of current offset 12 suffix
    idx_t j = SA0[p]; // pos of current offset 0  suffix
    if (SA12[t] < n0 ? 
        leq(begin[i],       s12[SA12[t] + n0], begin[j],       s12[j/3]) :
        leq(begin[i],begin[i+1],s12[SA12[t]-n0+1], begin[j],begin[j+1],s12[j/3+n0]))
    { // suffix from SA12 is smaller
      SA[k] = i;  t++;
      if (t == n02) { // done --- only SA0 suffixes left
        for (k++;  p < n0;  p++, k++) SA[k] = SA0[p];
      }
    } else { 
      SA[k] = j;  p++; 
      if (p == n0)  { // done --- only SA12 suffixes left
        for (k++;  t < n02;  t++, k++) SA[k] = GetI(); 
      }
    }  
  } 
};

