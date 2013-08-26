#ifndef VMATH_H
#define VMATH_H

#include <algorithm>
#include <vector>
#include <math.h>

class VMath {
public:
	template <typename T> 
	static void vpow (vector<T> &v, double exp) {
		for (int i = 0; i < v.size(); i++)	v[i] = (T) pow(v[i], exp);
	}
	template <typename T1, typename T2>
	static void add (vector<T1> &ans, vector<T2> &v, bool cont=false) {
		int sz, i;
		if (cont) {
			T1 l = v.size()?v[v.size()-1]:0;
			for (i = v.size(); i < ans.size(); i++)	v.push_back(l);
			T2 s = ans.size()?ans[ans.size()-1]:0;
			for (i = ans.size(); i < v.size(); i++)	ans.push_back(s);
		}
		sz = min (ans.size(), v.size());
		for (i = 0; i < sz; i++)	ans[i] += v[i];
		for (; i < v.size(); i++)	ans.push_back(v[i]);
	}
	template <typename T1, typename T2>
	static void sub (vector<T1> &ans, vector<T2> &v, bool cont=false) {
		int sz, i;
		if (cont) {
			T1 l = v.size()?v[v.size()-1]:0;
			for (i = v.size(); i < ans.size(); i++)	v.push_back(l);
			T2 s = ans.size()?ans[ans.size()-1]:0;
			for (i = ans.size(); i < v.size(); i++)	ans.push_back(s);
		}
		sz = min (ans.size(), v.size());
		for (i = 0; i < sz; i++)	ans[i] -= v[i];
		for (; i < v.size(); i++)	ans.push_back(-v[i]);
	}
	template <typename T1, typename T2>
	static void div (vector<T1> &ans, T2 val) {
		for (int i = 0; i < ans.size(); i++)	ans[i] /= val;
	}
	template <typename T>
	static T vmax (vector<T> v) {
		T ans = v[0];
		for (int i = 1; i < v.size(); i++)	ans = max (ans, v[i]);
		return ans;
	}
};
#endif
