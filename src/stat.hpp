#ifndef STAT_H
#define STAT_H
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <math.h>
#include <typeinfo>

using namespace std;

class Statistic {
public:
	template <typename T>
	static void normalize (map<T, double> &ans) {
		double sum = 0;
		typename map<T, double>::iterator it;
		for (it = ans.begin(); it != ans.end(); it++)
			sum += it->second;
		for (it = ans.begin(); it != ans.end(); it++)
			it->second /= sum;
	}
	template <typename T>
	static void getDistribution (const vector<T> &v, map<T, double> &ans) {
		typename map<T, double>::iterator it;
		for (int i = 0; i < v.size(); i++)	ans[v[i]]++;
		normalize(ans);
	}
	template <typename T>
	static void getFrequency (const vector<T> &v, map<T, double> &ans) {
		typename map<T, double>::iterator it;
		for (int i = 0; i < v.size(); i++)	ans[v[i]]++;
	}
	template <typename T>
	static void getDistribution (const vector<vector<T> > &M, map<T, double> &ans) {
		typename map<T, double>::iterator it;
		for (int i = 0; i < M.size(); i++)	
			for (int j = 0; j < M[i].size(); j++)	ans[M[i][j]]++;
		normalize(ans);
	}
	template <typename T>
	static void getFrequency (const vector<vector<T> > &M, map<T, double> &ans) {
		typename map<T, double>::iterator it;
		for (int i = 0; i < M.size(); i++)	
			for (int j = 0; j < M[i].size(); j++)	ans[M[i][j]]++;
	}
	
	template <typename T>
	static double getMoment (int order, map<T, double> & dist) {
		typename map<T, double>::iterator it;
		double ans = 0;
		for (it = dist.begin(); it != dist.end(); it++)
			ans += pow(it->first, order)*it->second;
		return ans;
	}
	template <typename T>
	static double getShannonEntropy (map<T, double> dist) {
		typename map<T, double>::iterator it;
		double ans = 0, log2 = log(2);
		for (it = dist.begin(); it != dist.end(); it++)
			ans += it->second*log(it->second);
		ans /= -log2;
		return ans;
	}
	static double unif (void) {
		return rand()/(double)RAND_MAX;
	}
};
#endif
