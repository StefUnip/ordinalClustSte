#pragma once
#include "Distribution.h"
#include "Bos.h"

// [[Rcpp::depends(RcppArmadillo)]] 
#include <armadillo>
#include <limits>
#include <cmath>
#include <list>
#include <typeinfo>
#include <iostream>
#include <initializer_list> 
#include <vector>
#include <numeric>


using namespace arma;
using namespace std;
class CoClusteringContext
{
public:
	CoClusteringContext(arma::mat x, std::vector< arma::urowvec > dlist, int kr, std::vector< int > kc,
		std::string init, int nbSEM, int nbSEMburn, int nbindmini, std::vector< int > m);
	CoClusteringContext();
	~CoClusteringContext();
	void missingValuesInit();
	bool initialization();
	void Mstep();
	void MstepVW();
	void SEstep();
	void SEstepRow();
	void SEstepCol();
	void sampleV();
	void sampleW();
	void sampleVW();
	void sampleVWStock();
	void imputeMissingData();
	bool verif();
	void fillParameters(int iteration);
	void fillLabels(int iteration);
	void getBurnedParameters();
	void printResults();
	void returnResults();
	void  putParamsToZero();
	S4 returnCoclustering();
	double computeICL();

protected:
	mat _x;
	int _Nr;
	vector<int> _Jc;
	vector<int>_m;
	vector<urowvec> _dlist;
	vector<Distribution*> _distrib_objects;
	int _number_distrib;
	int _kr;
	vector<int> _kc;
	vector<int> _zr;
	vector<vector<int>> _zc;
	mat _probaV;
	vector<mat> _probaW;
	mat _logprobaV;
	vector<mat> _logprobaW;
	mat _V;
	vector<mat> _W;
	vector<rowvec> _rho;
	rowvec _gamma;
	vector<vector<rowvec>> _allrho;
	vector<rowvec> _allgamma;
	vector<rowvec> _resrho;
	rowvec _resgamma;
	// params of co-clustering :
	string _init;
	int _nbSEM;
	int _nbSEMburn;
	int _nbindmini;

	// to sample distributions
	random_device _rd;

	//V and W init
	vector<int> _zrinit;
	vector<int> _zcinit;

	// zr and zc chains
	mat _zrchain;
	std::vector<mat> _zcchain;

	double _icl;

	//utils

	rowvec getMeans(mat VorW);
	double logsum(rowvec logx);
	mat kmeansi();
	double getDistance(vec &a, vec &b);

};

