#pragma once
#include "Distribution.h"
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;
using namespace arma;
using namespace Rcpp;
class BosPredict
{

public:
	BosPredict(int kr, int kc, int m, mat pis, mat mus);
	BosPredict();
	~BosPredict();
	mat missingValuesInit(mat& x);
	mat SEstep_predict(mat W, mat x);
	cube getCubeProbs();
	cube gettabpej();
protected:
	int _m;
	mat _pis;
	mat _mus;
	int _kr;
	int _kc;
	vector<vector<int>> _miss;
	cube _tab_pejs;
	random_device _rd; // to sample distributions
};

