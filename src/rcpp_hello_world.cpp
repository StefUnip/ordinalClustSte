// [[Rcpp::depends(RcppArmadillo)]] 
#include <RcppArmadillo.h>
#include <limits>
#include <cmath>

// with devtools install(local=F)

// [[Rcpp::plugins(cpp11)]]
double inf = std::numeric_limits<double>::infinity();

//[[Rcpp::plugins("cpp11")]]


using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}





// [[Rcpp::export]]
int unsigned_to_signed(unsigned x)
{
	if (x <= INT_MAX)
		return static_cast<int>(x);

	if (x >= INT_MIN)
		return static_cast<int>(x - INT_MIN) + INT_MIN;

	throw x; // Or whatever else you like
}
int unsigned_to_signed(unsigned x);
// [[Rcpp::export]]
bool compare_vec(arma::urowvec vec1, arma::rowvec vec2) {
	bool result = true;
	if (vec1.size() != vec2.size()) {
		return(false);
	}
	else {
		for (unsigned int i = 0; i < vec1.size(); ++i) {
			int signed_cast = unsigned_to_signed(vec1(i));
			if (signed_cast != vec2(i)) {
				result = false;
				break;
			}
		}
	}
	return(result);
}
bool compare_vec(arma::urowvec vec1, arma::rowvec vec2);

// [[Rcpp::export]]
arma::umat allej(int j, int m) { // TODO . to implement
	arma::umat result;
	if (j == 1) {
		result << 1 << m << arma::endr;
		return(result);
	}
	if (j > m) {
		//std::cout << "oups" << std::endl;
		return(result);
	}
	arma::rowvec indicesi = arma::linspace<arma::rowvec>(1, m - j + 1, m - j + 1);
	for (unsigned int i = 0; i < indicesi.size(); ++i) {
		int sizeej = indicesi(i);
		arma::rowvec indicesj = arma::linspace<arma::rowvec>(1, m - sizeej + 1, m - sizeej + 1);
		for (unsigned int j = 0; j < indicesj.size(); ++j) {
			int binf = indicesj(j);
			int bsup = binf + sizeej - 1;
			arma::umat rowvector;
			rowvector << binf << bsup << arma::endr;
			result = arma::join_vert(result, rowvector);
		}
	}
	return(result);
}
arma::umat allej(int j, int m);

// [[Rcpp::export]]
double pejp1_yjej(arma::urowvec ejp1, int yj, arma::urowvec ej, int mu, double p) {
	double proba = 0;
	arma::rowvec ejminus;
	ejminus << ej(0) << yj - 1;
	arma::rowvec ejequal;
	ejequal << yj << yj;
	arma::rowvec ejplus;
	ejplus << yj + 1 << ej(1);

	//pejp1_yjejzj0
	double pejp1_yjejzj0 = 0;
	if (compare_vec(ejp1, ejminus) || compare_vec(ejp1, ejequal) || compare_vec(ejp1, ejplus)) {
		pejp1_yjejzj0 = (double)(ejp1(1) - ejp1(0) + 1) / (ej(1)- ej(0) + 1);
	}
	//pejp1_yjejzj1
	double dmuejminus = 0;
	if (ejminus(0) > ejminus(1)) {
		dmuejminus = inf;
	}
	else {
		arma::rowvec ejminusbis;
		ejminusbis = ejminus;
		ejminusbis.for_each([mu](arma::rowvec::elem_type& val) {
			val -= mu; 
			val = std::abs(val);
		});
		dmuejminus = arma::min(ejminusbis);
	}
	double dmuejplus = 0;
	if (ejplus(0) > ejplus(1)) {
		dmuejplus = inf;
	}
	else {
		arma::rowvec ejplusbis;
		ejplusbis = ejplus;
		ejplusbis.for_each([mu](arma::rowvec::elem_type& val) {
			val -= mu; 
			val = std::abs(val);
		});
		dmuejplus = arma::min(ejplusbis);
	}

	arma::rowvec ejequalbis;
	ejequalbis = ejequal;
	ejequalbis.for_each([mu](arma::rowvec::elem_type& val) {
		val -= mu; val = std::abs(val);
	});
	double dmuejequal = arma::min(ejequalbis);

	arma::rowvec ejp1bis(ejp1.n_elem);
	for (int in = 0; in < ejp1.n_elem; in++) {
		ejp1bis(in) = ejp1(in);
	}
	ejp1bis.for_each([mu](arma::rowvec::elem_type& val) {
		val -= mu; val = std::abs(val);
	});
	double dmuejp1 = arma::min(ejp1bis);

	arma::rowvec all_dmu;
	all_dmu << dmuejminus << dmuejequal << dmuejplus;
	int pejp1_yjejzj1 = 0;
	if ((dmuejp1 == arma::min(all_dmu)) && ((compare_vec(ejp1, ejminus) || compare_vec(ejp1, ejequal) || compare_vec(ejp1, ejplus)))) {
		pejp1_yjejzj1 = 1;
	}
	else {
		pejp1_yjejzj1 = 0;
	}
	proba = p * pejp1_yjejzj1 + (1 - p) * pejp1_yjejzj0;
	return(proba);
}
double pejp1_yjej(arma::urowvec ejp1, int yj, arma::urowvec ej, int mu, double p);

// [[Rcpp::export]]
double pejp1zj1_yjej(arma::urowvec ejp1, unsigned int yj, arma::urowvec ej, int mu, double p) {
	double proba = 0;
	arma::rowvec ejminus;
	ejminus << ej(0) << yj - 1;
	arma::rowvec ejequal;
	ejequal << yj << yj;
	arma::rowvec ejplus;
	ejplus << yj + 1 << ej(1);

	double dmuejminus = 0;
	double dmuejplus = 0;
	if (ejminus(0) > ejminus(1)) {
		dmuejminus = inf;
	}
	else {
		arma::rowvec ejminusbis;
		ejminusbis = ejminus;
		ejminusbis.for_each([mu](arma::rowvec::elem_type& val) {
			val -= mu; val = std::abs(val);
		});
		dmuejminus = arma::min(ejminusbis);
	}
	if (ejplus(0) > ejplus(1)) {
		dmuejplus = inf;
	}
	else {
		arma::rowvec ejplusbis;
		ejplusbis = ejplus;
		ejplusbis.for_each([mu](arma::rowvec::elem_type& val) {
			val -= mu; val = std::abs(val);
		});
		dmuejplus = arma::min(ejplusbis);
	}

	arma::rowvec ejequalbis;
	ejequalbis = ejequal;
	ejequalbis.for_each([mu](arma::rowvec::elem_type& val) {
		val -= mu; val = std::abs(val);
	});
	double dmuejequal = arma::min(ejequalbis);

	arma::rowvec ejp1bis(ejp1.n_elem);
	for (int in = 0; in < ejp1.n_elem; in++) {
		ejp1bis(in) = (int) ejp1(in);
	}
	ejp1bis.for_each([mu](arma::rowvec::elem_type& val) {
		val -= mu; val = std::abs(val);
	});
	double dmuejp1 = arma::min(ejp1bis);

	arma::rowvec all_dmu;
	all_dmu << dmuejminus << dmuejequal << dmuejplus;
	int pejp1_yjejzj1 = 0;
	if (dmuejp1 == arma::min(all_dmu) && (compare_vec(ejp1, ejminus) || compare_vec(ejp1, ejequal) || compare_vec(ejp1, ejplus))) {
		pejp1_yjejzj1 = 1;
	}
	else {
		pejp1_yjejzj1 = 0;
	}
	proba = p * pejp1_yjejzj1;
	return(proba);
}
double pejp1zj1_yjej(arma::urowvec ejp1, unsigned int yj, arma::urowvec ej, int mu, double p);

// [[Rcpp::export]]
double pejp1zj1_ej(arma::urowvec ejp1, arma::urowvec ej, int mu, double p) {
	double proba = 0;
	arma::uvec ej_bounds = linspace<uvec>(ej(0), ej(1), ej(1) - ej(0) + 1);
	for (unsigned int i = 0; i < ej_bounds.size(); ++i) {
		int yj = ej_bounds(i);
		proba += pejp1zj1_yjej(ejp1, yj, ej, mu, p);
	}
	proba = (double)proba / (ej(1) - ej(0) + 1);
	return(proba);
}
double pejp1zj1_ej(arma::urowvec ejp1, arma::urowvec ej, int mu, double p);

// [[Rcpp::export]]
double pyj_ej(unsigned int yj, arma::urowvec ej) {
	double proba = 0;
	if (ej(0) <= yj && yj <= ej(1)) {
		proba = 1.0 / (ej(1) - ej(0) + 1);
	}
	else {
		proba = 0;
	}
	return(proba);
}
double pyj_ej(unsigned int yj, arma::urowvec ej);

// [[Rcpp::export]]
double pejp1_ej(arma::urowvec ejp1, arma::urowvec ej, int mu, double p) {  
	double proba = 0;
	arma::uvec allyj;
	// we suppose ejp1 included in ej (not checked here)
	if (ejp1(1) == ejp1(0)) { // |ejp1|=1
		if (ejp1(1) < ej(1) && ejp1(0) > ej(0)) { //ejp1 doesn't touch any bound
			allyj << ejp1(0);
		}
		else {
			if (ejp1(1) < ej(1)) { // ejp1 doesn't touch right bound
				allyj << ejp1(0) << ejp1(0) + 1;
			}
			else { // ejp1 doesn't touch left bound
				allyj << ejp1(0) - 1 << ejp1(0);
			}
		}
	}
	else { // |ejp1|>1
		if (ejp1(1) < ej(1)) {// ejp1 doesn't touch right bound
			allyj << ejp1(1) + 1;
		}
		else { // ejp1 doesn't touch left bound
			allyj << ejp1(0) - 1;
		}
	}
	for (unsigned int i = 0; i < allyj.size(); ++i){
		unsigned int yj = allyj(i);
		proba += pejp1_yjej(ejp1, yj, ej, mu, p) * pyj_ej(yj, ej);
	}
	return(proba);
}
double pejp1_ej(arma::urowvec ejp1, arma::urowvec ej, int mu, double p);


// [[Rcpp::export]]
double pej(arma::urowvec& ej,
	int j, int m, int mu, double p,
	arma::colvec& z1tozjm1)
{

	if (j == 1) {
		return(1);
	}

	if (ej.size() == 1) {
		unsigned int ejn = ej(0);
		ej = arma::urowvec(2);
		ej = ej.ones() * ejn;
	}

	arma::colvec z1tozjm2 = arma::vectorise(z1tozjm1, 0);
	double zjm1 = z1tozjm2(z1tozjm2.n_rows - 1);
	z1tozjm2.shed_row(z1tozjm2.n_rows-1);

	double proba = 0;
	if (zjm1) { // zjm1 is known
		arma::umat tabint = allej(j - 1, m); // may have a problem here
		int nbtabint = tabint.n_rows;
		for (int i = 0; i < nbtabint; ++i) {
			arma::urowvec ejm1 = tabint.row(i);
			if ((ej(0) >= ejm1(0)) && (ej(1) <= ejm1(1))) { //to accelerate, check if ejm is included in ej
				proba += pejp1zj1_ej(ej, ejm1, mu, p) * pej(ejm1, j - 1, m, mu, p, z1tozjm2);
			}
		}	
	}
	
	else { // zjm1 is unknown
		arma::umat tabint = allej(j - 1, m);
		int nbtabint = tabint.n_rows;
		for (int i = 0; i < nbtabint; ++i) {
			arma::urowvec ejm1 = tabint.row(i);
			if ((ej(0) >= ejm1(0)) && (ej(1) <= ejm1(1))) { //to accelerate, check if ejm is included in ej
				proba += pejp1_ej(ej, ejm1, mu, p) * pej(ejm1, j - 1, m, mu, p, z1tozjm2);
			}
		}
	}
	//std::cout << proba << std::endl;
	return(proba);
}

double pej(arma::urowvec& ej,
	int j, int m, int mu, double p,
	arma::colvec& z1tozjm1);

// [[Rcpp::export]]
List ordiemCpp(int m,
	const arma::cube& tab_pej, // don't forget to pass it 
	const arma::colvec& x,
	const arma::colvec& tabmu0,
	const arma::colvec& tabp0,
	double eps = 1,
	int iter_max = 100)
{

	double ml = -inf;
	double p_ml = tabp0(0);
	int mu_ml = tabmu0(0);

	int n = x.n_rows;
	arma::vec w = arma::vec(n); // in code, it is as a paramter, but never used as paramteer
	w.ones();
	

	double ntot = arma::sum(w);

	int a = 0;

	for (unsigned int imu = 0; imu < tabmu0.n_rows; ++imu) {
		int mu = tabmu0(imu);
		double mlold = -inf;
		for (unsigned int ip = 0; ip < tabp0.n_rows; ++ip) {
			double p = tabp0(ip);
			int nostop = 1;
			int iter = 0;
			while (nostop) {
				++iter;
				// -- E step ---
				// first: compute px for each modality
				arma::colvec pallx = arma::colvec(m);
				pallx.zeros();

				arma::colvec px = arma::colvec(n);
				px.zeros();

				for (int i = 0; i<m; ++i) {
					arma::colvec tmp = arma::colvec(m);
					tmp = tmp.ones() * p;
					for (int deg = 0; deg < m; ++deg) {
						tmp(deg) = std::pow(tmp(deg), deg);
					}
					//Q.subcube( first_row, first_col, first_slice, last_row, last_col, last_slice )
					arma::rowvec subcube = tab_pej.subcube(i, imu, 0, i, imu, tab_pej.n_slices - 1);
					
					arma::vec vecprod = trans(subcube) % tmp;
					pallx(i) = arma::sum(vecprod);

					arma::uvec ids = arma::find(x == (i+1));

					px.elem(ids).fill(std::max(1e-300,pallx(i)));
				}
				
				arma::vec vecprod = w % arma::log(px);

				double mlnew =  arma::sum(vecprod);
				
				// first: compute pxz1 for each modality
				arma::mat pallxz1 = arma::mat(m,m-1);
				pallxz1.zeros();


				for(int i=0; i < m; ++i){
					for(int j=0; j < (m-1); ++j){
						arma::colvec z1tozmm1 = arma::colvec(m-1);
						z1tozmm1.zeros();
						z1tozmm1(j) = 1;
						arma::urowvec ivec; // to gain time : put it out of this loop
						ivec << (i+1);

						double pejv = pej(ivec,m,m,mu,p,z1tozmm1); // TODO change that
						pallxz1(i,j) = (double)(std::max(1e-300, pejv))/std::max(1e-300,pallx(i));
					}
				}

				//second: affect each modality value to the corresponding units
				arma::mat pxz1 = arma::mat(n,m-1);
				pxz1.zeros(); 

				for(int i=0; i < m; ++i){
					arma::uvec whereisi = arma::find(x == (i+1));
					arma::uvec all_cols_pxz1 = arma::linspace<arma::uvec>(0, pxz1.n_cols-1, pxz1.n_cols);

					int sumwhereisi = whereisi.size();
					arma::colvec matsumwhereisi = arma::colvec(sumwhereisi);
					matsumwhereisi.ones();

					
					arma::rowvec pallxz1i = pallxz1.row(i);
					arma::mat submat = matsumwhereisi * pallxz1i;
					pxz1.submat(whereisi, all_cols_pxz1) = submat; // maybe a problem here

				}
				arma::mat temp1 = arma::mat(1,m-1);
				temp1.ones();
				arma::mat temp2 = w * temp1;
				pxz1 = temp2 % pxz1;

				

				// ---- M step ----
				double sum = arma::accu(pxz1);
				double pmean = (double)sum/(ntot*(m-1));
				//std::cout << pmean << std::endl;
				if(!(mlnew==-inf)){
					if((std::abs(mlnew-mlold)/ntot < eps) || (iter>(iter_max-1))){
						nostop = 0; 
						if (mlnew > ml) {
							ml = mlnew;
							p_ml = pmean;
							mu_ml = mu;
						}
					}
				}
				else{
					if(iter>(iter_max-1)){
						nostop = 0;
						if (mlnew > ml)  {
							ml = mlnew;
							p_ml = pmean;
							mu_ml = mu ;
						}
					}
				}
				mlold = mlnew;

			}// while
		} // p
		++a;
	}//mu

	List result = List::create(mu_ml, p_ml, ml);
	/*std::string result = "ok";
	std::cout << mu_ml << "  " << p_ml << "  " << ml << endl;*/
	return(result);
}
