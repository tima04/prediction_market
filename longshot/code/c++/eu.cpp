#include <iostream>
#include <armadillo>
#include <map>
#include <vector>
#include <set>
//#include <math>

/*
head(df) ==>
raceid winner numhorses winodds
  1      0        10    1.25
 */


using namespace std;
using namespace arma;

double cara(double x, double theta) {
  return (1 - exp(-theta*x))/theta;
}

vec solve_ps(double t, vec Rs, double epsilon=1e-10) {
  int n  =  Rs.n_rows; 
  vec utils(n);
  double u_min_1 = cara(-1, t);

  double sum = 0;
  for (int i = 0; i < n; i++) {
	utils(i) = cara(Rs(i) + 1, t); 
	sum += 1/(utils(i) - u_min_1);
	  }
  double w = 1/sum + u_min_1;

  vec ps(n);
  for (int i = 0; i < n; i++) 
    ps(i) = (w - u_min_1)/(utils(i) - u_min_1);
  return (ps);
}

double loglike(double t, map<int, mat> df) { 
  double rslt = 0;
  for (auto it = df.begin(); it != df.end(); ++it) {
	vec Rs = it->second.col(3);
	vec winner = it->second.col(1);
	vec ps = solve_ps(t, Rs);

	rslt += log(dot(ps, it->second.col(1)));

	double foo = log(dot(ps, it->second.col(1)));
	if (foo < -10e10) {
	  cout << it->second << endl << ps << endl;
	}
	
  }
  return rslt;
}

map<int, mat> make_df(string fname) {
  mat X;
  clock_t t0 = clock();
  X.load(fname);
  X = X.rows(1, X.n_rows-1); //first row is all zero, header is converted to 0s
  set<int> raceid;
  map<int, mat> df;

  for (int i = 0; i < X.n_rows; i++) {
	int id = X(i, 0); 
	raceid.insert(id);
  }

  for (auto it=raceid.begin(); it!=raceid.end(); ++it) {
	arma::uvec ids = find(X.col(0) == *it); // Find indices
	df[*it] = X.rows(ids);
  }
  return df;
}

// assumes df is sorted by raceid
map<int, mat> make_df2(string fname) {
  mat X;
  clock_t t0 = clock();
  X.load(fname);
  X = X.rows(1, X.n_rows-1); //first row is all zero, (header is converted to 0s)
  // cout << "time to load: " << clock() - t0 << endl;

  X.save("foo.txt", csv_ascii);
  set<int> raceid;
  map<int, mat> df;

  int old_id = X(0, 0); 
  int begin_race = 0, end_race = 0;

  for (int i = 0; i < X.n_rows; i++) {

	int id =  X(i, 0); 

	if (id != old_id) {
	  //debug: for the time being get rid of problematic races
	  // most likely rounding of raceid is creating problem.
	  mat temp = X.rows(begin_race, i-1);
	  if(sum(temp.col(1)) > 0.99) 
		df[old_id] = temp; //X.rows(begin_race, i-1);

	  old_id = id;
	  begin_race = i; 
	}

	//last race
	if (i == X.n_rows-1)
	  df[id] = X.rows(begin_race, i);
  }
  return df;
}

int main() {
  map<int, mat> df = make_df2("../../data/dataset_full.csv");
  cout << loglike(0.5, df) << endl;
  return 0;
}
