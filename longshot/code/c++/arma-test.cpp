#include <iostream>
#include <armadillo>
#include <map>
#include <vector>
#include <set>

using namespace std;
using namespace arma;

int main()
  {
  mat A = randu<mat>(4,5);
  mat B = randu<mat>(4,5);
  //cout << A*B.t() << endl;

  X.load("../../data/test_data.csv");
  X = X.rows(1, X.n_rows-1);
  set<int> raceid;
  map<int, mat> df;
  for (int i = 0; i < X.n_rows; i++) {
	int id = X(i, 0); 
	//cout << id << endl;
	raceid.insert(id);
  }
  for (auto it=raceid.begin(); it!=raceid.end(); ++it) {
	arma::uvec ids = find(X.col(0) == *it); // Find indices
	//cout << ids << endl;
	df[*it] = X.rows(ids);
  }
  cout << X.rows(0,3) << endl;
  return 0;
  }

