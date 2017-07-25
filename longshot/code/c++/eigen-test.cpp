#include <iostream>
#include <Eigen/Dense>
//using Eigen::MatrixXd;
using namespace Eigen;
using namespace std;

/*
void foo(int nrows, int ncols) {
	MatrixXf X = MatrixXf::Zero(nrows, ncols);
	ifstream fin ("../../data/dataset_full.csv");

	if (fin.is_open())
	{
		for (int row = 0; row < nrows; row++)
			for (int col = 0; col < ncols; col++)
			{
				float item = 0.0;
				fin >> item;
				X(row, col) = item;
			}
		fin.close();
	}
	cout << "X = " << endl << X << endl;
}
*/


//# g++ -I ~/.local/lib/eigen-eigen-67e894c6cd8f eigen-test.cpp -o my_program                                                         

#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define MAXBUFSIZE  ((int) 1e6)

MatrixXd readMatrix(const char *filename)
    {

    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];


    // Read numbers from file into buffer.
    ifstream infile;
    infile.open(filename);
    while (! infile.eof())
        {
        string line;
        getline(infile, line);

        int temp_cols = 0;
        stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
        }

    infile.close();

    rows--;

    // Populate matrix with numbers.
    MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    return result;
    };


int main()
{
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;
  char *filename = "../../data/test_data.csv";
  std::cout << "hello" << endl;
  MatrixXd foo = readMatrix(filename);
}
