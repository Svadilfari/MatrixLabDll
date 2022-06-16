#include "pch.h"
#include "RSA.h"

//Предобработка данных
Matrix RSA::center(const Matrix& A) const {
	Matrix B = A.Transpose();
	vector<vector<double>> matrix = B.get_matrix();
	for (auto& row : matrix) {
		double m = 0;
		for (double element : row) m += element;
		m = m / B.get_columns();
		for (double& element : row) element -= m;
	}
	Matrix C(matrix);
	return C.Transpose();
}

Matrix RSA::scaling(const Matrix& A, bool centered) const {
	Matrix B = A.Transpose();
	vector<vector<double>> matrix = B.get_matrix();
	for (auto& row : matrix) {
		if (centered) {
			double s = 0;
			for (double element : row) s += element * element;
			s = s / (B.get_columns() - 1.);
			s = pow(s, 0.5);
			if (fabs(s) < ACCURACY) continue;
			for (double& element : row) {
				element = element / s;
			}
		}
		else {
			double m = 0;
			for (double element : row) m += element;
			m = m / B.get_columns();
			double s = 0;
			for (double element : row) s += (element - m) * (element - m);
			s = s / (B.get_columns() - 1.);
			if (fabs(s) < ACCURACY) continue;
			for (double& element : row) element = element / s;
		}
	}
	Matrix C(matrix);
	return C.Transpose();
}

//Алгоритм NIPALS
vector<Matrix> RSA::NIPALS(int PC) const {
	int min_ij = mX.get_columns() < mX.get_rows() ? mX.get_columns() : mX.get_rows();
	if (PC > min_ij) throw IncorrectPC();
	if (PC < 0) PC = min_ij;
	vector<vector<double>> matrix_p;
	vector<vector<double>> matrix_t;

	Matrix D = mX;
	D = center(D);
	D = scaling(D);
	Matrix E = D;

	for (size_t h = 0; h < PC; ++h) {
		Matrix t = E.get_column_j(h);
		Matrix d;
		Matrix p;
		int a = 0;

		do {
			p = (t.Transpose() * E) / (scalarMult(t, t));
			p = p.Transpose();
			p = p / p.Norm();
			//cout << p << endl << endl;
			Matrix t_old = t;
			t = (E * p) / scalarMult(p, p);
			//cout << t << endl << endl;
			d = t_old - t;
			a++;
		} while (d.Norm() > EPS);
		//cout << a << endl;

		E = E - (t * p.Transpose());
		vector<double> p_col = ((p.Transpose()).get_matrix())[0];
		matrix_p.push_back(p_col);
		vector<double> t_col = ((t.Transpose()).get_matrix())[0];
		matrix_t.push_back(t_col);
	}
	Matrix P = Matrix(matrix_p).Transpose();
	Matrix T = Matrix(matrix_t).Transpose();
	vector<Matrix> matrices = { T, P, E };
	return matrices;
}

RSA::RSA(const Matrix& X, int PC) : mX(X) {
	vector<Matrix> matrices = NIPALS(PC);
	m_scores = matrices[0];
	m_loadings = matrices[1];
	mE = matrices[2];
}

//аксессоры
Matrix RSA::get_X() const { return mX; };
Matrix RSA::get_scores() const { return m_scores; };
Matrix RSA::get_loadings() const { return m_loadings; };
