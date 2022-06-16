#include "pch.h"
#include "RSA.h"

Residuals::Residuals(const Matrix& A, int PC) : RSA(A, PC) {
	Matrix L = m_scores.Transpose() * m_scores;
	for (int i = 0; i < PC; ++i) {
		m_lambda.push_back(L.get_matrix()[i][i]);
	}
	mPC = PC;
}

Residuals::Residuals(const RSA& rsa, int PC) : RSA(rsa) {
	vector<Matrix> matrices = NIPALS(PC);
	m_scores = matrices[0];
	Matrix L = m_scores.Transpose() * m_scores;
	for (int i = 0; i < PC; ++i) {
		m_lambda.push_back(L.get_matrix()[i][i]);
	}
	m_loadings = matrices[1];
	mE = matrices[2];
	mPC = PC;
}

//Анализ остатков
vector<double> Residuals::leverages() const {
	vector<double> h;
	vector<vector<double>> matrix = m_scores.get_matrix();
	for (vector<double> row : matrix) {
		double h_i = 0;
		for (int i = 0; i < mPC; ++i) h_i += (row[i] * row[i]) / m_lambda[i];
		h.push_back(h_i);
	}
	return h;
}

double Residuals::mean_leverage() const {
	vector<double> levs = leverages();
	double h0 = 0;
	for (double h_i : levs) h0 += h_i;
	h0 = h0 / levs.size();
	return h0;
}

vector<double> Residuals::variances() const {
	vector<double> v;
	vector<vector<double>> matrix = mE.get_matrix();
	for (vector<double> row : matrix) {
		double v_i = 0;
		for (double el : row) v_i += el * el;
		v.push_back(v_i);
	}
	return v;
}

double Residuals::mean_variance() const {
	vector<double> v = variances();
	double v0 = 0;
	for (double v_i : v) v0 += v_i;
	v0 = v0 / v.size();
	return v0;
}

double Residuals::TRV() const {
	int J = mE.get_columns();
	double trv = mean_variance() / J;
	return trv;
}

double Residuals::ERV() const {
	int I = mE.get_rows();
	double v0 = mean_variance();
	Matrix D = mX;
	D = center(D);
	D = scaling(D);
	double sum_x = 0;
	vector<vector<double>> matrix = D.get_matrix();
	for (vector<double> row : matrix) {
		for (double el : row) sum_x += el * el;
	}
	if (fabs(sum_x) < ACCURACY) throw ZeroDiv();

	double erv = 1 - (I * v0) / (sum_x);
	return erv;
}

//Аксессоры
Matrix Residuals::get_E() const { return mE; };
int Residuals::get_A() const { return mPC; };
