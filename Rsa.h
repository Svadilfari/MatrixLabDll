#pragma once
#include "Matrixclasses.h"
#include "RsaProblems.h"
#define EPS 0.00000001

#ifdef MATRIXDLL_EXPORTS
#define RSA_API __declspec(dllexport)
#else
#define RSA_API __declspec(dllimport)
#endif

class RSA_API RSA {
protected:
	Matrix mX;
	Matrix m_scores;
	Matrix m_loadings;
	Matrix mE;

public:
	Matrix center(const Matrix& A) const;
	Matrix scaling(const Matrix& A, bool centered = true) const;

	vector<Matrix> NIPALS(int PC = -1) const;
	RSA(const Matrix& X, int PC = -1);

	Matrix get_X() const;
	Matrix get_scores() const;
	Matrix get_loadings() const;

};

class RSA_API Residuals : public RSA {
protected:
	vector<double> m_lambda;
	int mPC;

public:
	Residuals(const Matrix& A, int PC);
	Residuals(const RSA& rsa, int PC);

	vector<double> leverages() const;
	double mean_leverage() const;
	vector<double> variances() const;
	double mean_variance() const;

	double TRV() const;
	double ERV() const;

	Matrix get_E() const;
	int get_A() const;
};