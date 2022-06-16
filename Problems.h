#pragma once
#include<iostream>
#include<exception>

#ifdef MATRIXDLL_EXPORTS
#define PROBLEMS_API __declspec(dllexport)
#else
#define PROBLEMS_API __declspec(dllimport)
#endif

class PROBLEMS_API MatrixProblems {
public:
    virtual char const* what() const noexcept = 0;
};

class PROBLEMS_API ZeroMN : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Matrix cannot have 0 rows or 0 columns!";
    }
};

class PROBLEMS_API NotSquare : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Matrix is not a square matrix!";
    }
};

class PROBLEMS_API NoInverse : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Inverse matrix doesn't exist!";
    }
};

class PROBLEMS_API ConcatProblem : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Matrix should have same amount of rows for concatenation!";
    }
};

class SizeProblem : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Matrix should have same sizes!";
    }
};

class PROBLEMS_API SizeProblemMult : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Number of columns of the 1st matrix should be equal to number of rows of the 2nd!";
    }
};

class PROBLEMS_API NotAVector : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Operand should be a vector!";
    }
};

class PROBLEMS_API RowsAndCols : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: All rows' sizes shoud equal M. All columns' sizes shoud equal N.";
    }
};

class PROBLEMS_API IndexOutOfRange : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Index is out of range!";
    }
};

class PROBLEMS_API ZeroDivision : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Division by zero!";
    }
};

