#pragma once
#include<iostream>
#include<exception>

#ifdef MATRIXDLL_EXPORTS
#define RSAPROBLEMS_API __declspec(dllexport)
#else
#define RSAPROBLEMS_API __declspec(dllimport)
#endif

class RSAPROBLEMS_API ProblemRSA {
public:
    virtual char const* what() const noexcept = 0;
};

class RSAPROBLEMS_API ZeroDiv : public ProblemRSA {
public:
    char const* what() const noexcept override
    {
        return "Error: Division by zero has occured!";
    }
};

class RSAPROBLEMS_API IncorrectPC : public ProblemRSA {
public:
    char const* what() const noexcept override
    {
        return "Error: Number of eigen values is incorrect!";
    }
};