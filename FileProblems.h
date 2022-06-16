#pragma once
#include<iostream>
#include<exception>

#ifdef MATRIXDLL_EXPORTS
#define FILEPROBLEMS_API __declspec(dllexport)
#else
#define FILEPROBLEMS_API __declspec(dllimport)
#endif


class FILEPROBLEMS_API FileProblem {
public:
    virtual char const* what() const noexcept = 0;
};

class FILEPROBLEMS_API OpenProblem : public FileProblem {
public:
    char const* what() const noexcept override
    {
        return "Error: Cannot open the file!";
    }
};

class FILEPROBLEMS_API FormatProblem : public FileProblem {
public:
    char const* what() const noexcept override
    {
        return "Error: Format of file is incorrect!";
    }
};

class FILEPROBLEMS_API UnknownExtention : public FileProblem {
public:
    char const* what() const noexcept override
    {
        return "Error: Extention of file is incorrect!";
    }
};
