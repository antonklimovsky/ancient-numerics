//--------------------------------------------------
// YASM - Yet Another Simplex Method implementation
// ver. 0.1, May 27, 2001
//
// Author: A. Klimovsky, root@ludus.kharkiv.com
// All rights reserved
//--------------------------------------------------
#include <math.h>
#include <fstream.h>
#include <stdlib.h>
#include <string>

#include "LinAlg.h"
#include "Rational.h"
#include "Log.h"

//-------------------------------------------
// Global data
//-------------------------------------------

//--- parser data

string stringValue;
int numberValue;
typedef enum {
    END, NUMBER, STRING, SLASH, OPENBRACE, CLOSEBRACE, ENDL
} lexState;
lexState state;
bool needGet;
int curLine;
ifstream in;

//--- parsed data

int n = 9;
int m = 9;
Matrix<FractInt> a;

//-----------------------------------
// Simple lexic analizer routines
//-----------------------------------

istream& myGet(char& ch)
{
    return in.get(ch);
}

lexState getToken()
{
    char ch = 0;

    if (!needGet) {
        needGet = true;
        return state;
    }
onceMore:
    do {
        if (!myGet(ch)) return END;
    } while (isspace(ch) && ch != '\n');

    switch (ch) {
        case '\n':
            curLine++;
            return ENDL;
        case '#': // one-line comments
            do {
                if (!myGet(ch)) return END;
            } while (ch != '\n');
            curLine++;
            goto onceMore;
        case '/':
            return SLASH;
        case '[':
            return OPENBRACE;
        case ']':
            return CLOSEBRACE;
        default:
            if (isdigit(ch) || ch == '-') {
                in.putback(ch);
                in >> numberValue;
                return NUMBER;
            }
            if (isalpha(ch)) {
                stringValue = ch;
                while (myGet(ch) && isalnum(ch))
                    stringValue.push_back(ch);
                in.putback(ch);
                return STRING;
            }
    }
}

void reportError(char* text)
{
    cout << "[Error] on line " << curLine << ": " << text;
    exit(1);
}

//----------------------------------------
// Second pass routines
// The main purpose of the second pass is
// to parse input data and fill
// necessary internal repesentation of the
// task.
//----------------------------------------

int secondPassGetVector(Vector<FractInt>& v)
{
    int i = 0;
    int nom;
    int denom;

    if ((state = getToken()) != ENDL)
        reportError("a new line expected");
    do {
        state = getToken();
        if (state != NUMBER) {
            needGet = false;
            break;
        }
        nom = numberValue;
        denom = 1;
        state = getToken();
        if (state == SLASH) {
            state = getToken();
            if (state != NUMBER)
                reportError("denominator expected");
            denom = numberValue;
            if (denom == 0)
                reportError("meaningless division by zero");
        } else
            needGet = false;
        v[i] = FractInt(nom, denom);
        i++;
    } while (true);
    return i;
}

void secondPassA()
{
    int oldN;
    int newN;

    oldN = secondPassGetVector(a[0]);
    if (oldN == 0)
        reportError("a matrix row expected");
    m = 1;
    while ((newN = secondPassGetVector(a[m])) > 0)
        if (oldN != newN)
            reportError("matrix should be rectangular");
        else
            m++;
    n = oldN;
}
void secondPassRoot()
{
    if ((state = getToken()) != OPENBRACE)
        reportError("'[' expected");
    if ((state = getToken()) != STRING || stringValue != "A")
        reportError("'A' expected");
    if ((state = getToken()) != CLOSEBRACE)
        reportError("']' expected");
    secondPassA();
}

void secondPass(char* fname)
{
    in.open(fname);
    if (!in) {
        cout << "Fatal: can't open " << fname;
        exit(1);
    }
    curLine = 1;
    needGet = true;
    secondPassRoot();
    in.close();
}

void main(int argc, char* argv[])
{
    Log<FractInt> log;

    cout << "YASM - Yet Another Simplex Method implementation" << "\n";
    cout << "ver. 0.1 (c) 2001 A.Klimovsky" << "\n\n";

    if (argc != 2) {
        cout << "Usage: yasm InputFile" << "\n";
        return;
    }

    a = Matrix<FractInt>(m, n);

    secondPass(argv[1]);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            a[j][i] += 1;
    log(a);
    log("(A+7)^t\n");
    log(a.transpose());
}