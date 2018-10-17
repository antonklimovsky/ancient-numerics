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
#include "SimplexMethod.h"
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

int n;
int m;
Matrix<FractInt> a;
Vector<FractInt> b;
Vector<FractInt> c;

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

//---------------------------------------
// First pass routines
// The main purpose of the first pass is
// to obtaine a priori unknown task
// dimensions
//---------------------------------------

int firstPassGetVector()
{
    int i = 0;

    if ((state = getToken()) != ENDL)
        reportError("a new line expected");
    do {
        state = getToken();
        if (state != NUMBER) {
            needGet = false;
            break;
        }
        state = getToken();
        if (state == SLASH) {
            state = getToken();
            if (state != NUMBER)
                reportError("denominator expected");
        } else
            needGet = false;
        i++;
    } while (true);
    return i;
}

void firstPassA()
{
    int oldN;
    int newN;

    oldN = firstPassGetVector();
    if (oldN == 0)
        reportError("a matrix row expected");
    m = 1;
    while ((newN = firstPassGetVector()) > 0)
        if (oldN != newN)
            reportError("matrix should be rectangular");
        else
            m++;
    n = oldN;
}

void firstPassB()
{
    int m;

    m = firstPassGetVector();
    if (m != ::m)
        reportError("dimensions of A and b are inconsistent");
    if ((state = getToken()) != ENDL)
        reportError("new line expected");
}

void firstPassC()
{
    int n;

    n = firstPassGetVector();
    if (n != ::n)
        reportError("dimensions of A and c are inconsistent");
    state = getToken();
    if (state != ENDL && state != END)
        reportError("new line or EOF expected");
}

void firstPassRoot()
{
    if ((state = getToken()) != OPENBRACE)
        reportError("'[' expected");
    if ((state = getToken()) != STRING || stringValue != "A")
        reportError("'A' expected");
    if ((state = getToken()) != CLOSEBRACE)
        reportError("']' expected");
    firstPassA();

    if ((state = getToken()) != OPENBRACE)
        reportError("'[' expected");
    if ((state = getToken()) != STRING || stringValue != "b")
        reportError("'b' expected");
    if ((state = getToken()) != CLOSEBRACE)
        reportError("']' expected");
    firstPassB();

    if ((state = getToken()) != OPENBRACE)
        reportError("'[' expected");
    if ((state = getToken()) != STRING || stringValue != "c")
        reportError("'c' expected");
    if ((state = getToken()) != CLOSEBRACE)
        reportError("']' expected");
    firstPassC();
}

void firstPass(char* fname)
{
    in.open(fname);
    if (!in) {
        cout << "Fatal: can't open " << fname;
        exit(1);
    }
    curLine = 1;
    needGet = true;
    firstPassRoot();
    in.close();
}

//----------------------------------------
// Second pass routines
// The main purpose of the second pass is
// to parse the input data and to fill
// the internal repesentation of the
// task.
//----------------------------------------

void secondPassGetVector(Vector<FractInt>& v)
{
    int j;
    int nom;
    int denom;
    int i;

    state = getToken();

    for (i = 0; i < v.get_size(); i++) {
        state = getToken();
        nom = numberValue;
        denom = 1;
        state = getToken();
        if (state == SLASH) {
            state = getToken();
            denom = numberValue;
        } else
            needGet = false;
        v[i] = FractInt(nom, denom);
    }
}

void secondPassA()
{
    int j;

    for (j = 0; j < m; j++)
        secondPassGetVector(a[j]);
    state = getToken();     
}

void secondPassB()
{
    int m;

    secondPassGetVector(b);
    state = getToken();
}

void secondPassC()
{
    secondPassGetVector(c);
    state = getToken(); 
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

    if ((state = getToken()) != OPENBRACE)
        reportError("'[' expected");
    if ((state = getToken()) != STRING || stringValue != "b")
        reportError("'b' expected");
    if ((state = getToken()) != CLOSEBRACE)
        reportError("']' expected");
    secondPassB();

    if ((state = getToken()) != OPENBRACE)
        reportError("'[' expected");
    if ((state = getToken()) != STRING || stringValue != "c")
        reportError("'c' expected");
    if ((state = getToken()) != CLOSEBRACE)
        reportError("']' expected");
    secondPassC();
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

int main(int argc, char* argv[])
{
    Log<FractInt> log;
    SimplexMethod<FractInt>* s;

    cout << "YASM - Yet Another Simplex Method implementation" << "\n";
    cout << "ver. 0.1 (c) 2001 A.Klimovsky" << "\n\n";

    if (argc != 2) {
        cout << "Usage: yasm InputFile" << "\n";
        return 1;
    }

    firstPass(argv[1]);

    a = Matrix<FractInt>(m, n);
    b = Vector<FractInt>(m);
    c = Vector<FractInt>(n);

    secondPass(argv[1]);

    s = new SimplexMethod<FractInt>(a, b, c, log);
    
    s->go();

    delete s;
}
