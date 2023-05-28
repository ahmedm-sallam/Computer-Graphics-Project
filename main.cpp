/*
صلي علي سيدنا محمد
   _____              _        _                   __  __
  / ____|     /\     | |      | |          /\     |  \/  |
 | (___      /  \    | |      | |         /  \    | \  / |
  \___ \    / /\ \   | |      | |        / /\ \   | |\/| |
  ____) |  / ____ \  | |____  | |____   / ____ \  | |  | |
 |_____/  /_/    \_\ |______| |______| /_/    \_\ |_|  |_|
*/
#include <tchar.h>
#include <windows.h>
#include<vector>
#include <iostream>
#include<cmath>
#include<fstream>
#include<sstream>
#include <list>
#include <stack>

#define MAXENTRIES 600

using namespace std; //to make combination between console and my window.

int Num_of_Points = 0;
int counter = 0;

struct Vertex {
    int x, y;

    explicit Vertex(int x1 = 0, int y1 = 0) {
        x = x1;
        y = y1;
    }
};

typedef vector<Vertex> VertexList;

typedef bool (*IsInFunc)(Vertex &v, int edge);

typedef Vertex (*IntersectFunc)(Vertex &v1, Vertex &v2, int edge);

LPCSTR curs = IDC_CROSS;  //Initial Mouse Shape Cross.
int Round(int x) {
    return (int) (x + 0.5);
}


void swapp(int &x1, int &y1, int &x2, int &y2) {
    int temp = x1;
    x1 = x2;
    x2 = temp;
    temp = y1;
    y1 = y2;
    y2 = temp;
}


struct point {
    int x, y;

    point(int x, int y) : x(x), y(y) {}

    point() {}
};


static point _leftBottom;
static point _rightTop;
static point _rightBottom;
static point _leftTop;

void Draw8Points(double xc, double yc, double x, double y, COLORREF c, HDC hdc) {
    SetPixel(hdc, xc + x, yc + y, c);
    SetPixel(hdc, xc - x, yc + y, c);
    SetPixel(hdc, xc - x, yc - y, c);
    SetPixel(hdc, xc + x, yc - y, c);
    SetPixel(hdc, xc + y, yc + x, c);
    SetPixel(hdc, xc + y, yc - x, c);
    SetPixel(hdc, xc - y, yc + x, c);
    SetPixel(hdc, xc - y, yc - x, c);
}

void DrawCircleMidPoint(HDC hdc, double xc, double yc, int R, COLORREF c) {
    //f(x,y)=x*x+y*y-R*R
    //d=f(x+1,y-0.5) = (x+1)*(x+1)+ (y-0.5)*(y-0.5) - R*R
    double x = 0;
    double y = R;
    Draw8Points(xc, yc, x, y, c, hdc);
    double d = 0;
    while (x < y) {
        d = (x + 1) * (x + 1) + (y - 0.5) * (y - 0.5) - R * R;
        if (d < 0) {
            x++;
        } else {
            x++;
            y--;
        }
        Draw8Points(xc, yc, x, y, c, hdc);
    }
}

//////////////DDA Algorithm////////////////////////////////////////////
void DrawLine_DDA(HDC hdc, int x1, int y1, int x2, int y2, COLORREF c) {
    int dx = x2 - x1;
    int dy = y2 - y1;
    double slope = (double) dy / (double) dx;
    if (abs(dx) >= abs(dy)) {
        if (x1 > x2)
            swapp(x1, y1, x2, y2);
        int x = x1;
        double y = y1;
        SetPixel(hdc, x1, y1, c);
        while (x < x2) {
            x++;
            y += slope;
            SetPixel(hdc, x, y, c);
        }
    } else {
        if (y1 > y2)
            swapp(x1, y1, x2, y2);
        int y = y1;
        double x = x1;
        SetPixel(hdc, x1, y1, c);
        while (y < y2) {
            y++;
            x += 1 / slope;
            SetPixel(hdc, x, y, c);
        }
    }
}

//////////////MIDPOINT ALGORITHM FOR LINE //////////////////////////
void BresenhamLine(HDC hdc, int x1, int y1, int x2, int y2, int c) {
    int x, y, dx, dy, Adx, Ady, dAdy, dAdx, x_1, y_1, x_2, y_2;
    dx = x2 - x1;
    dy = y2 - y1;
    Adx = abs(dx);
    Ady = abs(dy);
    dAdy = 2 * Ady;
    dAdx = 2 * Adx;
    x_1 = dAdy - Adx;
    y_1 = dAdx - Ady;
    if (Ady <= Adx) {
        if (dx >= 0) {
            x = x1;
            y = y1;
            x_2 = x2;
        } else {
            x = x2;
            y = y2;
            x_2 = x1;
        }
        SetPixel(hdc, x, y, c);
        for (int i = 0; x < x_2; i++) {
            x++;
            if (x_1 < 0)
                x_1 += dAdy;
            else {
                if ((dx < 0 && dy < 0) || (dx > 0 && dy > 0))
                    y++;
                else
                    y--;
                x_1 += 2 * (Ady - Adx);
            }
            SetPixel(hdc, x, y, c);
        }
    } else {
        if (dy >= 0) {
            x = x1;
            y = y1;
            y_2 = y2;
        } else {
            x = x2;
            y = y2;
            y_2 = y1;
        }
        SetPixel(hdc, x, y, c);
        for (int i = 0; y < y_2; i++) {
            y++;
            if (y_1 <= 0)
                y_1 += dAdx;
            else {
                if ((dx < 0 && dy < 0) || (dx > 0 && dy > 0))
                    x++;
                else
                    x--;
                y_1 += 2 * (Adx - Ady);
            }
            SetPixel(hdc, x, y, c);
        }
    }
}

/////////////////////////////////////////////////////////////////////////
//////////////PARAMETRIC ALGORITHM///////////////////////////////////////
void DrawLine_Parametric(HDC hdc, int x1, int y1, int x2, int y2, COLORREF c) {

    double dt = (double) 1.0 / std::max(abs(x2 - x1), abs(y2 - y1));
    for (double t = 0; t <= 1; t += dt) {
        double x = x1 + t * (x2 - x1);
        double y = y1 + t * (y2 - y1);
        SetPixel(hdc, Round(x), Round(y), c);
    }
}

/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Draw8point(HDC hdc, int xc, int yc, int x, int y, COLORREF c) {
    SetPixel(hdc, xc + x, yc + y, c);
    SetPixel(hdc, xc + x, yc - y, c);
    SetPixel(hdc, xc - x, yc + y, c);
    SetPixel(hdc, xc - x, yc - y, c);
    SetPixel(hdc, xc + y, yc - x, c);
    SetPixel(hdc, xc - y, yc + x, c);
    SetPixel(hdc, xc - y, yc - x, c);
    SetPixel(hdc, xc + y, yc + x, c);


}

///////////////////////////////////////////////////////////////////////
////////////////////////Direct circle/////////////////////////////////
void DrawCircle_Direct(HDC hdc, int xc, int yc, int R, COLORREF c) {
    double x = 0;
    double y = R;
    Draw8point(hdc, xc, yc, R, 0, c);
    while (x < y) {
        x++;
        y = std::sqrt(R * R - x * x);
        Draw8point(hdc, xc, yc, x, y, c);
    }
}

////////////////////////////////////////////////////////////////////
///////////////////Polar circle/////////////////////////////////////
void DrawCircle_polar(HDC hdc, int xc, int yc, int R, COLORREF color) {
    double dtheta = 1.0 / R;
    double c = std::cos(dtheta);
    double s = std::sin(dtheta);
    double x = R;
    double y = 0;
    Draw8point(hdc, xc, yc, R, 0, color);
    while (x > y) {
        double x1 = x * c - y * s;
        y = x * s + y * c;
        x = x1;
        Draw8point(hdc, xc, yc, Round(x), Round(y), color);
    }
}

/////////////////////////////////////////////////////////////////
void Draw_Midpoint_circle(HDC hdc, int xc, int yc, int r, COLORREF color) {
    int x = 0;
    int y = r;
    int d = 1 - r;
    Draw8point(hdc, xc, yc, x, y, color);
    while (x < y) {
        if (d < 0)
            d += 2 * x + 3;
        else {
            d += 2 * (x - y) + 5;
            y--;
        }
        x++;
        Draw8point(hdc, xc, yc, x, y, color);
    }
}

/////////////////////////////////////////////////////////////////////
/////////////////Midpoint circle modified ///////////////////////////
void DrawCircle_MidPoint_Modified(HDC hdc, int xc, int yc, int R, COLORREF color) {
    int x = 0;
    int y = R;
    int d = 1 - R;
    int d1 = 3;
    int d2 = 5 - 2 * R;
    Draw8point(hdc, xc, yc, x, y, color);
    while (x < y) {
        if (d < 0) {
            d += d1;
            d1 += 2;
            d2 += 2;
        } else {
            d += d2;
            d2 += 4;
            d1 += 2;
            y--;
        }
        x++;
        Draw8point(hdc, xc, yc, x, y, color);
    }
}

////////////////////////////////////////////////////////////////////////////////
////////////////filing quarter/////////////////////////////////////////////////
void Draw8pointFilling(HDC hdc, int xc, int yc, int x, int y, int quarter, COLORREF color) {
    SetPixel(hdc, xc + x, yc + y, color);
    SetPixel(hdc, xc - x, yc + y, color);
    SetPixel(hdc, xc + x, yc - y, color);
    SetPixel(hdc, xc - x, yc - y, color);
    SetPixel(hdc, xc + y, yc + x, color);
    SetPixel(hdc, xc - y, yc + x, color);
    SetPixel(hdc, xc + y, yc - x, color);
    SetPixel(hdc, xc - y, yc - x, color);
    if (quarter == 1) {
        DrawLine_DDA(hdc, xc, yc, xc + x, yc - y, color);
        DrawLine_DDA(hdc, xc, yc, xc + y, yc - x, color);
    } else if (quarter == 2) {
        DrawLine_DDA(hdc, xc, yc, xc - x, yc - y, color);
        DrawLine_DDA(hdc, xc, yc, xc - y, yc - x, color);

    } else if (quarter == 3) {
        DrawLine_DDA(hdc, xc, yc, xc - x, yc + y, color);
        DrawLine_DDA(hdc, xc, yc, xc - y, yc + x, color);

    } else if (quarter == 4) {
        DrawLine_DDA(hdc, xc, yc, xc + x, yc + y, color);
        DrawLine_DDA(hdc, xc, yc, xc + y, yc + x, color);
    }
}

void DrawCircleFilling(HDC hdc, int xc, int yc, int R, int quarter, COLORREF color) {
    double dtheta = 1.0 / R;
    double c = std::cos(dtheta);
    double s = std::sin(dtheta);
    double x = R;
    double y = 0;
    Draw8pointFilling(hdc, xc, yc, R, 0, quarter, color);
    while (x > y) {
        double x1 = x * c - y * s;
        y = x * s + y * c;
        x = x1;
        Draw8pointFilling(hdc, xc, yc, Round(x), Round(y), quarter, color);
    }
}

///////////////////////////////////Filter with circle////////////////
void Draw8Point(HDC hdc, int xc, int yc, int x, int y, int q, COLORREF c) {

    if (q == 1) {
        //First quarter
        SetPixel(hdc, xc + x, yc - y, c);
        SetPixel(hdc, xc + y, yc - x, c);
    }

    if (q == 2) {

        //Second quarter
        SetPixel(hdc, xc - x, yc - y, c);
        SetPixel(hdc, xc - y, yc - x, c);

    }
    if (q == 3) {
        //Third quarter
        SetPixel(hdc, xc - x, yc + y, c);
        SetPixel(hdc, xc - y, yc + x, c);
    }
    if (q == 4) {
        //Fourth quarter
        SetPixel(hdc, xc + x, yc + y, c);
        SetPixel(hdc, xc + y, yc + x, c);
    }
}

void DrawCircleCartes(HDC hdc, double xc, double yc, int R, int q, COLORREF c) {
    //x*x +y*y =R*R
    //y=+=sqrt(R*R-x*x)
    double x = 0;
    double y = R;
    Draw8Point(hdc, xc, yc, Round(x), Round(y), q, c);
    while (x < y) {
        x++;
        y = sqrt(R * R - x * x);
        Draw8Point(hdc, xc, yc, Round(x), Round(y), q, c);
    }
}

void DrawSolve(HDC hdc, double xc, double yc, int R, int q, COLORREF c) {

    double x = 0;
    double y = R;
    SetPixel(hdc, xc + x, yc + y, c);
    SetPixel(hdc, xc - x, yc + y, c);
    SetPixel(hdc, xc - x, yc - y, c);
    SetPixel(hdc, xc + x, yc - y, c);
    SetPixel(hdc, xc + y, yc + x, c);
    SetPixel(hdc, xc + y, yc - x, c);
    SetPixel(hdc, xc - y, yc + x, c);
    SetPixel(hdc, xc - y, yc - x, c);
    while (x < y) {
        x++;
        y = sqrt(R * R - x * x);
        SetPixel(hdc, xc + x, yc + y, RGB(0, 0, 0));
        SetPixel(hdc, xc - x, yc + y, RGB(0, 0, 0));
        SetPixel(hdc, xc - x, yc - y, RGB(0, 0, 0));
        SetPixel(hdc, xc + x, yc - y, RGB(0, 0, 0));
        SetPixel(hdc, xc + y, yc + x, RGB(0, 0, 0));
        SetPixel(hdc, xc + y, yc - x, RGB(0, 0, 0));
        SetPixel(hdc, xc - y, yc + x, RGB(0, 0, 0));
        SetPixel(hdc, xc - y, yc - x, RGB(0, 0, 0));
    }
    // R--;
    while (R != 0) {
        DrawCircleCartes(hdc, xc, yc, R, q, c);
        R--;
    }
}

///////////Direct Ellipse //////////////////////////////////////////////////
void draw4point(HDC hdc, int x, int y, int xc, int yc, COLORREF c) {
    SetPixel(hdc, xc + x, yc + y, c);
    SetPixel(hdc, xc - x, yc + y, c);
    SetPixel(hdc, xc + x, yc - y, c);
    SetPixel(hdc, xc - x, yc - y, c);
}

void directellipse(HDC hdc, int xc, int yc, int R, int r2, COLORREF c) {
    double x = R;
    int y = 0;

    draw4point(hdc, R, 0, xc, yc, c);
    while (x >= 0) {
        x -= 0.01;
        y = abs(r2 * sqrt(1 - (pow(x, 2) / pow(R, 2))));
        draw4point(hdc, x, y, xc, yc, c);
    }
}

//////////////////////////////////////////////////////////////////////
////////////////Polar Ellipse ///////////////////////////////////////
void DrawEllipsePolar(HDC hdc, int xc, int yc, int r, int r2, COLORREF c) {
    double dtheta = 1.0 / r;
    for (double theta = 0.0; theta < 6.28; theta += dtheta) {
        double x = xc + r2 * cos(theta);
        double y = yc + r * sin(theta);
        SetPixel(hdc, Round(x), Round(y), c);
    }
}


////////////////////////////////////////////////////////////////////////////////////
///////////////////////////CLIPPING LINE////////////////////////////////////////////

union OutCode {
    unsigned All: 4;
    struct {
        unsigned left: 1, top: 1, right: 1, bottom: 1;
    };
};

OutCode GetOutCode(double x, double y, int xleft, int ytop, int xright, int ybottom) {
    OutCode out;
    out.All = 0;
    if (x < xleft)
        out.left = 1;
    else if (x > xright)
        out.right = 1;
    if (y < ytop)
        out.top = 1;
    else if (y > ybottom)
        out.bottom = 1;
    return out;
}

void VIntersect(double xs, double ys, double xe, double ye, int x, double *xi, double *yi) {
    *xi = x;
    *yi = ys + (x - xs) * (ye - ys) / (xe - xs);
}

void HIntersect(double xs, double ys, double xe, double ye, int y, double *xi, double *yi) {
    *yi = y;
    *xi = xs + (y - ys) * (xe - xs) / (ye - ys);
}

void CohenSuth(HDC hdc, int xs, int ys, int xe, int ye, int xleft, int ytop, int xright, int ybottom, COLORREF colol) {
    double x1 = xs, y1 = ys, x2 = xe, y2 = ye;
    OutCode out1 = GetOutCode(x1, y1, xleft, ytop, xright, ybottom);
    OutCode out2 = GetOutCode(x2, y2, xleft, ytop, xright, ybottom);
    while ((out1.All || out2.All) && !(out1.All & out2.All)) {
        double xi, yi;
        if (out1.All) {
            if (out1.left)
                VIntersect(x1, y1, x2, y2, xleft, &xi, &yi);
            else if (out1.top)
                HIntersect(x1, y1, x2, y2, ytop, &xi, &yi);
            else if (out1.right)
                VIntersect(x1, y1, x2, y2, xright, &xi, &yi);
            else
                HIntersect(x1, y1, x2, y2, ybottom, &xi, &yi);
            x1 = xi;
            y1 = yi;
            out1 = GetOutCode(x1, y1, xleft, ytop, xright, ybottom);
        } else {
            if (out2.left)
                VIntersect(x1, y1, x2, y2, xleft, &xi, &yi);
            else if (out2.top)
                HIntersect(x1, y1, x2, y2, ytop, &xi, &yi);
            else if (out2.right)
                VIntersect(x1, y1, x2, y2, xright, &xi, &yi);
            else
                HIntersect(x1, y1, x2, y2, ybottom, &xi, &yi);
            x2 = xi;
            y2 = yi;
            out2 = GetOutCode(x2, y2, xleft, ytop, xright, ybottom);
        }
    }
    if (!out1.All && !out2.All) {

        MoveToEx(hdc, round(x1), round(y1), NULL);
        DrawLine_DDA(hdc, x1, y1, x2, y2, colol);
    }
}

void CliipingLineWithCircle(HDC hdc, int xs, int ys, int xe, int ye, int xc, int yc, int r, COLORREF color) {
    double dt = (double) 1.0 / max(abs(xe - xs), abs(ye - ys));
    for (double t = 0; t <= 1; t += dt) {
        double x = xs + t * (xe - xs);
        double y = ys + t * (ye - ys);
        int test = (x - xc) * (x - xc) + (y - yc) * (y - yc) - (r) * (r);
        if (test <= 0) {
            SetPixel(hdc, x, y, color);
        }
    }
}

//////////////////////////////////////////////////////////////////////////
////////////////////////CLIPPING POINT/////////////////////////////////////

void PointClipping(HDC hdc, int x, int y, int xleft, int ytop, int xright, int ybottom, COLORREF color) {
    if (x >= xleft && x <= xright && y >= ytop && y <= ybottom)
        SetPixel(hdc, x, y, color);
}

void PointClipping(HDC hdc, int x, int y, int xc, int yc, int r, COLORREF color) {
    int test = (x - xc) * (x - xc) + (y - yc) * (y - yc) - (r) * (r);
    //cout<<r<<endl;
    if (test <= 0) {
        SetPixel(hdc, x, y, color);
        //cout<<"Value is :"<<test<<endl;
    }
    //cout<<"Value is :"<<test<<endl;
}

//////////////////////////////////////////////////////////////////////////
///////////////////////////CLIPPING POLYGON//////////////////////////////

VertexList ClipWithEdge(VertexList p, int edge, IsInFunc In, IntersectFunc Intersect) {
    VertexList OutList;
    Vertex v1 = p[p.size() - 1];
    bool v1_in = In(v1, edge);
    for (int i = 0; i < (int) p.size(); i++) {
        Vertex v2 = p[i];
        bool v2_in = In(v2, edge);
        if (!v1_in && v2_in) {
            OutList.push_back(Intersect(v1, v2, edge));
            OutList.push_back(v2);
        } else if (v1_in && v2_in) OutList.push_back(v2);
        else if (v1_in) OutList.push_back(Intersect(v1, v2, edge));
        v1 = v2;
        v1_in = v2_in;
    }
    return OutList;
}

bool InLeft(Vertex &v, int edge) {
    return v.x >= edge;
}

bool InRight(Vertex &v, int edge) {
    return v.x <= edge;
}

bool InTop(Vertex &v, int edge) {
    return v.y >= edge;
}

bool InBottom(Vertex &v, int edge) {
    return v.y <= edge;
}

Vertex VIntersect(Vertex &v1, Vertex &v2, int xedge) {
    Vertex res;
    res.x = xedge;
    res.y = v1.y + (xedge - v1.x) * (v2.y - v1.y) / (v2.x - v1.x);
    return res;
}

Vertex HIntersect(Vertex &v1, Vertex &v2, int yedge) {
    Vertex res;
    res.y = yedge;
    res.x = v1.x + (yedge - v1.y) * (v2.x - v1.x) / (v2.y - v1.y);
    return res;
}

void PolygonClip(HDC hdc, POINT *P, int Size, int xleft, int ytop, int xright, int ybottom) {
    VertexList vlist;
    for (int i = 0; i < Size; i++)vlist.push_back(Vertex(P[i].x, P[i].y));
    vlist = ClipWithEdge(vlist, xleft, InLeft, VIntersect);
    vlist = ClipWithEdge(vlist, ytop, InTop, HIntersect);
    vlist = ClipWithEdge(vlist, xright, InRight, VIntersect);
    vlist = ClipWithEdge(vlist, ybottom, InBottom, HIntersect);
    Vertex v1 = vlist[vlist.size() - 1];
    for (int i = 0; i < (int) vlist.size(); i++) {
        Vertex v2 = vlist[i];
        MoveToEx(hdc, Round(v1.x), Round(v1.y), NULL);
        LineTo(hdc, Round(v2.x), Round(v2.y));
        v1 = v2;
    }
}

struct Vector {
    double v[2];

    Vector(double x = 0, double y = 0) {
        v[0] = x;
        v[1] = y;
    }

    double &operator[](int i) {
        return v[i];
    }
};

struct Vector2 {
    double x, y;

    Vector2(double a = 0, double b = 0) {
        x = a;
        y = b;
    }
};

class Vector4 {
    double v[4];
public:
    Vector4(double a = 0, double b = 0, double c = 0, double d = 0) {
        v[0] = a;
        v[1] = b;
        v[2] = c;
        v[3] = d;
    }

    Vector4(double a[]) {
        memcpy(v, a, 4 * sizeof(double));
    }

    double &operator[](int i) {
        return v[i];
    }
};

class Matrix4 {
    Vector4 M[4];
public:
    Matrix4(double A[]) {
        memcpy(M, A, 16 * sizeof(double));
    }

    Vector4 &operator[](int i) {
        return M[i];
    }
};

Vector4 operator*(Matrix4 M, Vector4 &b) // right multiplication of M by b
{
    Vector4 res;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            res[i] += M[i][j] * b[j];
    return res;
}

double DotProduct(Vector4 &a, Vector4 &b) //multiplying a raw vector by a column vector
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

Vector4 GetHermiteCoeff(double x0, double s0, double x1, double s1) {
    static double H[16] = {2, 1, -2, 1, -3, -2, 3, -1, 0, 1, 0, 0, 1, 0, 0, 0};
    static Matrix4 basis(H);
    Vector4 v(x0, s0, x1, s1);
    return basis * v;
}

void DrawHermiteCurve2(HDC hdc, Vector2 &P0, Vector2 &T0, Vector2 &P1, Vector2 &T1, int
numpoints) {
    Vector4 xcoeff = GetHermiteCoeff(P0.x, T0.x, P1.x, T1.x);
    Vector4 ycoeff = GetHermiteCoeff(P0.y, T0.y, P1.y, T1.y);
    if (numpoints < 2)return;
    double dt = 1.0 / (numpoints - 1);
    for (double t = 0; t <= 1; t += dt) {
        Vector4 vt;
        vt[3] = 1;
        for (int i = 2; i >= 0; i--)vt[i] = vt[i + 1] * t;
        int x = round(DotProduct(xcoeff, vt));
        int y = round(DotProduct(ycoeff, vt));
        if (t == 0)MoveToEx(hdc, x, y, NULL);
        else LineTo(hdc, x, y);
    }
}

void DrawHermiteCurve(HDC hdc, int px1, int py1, int Tx1, int Ty1, int px2, int py2, int Tx2, int Ty2, COLORREF c) {
    double a0 = px1, a1 = Tx1,
            a2 = -3 * px1 - 2 * Tx1 + 3 * px2 - Tx2,
            a3 = 2 * px1 + Tx1 - 2 * px2 + Tx2;
    double b0 = py1, b1 = Ty1,
            b2 = -3 * py1 - 2 * Ty1 + 3 * py2 - Ty2,
            b3 = 2 * py1 + Ty1 - 2 * py2 + Ty2;
    for (double t = 0; t <= 1; t += 0.001) {
        double t2 = t * t, t3 = t2 * t;
        double x = a0 + a1 * t + a2 * t2 + a3 * t3;
        double y = b0 + b1 * t + b2 * t2 + b3 * t3;
        SetPixel(hdc, Round(x), Round(y), c);
    }
}

void DrawLine1(HDC hdc, int x1, int y1, int x2, int y2, COLORREF c) {
    int dx = x2 - x1;
    int dy = y2 - y1;
    if (abs(dy) <= abs(dx)) {
        if (x1 > x2)swapp(x1, y1, x2, y2);
        SetPixel(hdc, x1, y1, c);
        int x = x1;
        while (x < x2) {
            x++;
            double y = y1 + (double) (x - x1) * dy / dx;
            SetPixel(hdc, x, Round(y), c);

        }
    } else {
        if (y1 > y2)swapp(x1, y1, x2, y2);
        SetPixel(hdc, x1, y1, c);
        int y = y1;
        while (y < y2) {
            y++;
            double x = x1 + (double) (y - y1) * dx / dy;
            SetPixel(hdc, Round(x), y, c);
        }
    }

}

void DrawCardinalSpline(HDC hdc, Vector2 P[], int n, double c, int numpix) {
    double c1 = 1 - c;
    Vector2 T0(c1 * (P[2].x - P[0].x), c1 * (P[2].y - P[0].y));
    for (int i = 2; i < n - 1; i++) {
        Vector2 T1(c1 * (P[i + 1].x - P[i - 1].x), c1 * (P[i + 1].y - P[i - 1].y));
        DrawHermiteCurve2(hdc, P[i - 1], T0, P[i], T1, numpix);
        T0 = T1;
    }
}

void DrawSquare(HDC hdc, int x1, int y1, int x2, int y2, COLORREF c, COLORREF c2) {
    int distance = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
    int x3 = x2, y3 = y2 + distance, x4 = x1, y4 = y1 + distance;
    DrawLine1(hdc, x1, y1, x2, y2, c);
    DrawLine1(hdc, x2, y2, x3, y3, c);
    DrawLine1(hdc, x3, y3, x4, y4, c);
    DrawLine1(hdc, x4, y4, x1, y1, c);
    int z = x1 + 1;
    while (z < x2) {
        DrawHermiteCurve(hdc, z, y1 + 1, 0, 0, z, y4 - 1, 0, 0, c2);
        z++;
    }

}

void DrawBezierCurve(HDC hdc, int px1, int py1, int px2, int py2, int px3, int py3, int px4, int py4, COLORREF c) {
    int tx1 = 0, tx2 = 0;
    int ty1 = 0, ty2 = 0;
    DrawHermiteCurve(hdc, px1, py1, tx1, ty1, px3, py3, tx2, ty2, c);


}

void DrawRectangle(HDC hdc, int x1, int y1, int x2, int y2, int x3, int y3, COLORREF c, COLORREF c2) {
    int length = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)), width = sqrt(pow(x3 - x2, 2) + pow(y3 - y2, 2));
    int x4 = x1, y4 = y1 + width;
    DrawLine1(hdc, x1, y1, x2, y2, c);
    DrawLine1(hdc, x2, y2, x3, y3, c);
    DrawLine1(hdc, x3, y3, x4, y4, c);
    DrawLine1(hdc, x4, y4, x1, y1, c);
    int z = y1 + 1;
    while (z < y4) {
        DrawBezierCurve(hdc, x1 + 1, z, 0, 0, x2 - 1, z, 0, 0, c2);
        z++;
    }
}

struct Entry {
    int xmin, xmax;
};
#define MAXENTRIES 600

void InitEntries(Entry table[]) {
    for (int i = 0; i < MAXENTRIES; i++) {
        table[i].xmin = INT_MAX;
        table[i].xmax = -INT_MAX;
    }
}

void ScanEdge(POINT v1, POINT v2, Entry table[]) {
    if (v1.y == v2.y)return;
    if (v1.y > v2.y)swap(v1, v2);
    double minv = (double) (v2.x - v1.x) / (v2.y - v1.y);
    double x = v1.x;
    int y = v1.y;
    while (y < v2.y) {
        if (x < table[y].xmin)table[y].xmin = (int) ceil(x);
        if (x > table[y].xmax)table[y].xmax = (int) floor(x);
        y++;
        x += minv;
    }
}

void DrawSanLines(HDC hdc, Entry table[], COLORREF color) {
    for (int y = 0; y < MAXENTRIES; y++)
        if (table[y].xmin < table[y].xmax)
            for (int x = table[y].xmin; x <= table[y].xmax; x++)
                SetPixel(hdc, x, y, color);
}

void ConvexFill(HDC hdc, POINT p[], int n, COLORREF color) {
    Entry *table = new Entry[MAXENTRIES];
    InitEntries(table);
    POINT v1 = p[n - 1];
    for (int i = 0; i < n; i++) {
        POINT v2 = p[i];
        ScanEdge(v1, v2, table);
        v1 = p[i];
    }
    DrawSanLines(hdc, table, color);
    delete table;
}

struct EdgeRec {
    double x;
    double minv;
    int ymax;

    bool operator<(EdgeRec r) {
        return x < r.x;
    }
};

typedef list<EdgeRec> EdgeList;

EdgeRec InitEdgeRec(POINT &v1, POINT &v2) {
    if (v1.y > v2.y)swap(v1, v2);
    EdgeRec rec;
    rec.x = v1.x;
    rec.ymax = v2.y;
    rec.minv = (double) (v2.x - v1.x) / (v2.y - v1.y);
    return rec;
}

void InitEdgeTable(POINT *polygon, int n, EdgeList table[]) {
    POINT v1 = polygon[n - 1];
    for (int i = 0; i < n; i++) {
        POINT v2 = polygon[i];
        if (v1.y == v2.y) {
            v1 = v2;
            continue;
        }
        EdgeRec rec = InitEdgeRec(v1, v2);
        table[v1.y].push_back(rec);
        v1 = polygon[i];
    }
}

void GeneralPolygonFill(HDC hdc, POINT *polygon, int n, COLORREF c) {
    EdgeList *table = new EdgeList[MAXENTRIES];
    InitEdgeTable(polygon, n, table);
    int y = 0;
    while (y < MAXENTRIES && table[y].size() == 0)y++;
    if (y == MAXENTRIES)return;
    EdgeList ActiveList = table[y];
    while (ActiveList.size() > 0) {
        ActiveList.sort();
        for (EdgeList::iterator it = ActiveList.begin(); it != ActiveList.end(); it++) {
            int x1 = (int) ceil(it->x);
            it++;
            int x2 = (int) floor(it->x);
            for (int x = x1; x <= x2; x++)SetPixel(hdc, x, y, c);
        }
        y++;
        EdgeList::iterator it = ActiveList.begin();
        while (it != ActiveList.end())
            if (y == it->ymax) it = ActiveList.erase(it);
            else it++;
        for (EdgeList::iterator it = ActiveList.begin(); it != ActiveList.end(); it++)
            it->x += it->minv;
        ActiveList.insert(ActiveList.end(), table[y].begin(), table[y].end());
    }
    delete[] table;
}

//////////////////////////////FLOODFILL////////////////////

void FloodFill(HDC hdc, int x, int y, COLORREF Cb, COLORREF Cf) {
    COLORREF C = GetPixel(hdc, x, y);
    if (C == Cb || C == Cf)return;
    SetPixel(hdc, x, y, Cf);
    FloodFill(hdc, x + 1, y, Cb, Cf);
    FloodFill(hdc, x - 1, y, Cb, Cf);
    FloodFill(hdc, x, y + 1, Cb, Cf);
    FloodFill(hdc, x, y - 1, Cb, Cf);
}

void NRFloodFill(HDC hdc, int x, int y, COLORREF Cb, COLORREF Cf) {
    stack<Vertex> S;
    S.emplace(x, y);
    while (!S.empty()) {
        Vertex v = S.top();
        S.pop();
        COLORREF c = GetPixel(hdc, v.x, v.y);
        if (c == Cb || c == Cf)continue;
        SetPixel(hdc, v.x, v.y, Cf);
        S.emplace(v.x + 1, v.y);
        S.emplace(v.x - 1, v.y);
        S.emplace(v.x, v.y + 1);
        S.emplace(v.x, v.y - 1);
    }
}


POINT P[5];
Vector2 Vec[5];

/////////////////////////////////////////////
////////////Save Points/////////////////////
struct Save_Point {
    string Function_Name;
    int x1 = 0, x2 = 0, x3 = 0, x4 = 0, x5 = 0;
    int y1 = 0, y2 = 0, y3 = 0, y4 = 0, y5 = 0;
    int R = 0, r2 = 0;
    int Rc, Gc, Bc;
    int Rc2, Gc2, Bc2;
    int quarter;
    point _leftBottom = point(0, 0);
    point _rightTop = point(0, 0);
    point _rightBottom = point(0, 0);
    point _leftTop = point(0, 0);
    int xleft = 0;
    int ytop = 0;
    int xright = 0;
    int ybottom = 0;

    Save_Point(string Function_Name, int x1, int y1, int x2, int y2, int Rc, int Gc, int Bc)  //Line
    {
        this->Function_Name = Function_Name;
        this->x1 = x1;
        this->x2 = x2;
        this->y1 = y1;
        this->y2 = y2;
        this->Rc = Rc;
        this->Gc = Gc;
        this->Bc = Bc;
    }

    Save_Point(string Function_Name, int x1, int y1, int x2, int y2, int Rc, int Gc, int Bc, int Rc2, int Gc2,
               int Bc2)   //Square
    {
        this->Function_Name = Function_Name;
        this->x1 = x1;
        this->x2 = x2;
        this->y1 = y1;
        this->y2 = y2;
        this->Rc = Rc;
        this->Gc = Gc;
        this->Bc = Bc;

        this->Rc2 = Rc2;
        this->Gc2 = Gc2;
        this->Bc2 = Bc2;
    }

    Save_Point(string Function_Name, int x1, int y1, int x2, int y2, int x3, int y3, int Rc, int Gc, int Bc, int Rc2,
               int Gc2, int Bc2)     //Rectangle
    {
        this->Function_Name = Function_Name;
        this->x1 = x1;
        this->x2 = x2;
        this->x3 = x3;

        this->y1 = y1;
        this->y2 = y2;
        this->y3 = y3;

        this->Rc = Rc;
        this->Gc = Gc;
        this->Bc = Bc;

        this->Rc2 = Rc2;
        this->Gc2 = Gc2;
        this->Bc2 = Bc2;
    }

    Save_Point(string Function_Name, int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4, int x5, int y5,
               int Rc, int Gc, int Bc)        //Convex and Non convex filling
    {
        this->Function_Name = Function_Name;
        this->x1 = x1;
        this->x2 = x2;
        this->x3 = x3;
        this->x4 = x4;
        this->x5 = x5;

        this->y1 = y1;
        this->y2 = y2;
        this->y3 = y3;
        this->y4 = y4;
        this->y5 = y5;

        this->Rc = Rc;
        this->Gc = Gc;
        this->Bc = Bc;


    }

    Save_Point(string Function_Name, int x1, int y1, int x2, int y2, int xleft, int ytop, int xright, int ybottom,
               int Rc, int Gc, int Bc)//clipping line
    {
        this->Function_Name = Function_Name;
        this->x1 = x1;
        this->x2 = x2;
        this->y1 = y1;
        this->y2 = y2;
        this->xleft = xleft;
        this->xright = xright;
        this->ytop = ytop;
        this->ybottom = ybottom;
        this->Rc = Rc;
        this->Gc = Gc;
        this->Bc = Bc;
    }

    Save_Point(string Function_Name, int x1, int y1, int R, int r2, int Rc, int Gc, int Bc, char e) //Ellipse
    {
        this->Function_Name = Function_Name;
        this->x1 = x1;
        this->R = R;
        this->r2 = r2;
        this->y1 = y1;
        this->Rc = Rc;
        this->Gc = Gc;
        this->Bc = Bc;
    }

    Save_Point(string Function_Name, int x1, int y1, int R, int Rc, int Gc, int Bc) //Circle
    {
        this->Function_Name = Function_Name;
        this->x1 = x1;
        this->R = R;
        this->y1 = y1;
        this->Rc = Rc;
        this->Gc = Gc;
        this->Bc = Bc;
    }

    Save_Point(string Function_Name, int x1, int y1, int R, int quarter, int Rc, int Gc, int Bc, string f) //Filling
    {
        this->Function_Name = Function_Name;
        this->x1 = x1;
        this->R = R;
        this->y1 = y1;
        this->Rc = Rc;
        this->Gc = Gc;
        this->Bc = Bc;
        this->quarter = quarter;
    }

    Save_Point(string Function_Name, int x1, int y1, int Rc, int Gc, int Bc, int Rc2, int Gc2,
               int Bc2) // For filling recursive and non recursive
    {
        this->Function_Name = Function_Name;
        this->x1 = x1;
        this->y1 = y1;
        this->Rc = Rc;
        this->Gc = Gc;
        this->Bc = Bc;
        this->Rc2 = Rc2;
        this->Gc2 = Gc2;
        this->Bc2 = Bc2;
    }

    Save_Point(string Function_Name, int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4, int x5, int y5,
               string shc)  // cardinal spline curve
    {
        this->Function_Name = Function_Name;
        this->x1 = x1;
        this->x2 = x2;
        this->x3 = x3;
        this->x4 = x4;
        this->x5 = x5;
        this->y1 = y1;
        this->y2 = y2;
        this->y3 = y3;
        this->y4 = y4;
        this->y5 = y5;
    }

};


vector<Save_Point> Arr_Save_Point;

void Save() {
    string data = "";

    for (auto it = Arr_Save_Point.begin(); it != Arr_Save_Point.end(); ++it) {
        /* Save Function of All Line */
        if (it->Function_Name == "DDLine") {
            data += it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->x2) +
                    "," + to_string(it->y2) + "," + to_string(it->Rc) + "," + to_string(it->Gc) + "," +
                    to_string(it->Bc) + "\n";
        } else if (it->Function_Name == "DPLine") {
            data += it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->x2) +
                    "," + to_string(it->y2) + "," + to_string(it->Rc) + "," + to_string(it->Gc) + "," +
                    to_string(it->Bc) + "\n";
        } else if (it->Function_Name == "DMLine") {
            data += it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->x2) +
                    "," + to_string(it->y2) + "," + to_string(it->Rc) + "," + to_string(it->Gc) + "," +
                    to_string(it->Bc) + "\n";
        }


        /* save Function of ALL Circle */

        if (it->Function_Name == "DDCircle") {
            data += it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->R) +
                    "," + to_string(it->Rc) + "," + to_string(it->Gc) + "," + to_string(it->Bc) + "\n";
        } else if (it->Function_Name == "DPCircle") {
            data += it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->R) +
                    "," + to_string(it->Rc) + "," + to_string(it->Gc) + "," + to_string(it->Bc) + "\n";
        } else if (it->Function_Name == "DMCircle") {
            data += it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->R) +
                    "," + to_string(it->Rc) + "," + to_string(it->Gc) + "," + to_string(it->Bc) + "\n";
        }


        /*save Function of All Filling */

        if (it->Function_Name == "DPCircleFilling") {
            data += it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->R) +
                    "," + to_string(it->quarter) + "," + to_string(it->Rc) + "," + to_string(it->Gc) + "," +
                    to_string(it->Bc) + "\n";
        } else if (it->Function_Name == "DSqHer") {
            data += it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->x2) +
                    "," + to_string(it->y2) + "," + to_string(it->Rc) + "," + to_string(it->Gc) + "," +
                    to_string(it->Bc) +
                    "," + to_string(it->Rc2) + "," + to_string(it->Gc2) + "," + to_string(it->Bc2) + "\n";
        } else if (it->Function_Name == "DRecBez") {
            data += it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->x2) +
                    "," + to_string(it->y2) + "," + to_string(it->x3) + "," + to_string(it->y3) + ","
                    + to_string(it->Rc) + "," + to_string(it->Gc) + "," + to_string(it->Bc) + "," + to_string(it->Rc2) +
                    "," + to_string(it->Gc2) + "," + to_string(it->Bc2) + "\n";
        } else if (it->Function_Name == "ConvexFill") {
            data += (it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->x2) +
                     "," + to_string(it->y2) + "," + to_string(it->x3) + "," + to_string(it->y3) + "," +
                     to_string(it->x4) + "," + to_string(it->y4) + "," + to_string(it->x5) + "," + to_string(it->y5)
                     + "," + to_string(it->Rc) + "," + to_string(it->Gc) + "," + to_string(it->Bc) + "\n");
        } else if (it->Function_Name == "GenFill") {
            data += (it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->x2) +
                     "," + to_string(it->y2) + "," + to_string(it->x3) + "," + to_string(it->y3) + "," +
                     to_string(it->x4) + "," + to_string(it->y4) + "," + to_string(it->x5) + "," + to_string(it->y5)
                     + "," + to_string(it->Rc) + "," + to_string(it->Gc) + "," + to_string(it->Bc) + "\n");
        }


        /* save Function of ALL Ellipse */
        if (it->Function_Name == "DDEllipse") {
            data += it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->R) +
                    "," + to_string(it->r2) + "," + to_string(it->Rc) + "," + to_string(it->Gc) + "," +
                    to_string(it->Bc) + "\n";
        } else if (it->Function_Name == "DPEllipse") {
            data += it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->R) +
                    "," + to_string(it->r2) + "," + to_string(it->Rc) + "," + to_string(it->Gc) + "," +
                    to_string(it->Bc) + "\n";
        }

        /* save line clipped on rectangle using line*/
        if (it->Function_Name == "CohenSuth") {
            data += it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->x2) +
                    "," + to_string(it->y2) + "," + to_string(it->xleft) + "," + to_string(it->ytop) + "," +
                    to_string(it->xright) + "," + to_string(it->ybottom) + "," + to_string(it->Rc) + "," +
                    to_string(it->Gc) + "," + to_string(it->Bc) + "\n";
        }
        /*save Function of  filling recursive and non-recursive  */
        if (it->Function_Name == "FillRecursive") {
            data += it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->Rc) +
                    "," + to_string(it->Gc) + "," + to_string(it->Bc) + "\n";
        } else if (it->Function_Name == "NonFillRecursive") {
            data += it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->Rc) +
                    "," + to_string(it->Gc) + "," + to_string(it->Bc) + "," + to_string(it->Rc2) + "," +
                    to_string(it->Gc2) + "," + to_string(it->Bc2) + "\n";
        } else if (it->Function_Name == "CardSpline") {
            data += (it->Function_Name + "," + to_string(it->x1) + "," + to_string(it->y1) + "," + to_string(it->x2) +
                     "," + to_string(it->y2) + "," + to_string(it->x3) + "," + to_string(it->y3) + "," +
                     to_string(it->x4) + "," + to_string(it->y4) + "," + to_string(it->x5) + "," + to_string(it->y5)
                     + "," + "\n");
        }
    }
    ofstream myfile;
    myfile.open("Graphics.txt");
    myfile << data;
    myfile.close();

}


//Load Points//
void Load(HDC hdc) {
    string Line;
    ifstream LoadFile;
    LoadFile.open("Graphics.txt");
    if (!LoadFile) {
        cout << "Unable to open file";
        return;
    }
    while (getline(LoadFile, Line)) {
        vector<string> Fun_Load;
        string buff;
        for (auto n: Line) {
            if (n != ',')
                buff += n;
            else if (n == ',' && buff != "") {
                Fun_Load.push_back(buff);
                buff = "";
            }
        }
        if (buff != "") Fun_Load.push_back(buff);

        /* Load Function of All Line */
        if (Fun_Load[0] == "DDLine") {
            //stoi() is a function to covert string to integer
            COLORREF c = RGB(stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7])); //the index of colors
            DrawLine_DDA(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]), c);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]));
            Arr_Save_Point.push_back(x);


        } else if (Fun_Load[0] == "DPLine") {
            COLORREF c = RGB(stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]));
            DrawLine_Parametric(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]), c);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]));
            Arr_Save_Point.push_back(x);
        } else if (Fun_Load[0] == "DMLine") {
            COLORREF c = RGB(stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]));
            BresenhamLine(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]), c);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]));
            Arr_Save_Point.push_back(x);
        }


        /* Load Function of ALL Circle */

        if (Fun_Load[0] == "DDCircle") {
            COLORREF c = RGB(stoi(Fun_Load[4]), stoi(Fun_Load[5]), stoi(Fun_Load[6]));
            DrawCircle_Direct(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), c);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]));
            Arr_Save_Point.push_back(x);
        } else if (Fun_Load[0] == "DPCircle") {
            COLORREF c = RGB(stoi(Fun_Load[4]), stoi(Fun_Load[5]), stoi(Fun_Load[6]));
            DrawCircle_polar(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), c);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]));
            Arr_Save_Point.push_back(x);
        } else if (Fun_Load[0] == "DMCircle") {
            COLORREF c = RGB(stoi(Fun_Load[4]), stoi(Fun_Load[5]), stoi(Fun_Load[6]));
            Draw_Midpoint_circle(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), c);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]));
            Arr_Save_Point.push_back(x);
        } else if (Fun_Load[0] == "DMMCircle") {
            COLORREF c = RGB(stoi(Fun_Load[4]), stoi(Fun_Load[5]), stoi(Fun_Load[6]));
            DrawCircle_MidPoint_Modified(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), c);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]));
            Arr_Save_Point.push_back(x);

        }


        /*Load Function of All Filling */
        if (Fun_Load[0] == "DPCircleFilling") {
            COLORREF c = RGB(stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]));
            DrawCircleFilling(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]), c);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]), "e");
            Arr_Save_Point.push_back(x);
        } else if (Fun_Load[0] == "DSqHer") {
            //stoi() is a function to covert string to integer
            COLORREF c = RGB(stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]));
            COLORREF c2 = RGB(stoi(Fun_Load[8]), stoi(Fun_Load[9]), stoi(Fun_Load[10])); //the index of colors
            DrawSquare(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]), c, c2);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]), stoi(Fun_Load[8]), stoi(Fun_Load[9]),
                         stoi(Fun_Load[10]));
            Arr_Save_Point.push_back(x);


        } else if (Fun_Load[0] == "DRecBez") {
            //stoi() is a function to covert string to integer
            COLORREF c = RGB(stoi(Fun_Load[7]), stoi(Fun_Load[8]), stoi(Fun_Load[9]));
            COLORREF c2 = RGB(stoi(Fun_Load[10]), stoi(Fun_Load[11]), stoi(Fun_Load[12])); //the index of colors
            DrawRectangle(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                          stoi(Fun_Load[5]), stoi(Fun_Load[6]), c, c2);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]), stoi(Fun_Load[8]), stoi(Fun_Load[9]),
                         stoi(Fun_Load[10]), stoi(Fun_Load[11]), stoi(Fun_Load[12]));
            Arr_Save_Point.push_back(x);


        } else if (Fun_Load[0] == "ConvexFill") {
            //stoi() is a function to covert string to integer
            COLORREF c = RGB(stoi(Fun_Load[11]), stoi(Fun_Load[12]), stoi(Fun_Load[13]));
            P[0].x = stoi(Fun_Load[1]);
            P[0].y = stoi(Fun_Load[2]);
            P[1].x = stoi(Fun_Load[3]);
            P[1].y = stoi(Fun_Load[4]);
            P[2].x = stoi(Fun_Load[5]);
            P[2].y = stoi(Fun_Load[6]);
            P[3].x = stoi(Fun_Load[7]);
            P[3].y = stoi(Fun_Load[8]);
            P[4].x = stoi(Fun_Load[9]);
            P[4].y = stoi(Fun_Load[10]);
            Polygon(hdc, P, 5);
            ConvexFill(hdc, P, 5, c);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]), stoi(Fun_Load[8]), stoi(Fun_Load[9]),
                         stoi(Fun_Load[10]), stoi(Fun_Load[11]), stoi(Fun_Load[12]), stoi(Fun_Load[13]));
            Arr_Save_Point.push_back(x);

        } else if (Fun_Load[0] == "GenFill") {
            //stoi() is a function to covert string to integer
            COLORREF c = RGB(stoi(Fun_Load[11]), stoi(Fun_Load[12]), stoi(Fun_Load[13]));
            P[0].x = stoi(Fun_Load[1]);
            P[0].y = stoi(Fun_Load[2]);
            P[1].x = stoi(Fun_Load[3]);
            P[1].y = stoi(Fun_Load[4]);
            P[2].x = stoi(Fun_Load[5]);
            P[2].y = stoi(Fun_Load[6]);
            P[3].x = stoi(Fun_Load[7]);
            P[3].y = stoi(Fun_Load[8]);
            P[4].x = stoi(Fun_Load[9]);
            P[4].y = stoi(Fun_Load[10]);
            Polygon(hdc, P, 5);
            GeneralPolygonFill(hdc, P, 5, c);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]), stoi(Fun_Load[8]), stoi(Fun_Load[9]),
                         stoi(Fun_Load[10]), stoi(Fun_Load[11]), stoi(Fun_Load[12]), stoi(Fun_Load[13]));
            Arr_Save_Point.push_back(x);

        }
        /* Load Function of  filling recursive and non-recursive */
        if (Fun_Load[0] == "FillRecursive") {
            COLORREF c = RGB(stoi(Fun_Load[3]), stoi(Fun_Load[4]), stoi(Fun_Load[5]));
            COLORREF c2 = RGB(stoi(Fun_Load[6]), stoi(Fun_Load[7]), stoi(Fun_Load[8]));
            FloodFill(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), c, c2);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]), stoi(Fun_Load[8]));
            Arr_Save_Point.push_back(x);

        }
        if (Fun_Load[0] == "NonFillRecursive") {
            COLORREF c = RGB(stoi(Fun_Load[3]), stoi(Fun_Load[4]), stoi(Fun_Load[5]));
            COLORREF c2 = RGB(stoi(Fun_Load[6]), stoi(Fun_Load[7]), stoi(Fun_Load[8]));
            NRFloodFill(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), c, c2);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]), stoi(Fun_Load[8]));
            Arr_Save_Point.push_back(x);

        } else if (Fun_Load[0] == "CardSpline") {
            //stoi() is a function to covert string to integer
            //COLORREF c =RGB(stoi(Fun_Load[11]),stoi(Fun_Load[12]),stoi(Fun_Load[13]));
            Vector2 Vec[5];
            Vec[0].x = stoi(Fun_Load[1]);
            Vec[0].y = stoi(Fun_Load[2]);
            Vec[1].x = stoi(Fun_Load[3]);
            Vec[1].y = stoi(Fun_Load[4]);
            Vec[2].x = stoi(Fun_Load[5]);
            Vec[2].y = stoi(Fun_Load[6]);
            Vec[3].x = stoi(Fun_Load[7]);
            Vec[3].y = stoi(Fun_Load[8]);
            Vec[4].x = stoi(Fun_Load[9]);
            Vec[4].y = stoi(Fun_Load[10]);

            DrawCardinalSpline(hdc, Vec, 5, 0.5, 50);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]), stoi(Fun_Load[8]), stoi(Fun_Load[9]),
                         stoi(Fun_Load[10]), "csc");
            Arr_Save_Point.push_back(x);

        }



        /* Load Function of ALL Ellipse */

        if (Fun_Load[0] == "DDEllipse") {
            COLORREF c = RGB(stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]));
            directellipse(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]), c);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]), 'e');
            Arr_Save_Point.push_back(x);
        } else if (Fun_Load[0] == "DPEllipse") {
            COLORREF c = RGB(stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]));
            DrawEllipsePolar(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]), c);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]), 'e');
            Arr_Save_Point.push_back(x);
        }

        /* Load Function of Draw Rectangle */
        if (Fun_Load[0] == "CohenSuth") {
            COLORREF c = RGB(stoi(Fun_Load[9]), stoi(Fun_Load[10]), stoi(Fun_Load[11]));
            CohenSuth(hdc, stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                      stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]), stoi(Fun_Load[8]), c);
            Save_Point x(Fun_Load[0], stoi(Fun_Load[1]), stoi(Fun_Load[2]), stoi(Fun_Load[3]), stoi(Fun_Load[4]),
                         stoi(Fun_Load[5]), stoi(Fun_Load[6]), stoi(Fun_Load[7]), stoi(Fun_Load[8]), stoi(Fun_Load[9]),
                         stoi(Fun_Load[10]), stoi(Fun_Load[11]));
            Arr_Save_Point.push_back(x);
        }

    }
}


//Clean Points//
void Clear() {
    Arr_Save_Point.clear();
    counter = 0;
    Num_of_Points = 0;
    //InvalidateRect(hWnd, NULL, TRUE);
}


/*  Declare Windows procedure  */
LRESULT CALLBACK WindowProcedure(HWND, UINT, WPARAM, LPARAM);

/*  Make the class name into a global variable  */
TCHAR szClassName[] = _T("Graphics_Project");

int WINAPI WinMain(HINSTANCE hThisInstance,
                   HINSTANCE hPrevInstance,
                   LPSTR lpszArgument,
                   int nCmdShow) {
    HWND hwnd;               /* This is the handle for our window */
    MSG messages;            /* Here messages to the application are saved */
    WNDCLASSEX wincl;        /* Data structure for the windowclass */

/* The Window structure */
    wincl.hInstance = hThisInstance;
    wincl.lpszClassName = szClassName;
    wincl.lpfnWndProc = WindowProcedure;      /* This function is called by windows */
    wincl.style = CS_DBLCLKS;                 /* Catch double-clicks */
    wincl.cbSize = sizeof(WNDCLASSEX);

/* Use default icon and mouse-pointer */
    wincl.hIcon = LoadIcon(NULL, IDI_APPLICATION);
    wincl.hIconSm = LoadIcon(NULL, IDI_APPLICATION);
    wincl.hCursor = LoadCursor(NULL, curs);
    wincl.lpszMenuName = NULL;                 /* No menu */
    wincl.cbClsExtra = 0;                      /* No extra bytes after the window class */
    wincl.cbWndExtra = 0;                      /* structure or the window instance */
/* Use Windows's default colour as the background of the window */
    wincl.hbrBackground = CreateSolidBrush(RGB(255, 255, 255));  /// background color white

/* Register the window class, and if it fails quit the program */
    if (!RegisterClassEx(&wincl))
        return 0;

/* The class is registered, let's create the program*/
    hwnd = CreateWindowEx(
            0,                   /* Extended possibilites for variation */
            szClassName,         /* Classname */
            _T("Graphics project App"),       /* Title Text */
            WS_OVERLAPPEDWINDOW, /* default window */
            CW_USEDEFAULT,       /* Windows decides the position */
            CW_USEDEFAULT,       /* where the window ends up on the screen */
            544,                 /* The programs width */
            375,                 /* and height in pixels */
            HWND_DESKTOP,        /* The window is a child-window to desktop */
            NULL,                /* No menu */
            hThisInstance,       /* Program Instance handler */
            NULL                 /* No Window Creation data */
    );

/* Make the window visible on the screen */
    ShowWindow(hwnd, nCmdShow);

/* Run the message loop. It will run until GetMessage() returns 0 */
    while (GetMessage(&messages, NULL, 0, 0)) {
/* Translate virtual-key messages into character messages */
        TranslateMessage(&messages);
/* Send message to WindowProcedure */
        DispatchMessage(&messages);
    }

/* The program return-value is 0 - The value that PostQuitMessage() gave */
    return messages.wParam;
}


/*  This function is called by the Windows function DispatchMessage()  */
HMENU hMenu;

void addMenu(HWND hwnd) {
    hMenu = CreateMenu();

    HMENU hFileMenu = CreateMenu();
    AppendMenu(hFileMenu, MF_STRING, 0, "Save");
    AppendMenu(hFileMenu, MF_STRING, 1, "Load");
    AppendMenu(hFileMenu, MF_STRING, 25, "Clean");
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) hFileMenu, "File");

    HMENU hLineMenu = CreateMenu();
    AppendMenu(hLineMenu, MF_STRING, 2, "DDA");
    AppendMenu(hLineMenu, MF_STRING, 3, "Mid Point");
    AppendMenu(hLineMenu, MF_STRING, 4, "Parametric");
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) hLineMenu, "Line");

    HMENU hEllipseMenu = CreateMenu();
    AppendMenu(hEllipseMenu, MF_STRING, 5, "Direct");
    AppendMenu(hEllipseMenu, MF_STRING, 6, "Polar");
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) hEllipseMenu, "Ellipse");

    HMENU hCircleMenu = CreateMenu();
    AppendMenu(hCircleMenu, MF_STRING, 7, "Direct,");
    AppendMenu(hCircleMenu, MF_STRING, 8, "Polar iterative Polar");
    AppendMenu(hCircleMenu, MF_STRING, 9, "Midpoint");
    AppendMenu(hCircleMenu, MF_STRING, 26, "Modified MidPoint");
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) hCircleMenu, "Circle");


    HMENU hColorMenu = CreateMenu();
    AppendMenu(hColorMenu, MF_STRING, 10, "Black");
    AppendMenu(hColorMenu, MF_STRING, 11, "Red");
    AppendMenu(hColorMenu, MF_STRING, 12, "Blue");
    AppendMenu(hColorMenu, MF_STRING, 13, "Green");
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) hColorMenu, "Colors");

    HMENU hFilling = CreateMenu();
    AppendMenu(hFilling, MF_STRING, 14, "Black");
    AppendMenu(hFilling, MF_STRING, 15, "Red");
    AppendMenu(hFilling, MF_STRING, 16, "Blue");
    AppendMenu(hFilling, MF_STRING, 17, "Green");
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) hFilling, "Fill Circle");


    HMENU hquarter = CreateMenu();
    AppendMenu(hquarter, MF_STRING, 18, "First");
    AppendMenu(hquarter, MF_STRING, 19, "Second");
    AppendMenu(hquarter, MF_STRING, 20, "Third");
    AppendMenu(hquarter, MF_STRING, 21, "Fourth");
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) hquarter, "Quarter Circle");

    HMENU hwindow = CreateMenu();
    AppendMenu(hwindow, MF_STRING, 50, "Clipping by Line(Rectangle)");
    AppendMenu(hwindow, MF_STRING, 51, "Clipping by Line(Square)");
    AppendMenu(hwindow, MF_STRING, 52, "Clipping by Line(Circle)");
    AppendMenu(hwindow, MF_STRING, 53, "Clipping by point(Rectangle)");
    AppendMenu(hwindow, MF_STRING, 54, "Clipping by Point(Square)");
    AppendMenu(hwindow, MF_STRING, 55, "Clipping by Point(Circle)");
    AppendMenu(hwindow, MF_STRING, 56, "Clipping by Polygon(Rectangle)");
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) hwindow, "Clipping");


    HMENU hSq = CreateMenu();
    AppendMenu(hSq, MF_STRING, 57, "Black");
    AppendMenu(hSq, MF_STRING, 58, "Red");
    AppendMenu(hSq, MF_STRING, 59, "Blue");
    AppendMenu(hSq, MF_STRING, 60, "Green");
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) hSq, "Fill Square Hermite");

    HMENU hRec = CreateMenu();
    AppendMenu(hRec, MF_STRING, 61, "Black");
    AppendMenu(hRec, MF_STRING, 62, "Red");
    AppendMenu(hRec, MF_STRING, 63, "Blue");
    AppendMenu(hRec, MF_STRING, 64, "Green");
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) hRec, "Fill Rectangle Bezier");

    HMENU hConvex = CreateMenu();
    AppendMenu(hConvex, MF_STRING, 65, "Black");
    AppendMenu(hConvex, MF_STRING, 67, "Green");
    AppendMenu(hConvex, MF_STRING, 68, "Blue");
    AppendMenu(hConvex, MF_STRING, 66, "RED");
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) hConvex, "Convex Filling");

    HMENU hGenFill = CreateMenu();
    AppendMenu(hGenFill, MF_STRING, 70, "Black");
    AppendMenu(hGenFill, MF_STRING, 71, "Red");
    AppendMenu(hGenFill, MF_STRING, 72, "Green");
    AppendMenu(hGenFill, MF_STRING, 73, "Blue");
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) hGenFill, "Non Convex Filling");


    HMENU hFillingByCircles = CreateMenu();
    AppendMenu(hFillingByCircles, MF_STRING, 27, "Black");
    AppendMenu(hFillingByCircles, MF_STRING, 28, "Red");
    AppendMenu(hFillingByCircles, MF_STRING, 29, "Blue");
    AppendMenu(hFillingByCircles, MF_STRING, 30, "Green");
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) hFillingByCircles, "Fill Circle by other circles");


    HMENU cursorMenu = CreateMenu();
    AppendMenu(cursorMenu, MF_STRING, 31, _T("Hand"));
    AppendMenu(cursorMenu, MF_STRING, 32, _T("Standard arrow"));
    AppendMenu(cursorMenu, MF_STRING, 33, _T("Crosshair"));
    AppendMenu(cursorMenu, MF_STRING, 34, _T("Arrow and question mark"));
    AppendMenu(cursorMenu, MF_STRING, 35, _T("Vertical arrow"));
    AppendMenu(cursorMenu, MF_STRING, 36, _T("I-beam"));
    AppendMenu(cursorMenu, MF_STRING, 37, _T("Hourglass"));
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) cursorMenu, _T("Cursor"));

    HMENU FloodMenu = CreateMenu();
    AppendMenu(FloodMenu, MF_STRING, 74, _T("Recursive Flood Fill"));
    AppendMenu(FloodMenu, MF_STRING, 75, _T("non-Recursive Flood Fill"));
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) FloodMenu, _T("Flood Fill"));
    HMENU CardinalSplines = CreateMenu();
    AppendMenu(CardinalSplines, MF_STRING, 76, _T("Cardinal Splines"));
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR) CardinalSplines, _T("Cardinal Splines"));
    SetMenu(hwnd, hMenu);
}

int Rc, Gc, Bc;
COLORREF c = 0;
COLORREF cF = 0;
int x1, x2, x3, x4, x5;
int xz;
int y_1, y2, y3, y4, y5, y = 0;
int r, r2;
int R;
int quarter;
int m;
int X_start, X_end, Y_start, Y_end, X_left, Y_top, X_right, Y_bottom;
int counter_ell = 0;
int Num_of_Points2 = 0, ind = 0;

LRESULT CALLBACK WindowProcedure(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam) {
    HDC hdc = GetDC(hwnd);
    switch (message)                  /* handle the messages */
    {
        case WM_LBUTTONDOWN:
            x1 = LOWORD(lParam);
            y_1 = HIWORD(lParam);

            if (ind == 0 && (m == 57 || m == 58 || m == 59 || m == 60 || m == 61 || m == 62 || m == 63 || m == 64)) {
                xz = LOWORD(lParam);
                y = HIWORD(lParam);
                ind = 1;
            }


            if (m == 50) //Line clipping (Rectangle)
            {
                if (Num_of_Points == 0) {
                    cout << "Num of Points : " << Num_of_Points << endl;
                    X_left = LOWORD(lParam);
                    Y_top = HIWORD(lParam);
                    Num_of_Points++;
                } else if (Num_of_Points == 1) {
                    cout << "Num of Points : " << Num_of_Points << endl;
                    X_right = LOWORD(lParam);
                    Y_bottom = HIWORD(lParam);
                    Rectangle(hdc, X_left, Y_top, X_right, Y_bottom);
                    Num_of_Points++;
                } else if (Num_of_Points == 2) {
                    cout << " Start Point Of line " << endl;
                    X_start = LOWORD(lParam);
                    Y_start = HIWORD(lParam);
                    Num_of_Points++;
                } else if (Num_of_Points == 3) {
                    cout << " End Point Of line " << endl;
                    X_end = LOWORD(lParam);
                    Y_end = HIWORD(lParam);
                    CohenSuth(hdc, X_start, Y_start, X_end, Y_end, X_left, Y_top, X_right, Y_bottom, c);
                    Num_of_Points = 2;
                }
                Save_Point x("CohenSuth", X_start, Y_start, X_end, Y_end, X_left, Y_top, X_right, Y_bottom, 0, 0, 0);
                Arr_Save_Point.push_back(x);
            } else if (m == 51) //Line clipping (Square)
            {
                if (Num_of_Points == 0) {
                    cout << "Counter " << Num_of_Points << endl;
                    X_left = LOWORD(lParam);
                    Y_top = HIWORD(lParam);
                    X_right = X_left + 200;
                    Y_bottom = Y_top + 200;
                    Rectangle(hdc, X_left, Y_top, X_right, Y_bottom);
                    cout << "Square Done.\n";
                    Num_of_Points++;
                } else if (Num_of_Points == 1) {
                    cout << "Start EndPoint.\n";
                    X_start = LOWORD(lParam);
                    Y_start = HIWORD(lParam);
                    Num_of_Points++;
                } else if (Num_of_Points == 2) {
                    cout << "Line Clipping Done.\n";
                    X_end = LOWORD(lParam);
                    Y_end = HIWORD(lParam);
                    CohenSuth(hdc, X_start, Y_start, X_end, Y_end, X_left, Y_top, X_right, Y_bottom, c);
                    Num_of_Points = 1;
                }
                Save_Point x("CohenSuth", X_start, Y_start, X_end, Y_end, X_left, Y_top, X_right, Y_bottom, 0, 0, 0);
                Arr_Save_Point.push_back(x);
            } else if (m == 52)  //Line Clipping (Circle)
            {
                if (counter == 0) {
                    cout << "Center Of Circle. " << endl;
                    X_left = LOWORD(lParam);      //Center
                    Y_top = HIWORD(lParam);
                    counter++;
                } else if (counter == 1) {
                    cout << "Circle Window Done " << endl;
                    X_right = LOWORD(lParam);
                    Y_bottom = HIWORD(lParam);
                    r = sqrt(pow((X_right - X_left), 2) + pow((Y_bottom - Y_top), 2));
                    DrawCircleMidPoint(hdc, X_left, Y_top, r, RGB(0, 0, 0));     //radius.
                    counter++;
                } else if (counter == 2) {
                    cout << "Start Point" << endl;
                    X_start = LOWORD(lParam);
                    Y_start = HIWORD(lParam);
                    counter++;
                } else if (counter == 3) {
                    cout << "End Point" << endl;
                    X_end = LOWORD(lParam);
                    Y_end = HIWORD(lParam);
                    CliipingLineWithCircle(hdc, X_start, Y_start, X_end, Y_end, X_left, Y_top, r, c);
//PointClipping(hdc,X_start,Y_start, X_left , Y_top , r , c);
                    counter = 2;
                }

            } else if (m == 53) //Point clipping (Rectangle).
            {
//c = RGB(0, 0,255);

                if (counter == 0) {
                    cout << "Counter " << counter << endl;
                    X_left = LOWORD(lParam);
                    Y_top = HIWORD(lParam);
                    counter++;
                } else if (counter == 1) {
                    cout << "Counter " << counter << endl;
                    X_right = LOWORD(lParam);
                    Y_bottom = HIWORD(lParam);
                    Rectangle(hdc, X_left, Y_top, X_right, Y_bottom);
                    counter++;
                } else if (counter == 2) {
// cout<<"Counter " <<counter<<endl;
                    X_start = LOWORD(lParam);
                    Y_start = HIWORD(lParam);
                    PointClipping(hdc, X_start, Y_start, X_left, Y_top, X_right, Y_bottom, c);
                    counter = 2;
                }
            } else if (m == 54) //Point clipping (Square)
            {
                if (Num_of_Points == 0) {
                    cout << "Counter " << Num_of_Points << endl;
                    X_left = LOWORD(lParam);
                    Y_top = HIWORD(lParam);
                    X_right = X_left + 200;
                    Y_bottom = Y_top + 200;
                    Rectangle(hdc, X_left, Y_top, X_right, Y_bottom);
                    cout << "Square Done.\n";
                    Num_of_Points++;
                } else if (Num_of_Points == 1) {
                    cout << "Point Clipping Done " << endl;
                    X_start = LOWORD(lParam);
                    Y_start = HIWORD(lParam);
                    PointClipping(hdc, X_start, Y_start, X_left, Y_top, X_right, Y_bottom, c);
                    Num_of_Points = 1;
                }
                Save_Point x("CohenSuth", X_start, Y_start, X_end, Y_end, X_left, Y_top, X_right, Y_bottom, 0, 0, 0);
                Arr_Save_Point.push_back(x);
            } else if (m == 55)   //Point clipping (Circle)
            {
                if (counter == 0) {
                    cout << "Center Of Circle. " << endl;
                    X_left = LOWORD(lParam);      //Center
                    Y_top = HIWORD(lParam);
                    counter++;
                } else if (counter == 1) {
                    cout << "Circle Window Done " << endl;
                    X_right = LOWORD(lParam);
                    Y_bottom = HIWORD(lParam);
                    r = sqrt(pow((X_right - X_left), 2) + pow((Y_bottom - Y_top), 2));
                    DrawCircleMidPoint(hdc, X_left, Y_top, r, RGB(0, 0, 0));     //radius.
                    counter++;
                } else if (counter == 2) {
                    cout << "Clipping Point (Circle) " << endl;
                    X_start = LOWORD(lParam);
                    Y_start = HIWORD(lParam);
//cout<<r1<<endl;
                    PointClipping(hdc, X_start, Y_start, X_left, Y_top, r, c);
                    counter = 2;
                }

            } else if (m == 56) //clipping Polygon (Rectangle)
            {

//P[5];

                if (counter == 0) {
                    cout << "Counter " << counter << endl;
                    X_left = LOWORD(lParam);
                    Y_top = HIWORD(lParam);
                    counter++;
                } else if (counter == 1) {
                    cout << "Counter " << counter << endl;
                    X_right = LOWORD(lParam);
                    Y_bottom = HIWORD(lParam);
                    Rectangle(hdc, X_left, Y_top, X_right, Y_bottom);
                    cout << "Clipping Window done " << counter << endl;
                    counter++;
                } else if (counter == 2) {
                    cout << "Counter " << counter << endl;
                    X_start = LOWORD(lParam);   //P[0]
                    Y_start = HIWORD(lParam);
                    P[0].x = X_start;
                    P[0].y = Y_start;
                    counter++;
                } else if (counter == 3) {
                    cout << "Counter " << counter << endl;
                    X_end = LOWORD(lParam);  //P[1]
                    Y_end = HIWORD(lParam);
                    P[1].x = X_end;
                    P[1].y = Y_end;
                    counter++;
                } else if (counter == 4) {
                    cout << "Counter " << counter << endl;
                    X_end = LOWORD(lParam);  //P[1]
                    Y_end = HIWORD(lParam);
                    P[2].x = X_end;
                    P[2].y = Y_end;
                    counter++;
                } else if (counter == 5) {
                    cout << "Counter " << counter << endl;
                    X_end = LOWORD(lParam);  //P[1]
                    Y_end = HIWORD(lParam);
                    P[3].x = X_end;
                    P[3].y = Y_end;
                    counter++;
                } else if (counter == 6) {

                    X_end = LOWORD(lParam);  //P[1]
                    Y_end = HIWORD(lParam);
                    P[4].x = X_end;
                    P[4].y = Y_end;
                    PolygonClip(hdc, P, 5, X_left, Y_top, X_right, Y_bottom);
                    counter = 0;
                }
            }

//////////////////////////////////////////////////////
            else if (m == 5 || m == 6) //Draw ellipse (Direct-polar)
            {
                if (m == 5) {
                    if (counter_ell == 0) {
                        x2 = LOWORD(lParam);
                        y2 = HIWORD(lParam);
                        counter_ell++;
                    } else if (counter_ell == 1) {
                        x3 = LOWORD(lParam);
                        y3 = HIWORD(lParam);
                        counter_ell++;
                    } else if (counter_ell == 2) {
                        r = Round(sqrt(pow(x2 - x1, 2.0) + pow(y2 - y_1, 2.0)));
                        r2 = Round(sqrt(pow(x3 - x1, 2.0) + pow(y3 - y_1, 2.0)));
                        cout << "Draw Direct Ellipse " << endl;
                        directellipse(hdc, x1, y_1, r, r2, c);
                        counter_ell = 0;
                    }

                    if (c == RGB(255, 0, 0)) {
                        Save_Point x("DDEllipse", x1, y_1, r, r2, 255, 0, 0, 'e');
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 255, 0)) {
                        Save_Point x("DDEllipse", x1, y_1, r, r2, 0, 255, 0, 'e');
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 0, 255)) {
                        Save_Point x("DDEllipse", x1, y_1, r, r2, 0, 0, 255, 'e');
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 0, 0)) {
                        Save_Point x("DDEllipse", x1, y_1, r, r2, 0, 0, 0, 'e');
                        Arr_Save_Point.push_back(x);
                    } else {
                        Save_Point x("DDEllipse", x1, y_1, r, r2, Rc, Gc, Bc, 'e');
                        Arr_Save_Point.push_back(x);
                    }

                    cout << "Draw Direct Ellipse Done !" << endl;
                } else if (m == 6)//Polar Ellipse
                {
                    if (counter_ell == 0) {
                        x2 = LOWORD(lParam);
                        y2 = HIWORD(lParam);
                        counter_ell++;
                    } else if (counter_ell == 1) {
                        x3 = LOWORD(lParam);
                        y3 = HIWORD(lParam);
                        counter_ell++;
                    } else if (counter_ell == 2) {
                        r = Round(sqrt(pow(x2 - x1, 2.0) + pow(y2 - y_1, 2.0)));
                        r2 = Round(sqrt(pow(x3 - x1, 2.0) + pow(y3 - y_1, 2.0)));
                        cout << "Draw polar Ellipse " << endl;
                        DrawEllipsePolar(hdc, x1, y_1, r, r2, c);
                        counter_ell = 0;
                    }

                    if (c == RGB(255, 0, 0)) {
                        Save_Point x("DPEllipse", x1, y_1, r, r2, 255, 0, 0, 'e');
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 255, 0)) {
                        Save_Point x("DPEllipse", x1, y_1, r, r2, 0, 255, 0, 'e');
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 0, 255)) {
                        Save_Point x("DPEllipse", x1, y_1, r, r2, 0, 0, 255, 'e');
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 0, 0)) {
                        Save_Point x("DPEllipse", x1, y_1, r, r2, 0, 0, 0, 'e');
                        Arr_Save_Point.push_back(x);
                    } else {
                        Save_Point x("DPEllipse", x1, y_1, r, r2, Rc, Gc, Bc, 'e');
                        Arr_Save_Point.push_back(x);
                    }
                    cout << "Draw Polar Ellipse Done !" << endl;

                }
            }
            if (ind == 2) {
                x3 = LOWORD(lParam);
                y3 = HIWORD(lParam);

                if (m == 61) {
                    Save_Point x("DRecBez", xz, y, x2, y, x2, y3, 0, 0, 0, 0, 0, 0);
                    Arr_Save_Point.push_back(x);

                    DrawRectangle(hdc, xz, y, x2, y, x2, y3, RGB(0, 0, 0), RGB(0, 0, 0));
                    cout << "Draw Rectangle and filling with Bezier curve (Black) is done. " << endl;
                    ind = 0;

                } else if (m == 62) {
                    Save_Point x("DRecBez", xz, y, x2, y, x2, y3, 0, 0, 0, 255, 0, 0);
                    Arr_Save_Point.push_back(x);
                    DrawRectangle(hdc, xz, y, x2, y, x2, y3, RGB(0, 0, 0), RGB(255, 0, 0));
                    cout << "Draw Rectangle and filling with Bezier curve (Red) is done. " << endl;
                    ind = 0;
                } else if (m == 64) {
                    Save_Point x("DRecBez", xz, y, x2, y, x2, y3, 0, 0, 0, 0, 255, 0);
                    Arr_Save_Point.push_back(x);
                    DrawRectangle(hdc, xz, y, x2, y, x2, y3, RGB(0, 0, 0), RGB(0, 255, 0));
                    cout << "Draw Rectangle and filling with Bezier curve (Green) is done. " << endl;
                    ind = 0;
                } else if (m == 63) {
                    Save_Point x("DRecBez", xz, y, x2, y, x2, y3, 0, 0, 0, 0, 0, 255);
                    Arr_Save_Point.push_back(x);
                    DrawRectangle(hdc, xz, y, x2, y, x2, y3, RGB(0, 0, 0), RGB(0, 0, 255));
                    cout << "Draw Rectangle and filling with Bezier curve (Blue) is done " << endl;
                    ind = 0;
                }


//ReleaseDC(hwnd, hdc);
            } else if (m == 65 || m == 66 || m == 67 || m == 68 || m == 69 || m == 70 || m == 71 || m == 72) {
                if (Num_of_Points2 == 0) {
                    x1 = LOWORD(lParam);
                    y_1 = HIWORD(lParam);
                    Num_of_Points2++;
                } else if (Num_of_Points2 == 1) {
                    x2 = LOWORD(lParam);
                    y2 = HIWORD(lParam);
                    Num_of_Points2++;
                } else if (Num_of_Points2 == 2) {
                    x3 = LOWORD(lParam);
                    y3 = HIWORD(lParam);
                    Num_of_Points2++;
                } else if (Num_of_Points2 == 3) {
                    x4 = LOWORD(lParam);
                    y4 = HIWORD(lParam);
                    Num_of_Points2++;
                } else if (Num_of_Points2 == 4) {
                    x5 = LOWORD(lParam);
                    y5 = HIWORD(lParam);
                    Polygon(hdc, P, 5);
                    Num_of_Points2++;
                } else if (Num_of_Points2 == 5) {
                    cout << x1 << y_1 << x2 << y2 << x5 << y5 << "\n";
                    P[0].x = x1;
                    P[0].y = y_1;
                    P[1].x = x2;
                    P[1].y = y2;
                    P[2].x = x3;
                    P[2].y = y3;
                    P[3].x = x4;
                    P[3].y = y4;
                    P[4].x = x5;
                    P[4].y = y5;

                    if (m == 65) {
                        Save_Point x("ConvexFill", x1, y_1, x2, y2, x3, y3, x4, y4, x5, y5, 0, 0, 0);
                        Arr_Save_Point.push_back(x);
//Polygon(hdc, P, 5);
                        ConvexFill(hdc, P, 5, RGB(0, 0, 0));
                        cout << "Convex filling (Black) is done " << endl;
                    } else if (m == 66) {
                        Save_Point x("ConvexFill", x1, y_1, x2, y2, x3, y3, x4, y4, x5, y5, 255, 0, 0);
                        Arr_Save_Point.push_back(x);
//Polygon(hdc, P, 5);
                        ConvexFill(hdc, P, 5, RGB(255, 0, 0));
                        cout << "Convex filling (Red) is done " << endl;
                    } else if (m == 67) {
                        Save_Point x("ConvexFill", x1, y_1, x2, y2, x3, y3, x4, y4, x5, y5, 0, 255, 0);
                        Arr_Save_Point.push_back(x);
//Polygon(hdc, P, 5);
                        ConvexFill(hdc, P, 5, RGB(0, 255, 0));
                        cout << "Convex filling (Green) is done " << endl;
                    } else if (m == 68) {
                        Save_Point x("ConvexFill", x1, y_1, x2, y2, x3, y3, x4, y4, x5, y5, 0, 0, 255);
                        Arr_Save_Point.push_back(x);
//Polygon(hdc, P, 5);
                        ConvexFill(hdc, P, 5, RGB(0, 0, 255));
                        cout << "Convex filling (Blue) is done " << endl;
                    } else if (m == 69) {
                        Save_Point x("GenFill", x1, y_1, x2, y2, x3, y3, x4, y4, x5, y5, 0, 0, 0);
                        Arr_Save_Point.push_back(x);
//Polygon(hdc, P, 5);
                        GeneralPolygonFill(hdc, P, 5, RGB(0, 0, 0));
                        cout << "Non Convex filling (Black) is done " << endl;
                    } else if (m == 70) {
                        Save_Point x("GenFill", x1, y_1, x2, y2, x3, y3, x4, y4, x5, y5, 255, 0, 0);
                        Arr_Save_Point.push_back(x);
//Polygon(hdc, P, 5);
                        GeneralPolygonFill(hdc, P, 5, RGB(255, 0, 0));
                        cout << "Non Convex filling (Red) is done " << endl;
                    } else if (m == 71) {
                        cout << "ll" << '\n';
                        Save_Point x("GenFill", x1, y_1, x2, y2, x2, y3, x4, y4, x5, y5, 0, 255, 0);
                        Arr_Save_Point.push_back(x);
//Polygon(hdc, P, 5);
                        GeneralPolygonFill(hdc, P, 5, RGB(0, 255, 0));
                        cout << "Non Convex filling (Green) is done " << endl;
                    } else if (m == 72) {
                        Save_Point x("GenFill", x1, y_1, x2, y2, x3, y3, x4, y4, x5, y5, 0, 0, 255);
                        Arr_Save_Point.push_back(x);
//Polygon(hdc, P, 5);
                        GeneralPolygonFill(hdc, P, 5, RGB(0, 0, 255));
                        cout << "Non Convex filling (Blue) is done " << endl;
                    }
                    Num_of_Points2 = 0;
                }
            }
            ReleaseDC(hwnd, hdc);
            break;

        case WM_RBUTTONDOWN:
            x2 = LOWORD(lParam);
            y2 = HIWORD(lParam);

            if (ind == 1 && (m == 57 || m == 58 || m == 59 || m == 60 || m == 61 || m == 62 || m == 63 || m == 64)) {
                x2 = LOWORD(lParam);
                y2 = HIWORD(lParam);
                if (m == 57) {

                    Save_Point x("DSqHer", x1, y_1, x2, y_1, 0, 0, 0, 0, 0, 0);
                    Arr_Save_Point.push_back(x);

                    DrawSquare(hdc, x1, y_1, x2, y_1, RGB(0, 0, 0), RGB(0, 0, 0));
                    cout << "Draw Square and filling with Hermite curve (Black) is done " << endl;
                    ind = 0;

                } else if (m == 58) {
                    Save_Point x("DSqHer", x1, y_1, x2, y_1, 0, 0, 0, 255, 0, 0);
                    Arr_Save_Point.push_back(x);
                    DrawSquare(hdc, x1, y_1, x2, y_1, RGB(0, 0, 0), RGB(255, 0, 0));
                    cout << "Draw Square and filling with Hermite curve (Red) is done " << endl;
                    ind = 0;
                } else if (m == 60) {
                    Save_Point x("DSqHer", x1, y_1, x2, y_1, 0, 0, 0, 0, 255, 0);
                    Arr_Save_Point.push_back(x);
                    DrawSquare(hdc, x1, y_1, x2, y_1, RGB(0, 0, 0), RGB(0, 255, 0));
                    cout << "Draw Square and filling with Hermite curve (Green) is done " << endl;
                    ind = 0;
                } else if (m == 59) {
                    Save_Point x("DSqHer", x1, y_1, x2, y_1, 0, 0, 0, 0, 0, 255);
                    Arr_Save_Point.push_back(x);
                    DrawSquare(hdc, x1, y_1, x2, y_1, RGB(0, 0, 0), RGB(0, 0, 255));
                    cout << "Draw Square and filling with Hermite curve (Blue) is done " << endl;
                    ind = 0;
                } else {
                    ind = 2;
                }
            }


            if (m == 7 || m == 8 || m == 9 || m == 26) {
                R = Round(std::sqrt(std::pow((x2 - x1), 2.0) + pow((y2 - y_1), 2.0)));
                if (m == 7) //Circle(Direct)
                {
                    cout << "Draw Direct Circle" << endl;
                    DrawCircle_Direct(hdc, x1, y_1, R, c);
                    if (c == RGB(255, 0, 0)) {
                        Save_Point x("DDCircle", x1, y_1, R, 255, 0, 0);
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 255, 0)) {

                        Save_Point x("DDCircle", x1, y_1, R, 0, 255, 0);
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 0, 255)) {

                        Save_Point x("DDCircle", x1, y_1, R, 0, 0, 255);
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 0, 0)) {

                        Save_Point x("DDCircle", x1, y_1, R, 0, 0, 0);
                        Arr_Save_Point.push_back(x);
                    } else {

                        Save_Point x("DDCircle", x1, y_1, R, Rc, Gc, Bc);
                        Arr_Save_Point.push_back(x);
                    }
                    cout << "Draw Direct Circle Done !" << endl;
                } else if (m == 8) //polar circle
                {
                    cout << "Draw Polar Circle" << endl;
                    DrawCircle_polar(hdc, x1, y_1, R, c);
                    if (c == RGB(255, 0, 0)) {
                        Save_Point x("DPCircle", x1, y_1, R, 255, 0, 0);
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 255, 0)) {

                        Save_Point x("DPCircle", x1, y_1, R, 0, 255, 0);
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 0, 255)) {

                        Save_Point x("DPCircle", x1, y_1, R, 0, 0, 255);
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 0, 0)) {

                        Save_Point x("DPCircle", x1, y_1, R, 0, 0, 0);
                        Arr_Save_Point.push_back(x);
                    } else {

                        Save_Point x("DPCircle", x1, y_1, R, Rc, Gc, Bc);
                        Arr_Save_Point.push_back(x);
                    }
                    cout << "Draw Polar Circle Done !" << endl;
                } else if (m == 9)//Circle(MidPooint)
                {
                    cout << "Draw MidPoint Circle " << endl;
                    Draw_Midpoint_circle(hdc, x1, y_1, R, c);
                    if (c == RGB(255, 0, 0)) {
                        Save_Point x("DMCircle", x1, y_1, R, 255, 0, 0);
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 255, 0)) {

                        Save_Point x("DMCircle", x1, y_1, R, 0, 255, 0);
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 0, 255)) {

                        Save_Point x("DMCircle", x1, y_1, R, 0, 0, 255);
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 0, 0)) {
                        Save_Point x("DMCircle", x1, y_1, R, 0, 0, 0);
                        Arr_Save_Point.push_back(x);
                    } else {
                        Save_Point x("DMCircle", x1, y_1, R, Rc, Gc, Bc);
                        Arr_Save_Point.push_back(x);
                    }
                    cout << "Draw MidPoint Circle Done !" << endl;
                } else if (m == 26) //Modified midpoint for circle
                {
                    cout << "Draw Modified MidPoint Circle " << endl;
                    DrawCircle_MidPoint_Modified(hdc, x1, y_1, R, c);
                    if (c == RGB(255, 0, 0)) {
                        Save_Point x("DMCircle", x1, y_1, R, 255, 0, 0);
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 255, 0)) {

                        Save_Point x("DMCircle", x1, y_1, R, 0, 255, 0);
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 0, 255)) {

                        Save_Point x("DMCircle", x1, y_1, R, 0, 0, 255);
                        Arr_Save_Point.push_back(x);
                    } else if (c == RGB(0, 0, 0)) {
                        Save_Point x("DMCircle", x1, y_1, R, 0, 0, 0);
                        Arr_Save_Point.push_back(x);
                    } else {
                        Save_Point x("DMCircle", x1, y_1, R, Rc, Gc, Bc);
                        Arr_Save_Point.push_back(x);
                    }
                    cout << "Draw Modified MidPoint Circle Done !" << endl;

                }

            } else if (m == 14 || m == 15 || m == 16 || m == 17) //Filling Circle
            {
//1-Circle
//2-color of Circle
//3-draw xc and yc click left button
//4-choose quarter
//5-choose filling color
//6-click right button
                x2 = LOWORD(lParam);
                y2 = HIWORD(lParam);
                R = Round(std::sqrt(std::pow(x2 - x1, 2.0) + pow(y2 - y_1, 2.0)));
                if (m == 14)//Filling Circle(Black color)
                {

                    cout << "Filling  Circle with black Color" << endl;
                    DrawCircleFilling(hdc, x1, y_1, R, quarter, RGB(0, 0, 0));
                    Save_Point x("DPCircleFilling", x1, y_1, R, quarter, 0, 0, 0, "e");
                    Arr_Save_Point.push_back(x);
                    cout << "Filling  Circle with black Color Done" << endl;


                } else if (m == 15)//Filling Circle(Red color)
                {
                    cout << "Filling  Circle with Red Color" << endl;
                    DrawCircleFilling(hdc, x1, y_1, R, quarter, RGB(255, 0, 0));
                    Save_Point x("DPCircleFilling", x1, y_1, R, quarter, 255, 0, 0, "e");
                    Arr_Save_Point.push_back(x);
                    cout << "Filling  Circle with Red Color Done" << endl;
                } else if (m == 16)//Filling Circle(Blue color)
                {
                    cout << "Filling  Circle with Blue Color" << endl;
                    DrawCircleFilling(hdc, x1, y_1, R, quarter, RGB(0, 0, 255));
                    Save_Point x("DPCircleFilling", x1, y_1, R, quarter, 0, 0, 255, "e");
                    Arr_Save_Point.push_back(x);
                    cout << "Filling  Circle with Blue Color Done" << endl;
                } else if (m == 17) //Filling Circle(Green color)
                {
                    cout << "Filling  Circle with Green Color" << endl;
                    DrawCircleFilling(hdc, x1, y_1, R, quarter, RGB(0, 255, 0));
                    Save_Point x("DPCircleFilling", x1, y_1, R, quarter, 0, 255, 0, "e");
                    Arr_Save_Point.push_back(x);
                    cout << "Filling  Circle with Green Color Done" << endl;
                }

            } else if (m == 27 || m == 28 || m == 29 || m == 30) //Filling Circle
            {
//1-Circle
//2-color of Circle
//3-draw xc and yc click left button
//4-choose quarter
//5-choose filling color
//6-click right button
                x2 = LOWORD(lParam);
                y2 = HIWORD(lParam);
                R = Round(std::sqrt(std::pow(x2 - x1, 2.0) + pow(y2 - y_1, 2.0)));
                if (m == 27)//Filling Circle(Black color)
                {

                    cout << "Filling  Circle with black Color" << endl;
                    DrawSolve(hdc, x1, y_1, R, quarter, RGB(0, 0, 0));
                    Save_Point x("DPCircleFilling", x1, y_1, R, quarter, 0, 0, 0, "e");
                    Arr_Save_Point.push_back(x);
                    cout << "Filling  Circle with black Color Done" << endl;


                } else if (m == 28)//Filling Circle(Red color)
                {
                    cout << "Filling  Circle with Red Color" << endl;
                    DrawSolve(hdc, x1, y_1, R, quarter, RGB(255, 0, 0));
                    Save_Point x("DPCircleFilling", x1, y_1, R, quarter, 255, 0, 0, "e");
                    Arr_Save_Point.push_back(x);
                    cout << "Filling  Circle with Red Color Done" << endl;
                } else if (m == 29)//Filling Circle(Blue color)
                {
                    cout << "Filling  Circle with Blue Color" << endl;
                    DrawSolve(hdc, x1, y_1, R, quarter, RGB(0, 0, 255));
                    Save_Point x("DPCircleFilling", x1, y_1, R, quarter, 0, 0, 255, "e");
                    Arr_Save_Point.push_back(x);
                    cout << "Filling  Circle with Blue Color Done" << endl;
                } else if (m == 30) //Filling Circle(Green color)
                {
                    cout << "Filling  Circle with Green Color" << endl;
                    DrawSolve(hdc, x1, y_1, R, quarter, RGB(0, 255, 0));
                    Save_Point x("DPCircleFilling", x1, y_1, R, quarter, 0, 255, 0, "e");
                    Arr_Save_Point.push_back(x);
                    cout << "Filling  Circle with Green Color Done" << endl;
                }

            } else if (m == 74) {
                int start_x = LOWORD(lParam);
                int start_y = HIWORD(lParam);
                cout << "Recursive Flodding Done!\n";
                int R = 0, G = 0, B = 0;
                if (c == RGB(0, 255, 0)) {
                    G = 255;
                }
                if (c == RGB(0, 0, 255)) {
                    B = 255;
                }
                if (c == RGB(255, 0, 0)) {
                    R = 255;
                }
                Save_Point x("FillRecursive", start_x, start_y, R, G, B, 216, 191, 216);
                Arr_Save_Point.push_back(x);
                FloodFill(hdc, start_x, start_y, c, RGB(216, 191, 216));

            } else if (m == 75) {
                int start_x = LOWORD(lParam);
                int start_y = HIWORD(lParam);
                cout << "Non-Recursive Flooding Done!\n";
                int R = 0, G = 0, B = 0;
                if (c == RGB(0, 255, 0)) {
                    G = 255;
                }
                if (c == RGB(0, 0, 255)) {
                    B = 255;
                }
                if (c == RGB(255, 0, 0)) {
                    R = 255;
                }
                Save_Point x("NonFillRecursive", start_x, start_y, R, G, B, 216, 191, 216);
                Arr_Save_Point.push_back(x);
                NRFloodFill(hdc, start_x, start_y, c, RGB(216, 191, 216));

            } else if (m == 76) {

                if (counter == 0) {
                    cout << "C:" << counter << endl;
                    X_start = LOWORD(lParam);
                    Y_start = HIWORD(lParam);
                    Vec[0].x = X_start;
                    Vec[0].y = Y_start;
                    counter++;
                } else if (counter == 1) {
                    cout << "C:" << counter << endl;
                    X_start = LOWORD(lParam);
                    Y_start = HIWORD(lParam);
                    Vec[1].x = X_start;
                    Vec[1].y = Y_start;
                    counter++;
                } else if (counter == 2) {
                    cout << "C:" << counter << endl;
                    X_start = LOWORD(lParam);
                    Y_start = HIWORD(lParam);
                    Vec[2].x = X_start;
                    Vec[2].y = Y_start;
                    counter++;
                } else if (counter == 3) {
                    cout << "C:" << counter << endl;

                    X_start = LOWORD(lParam);
                    Y_start = HIWORD(lParam);
                    Vec[3].x = X_start;
                    Vec[3].y = Y_start;
                    counter++;
                } else if (counter == 4) {
                    cout << "C:" << counter << endl;
                    X_start = LOWORD(lParam);
                    Y_start = HIWORD(lParam);
                    Vec[4].x = X_start;
                    Vec[4].y = Y_start;
                    counter = 0;
                    DrawCardinalSpline(hdc, Vec, 5, 0.5, 20);
                }
            } else if (m == 2) //DDA Line
            {
                cout << "Draw DDA line " << endl;
                DrawLine_DDA(hdc, x1, y_1, x2, y2, c);
                if (c == RGB(255, 0, 0)) {
                    Save_Point x("DDLine", x1, y_1, x2, y2, 255, 0, 0);
                    Arr_Save_Point.push_back(x);
                } else if (c == RGB(0, 255, 0)) {
                    Save_Point x("DDLine", x1, y_1, x2, y2, 0, 255, 0);
                    Arr_Save_Point.push_back(x);
                } else if (c == RGB(0, 0, 255)) {

                    Save_Point x("DDLine", x1, y_1, x2, y2, 0, 0, 255);
                    Arr_Save_Point.push_back(x);
                } else if (c == RGB(0, 0, 0)) {
                    Save_Point x("DDLine", x1, y_1, x2, y2, 0, 0, 0);
                    Arr_Save_Point.push_back(x);
                } else {
                    Save_Point x("DDLine", x1, y_1, x2, y2, Rc, Gc, Bc);
                    Arr_Save_Point.push_back(x);
                }
                cout << "Draw DDA line Done !" << endl;
            } else if (m == 3) //Midpoint Line
            {
                cout << "Draw MidPoint Line " << endl;
                BresenhamLine(hdc, x1, y_1, x2, y2, c);
                if (c == RGB(255, 0, 0)) {
                    Save_Point x("DMLine", x1, y_1, x2, y2, 255, 0, 0);
                    Arr_Save_Point.push_back(x);
                } else if (c == RGB(0, 255, 0)) {
                    Save_Point x("DMLine", x1, y_1, x2, y2, 0, 255, 0);
                    Arr_Save_Point.push_back(x);
                } else if (c == RGB(0, 0, 255)) {
                    Save_Point x("DMLine", x1, y_1, x2, y2, 0, 0, 255);
                    Arr_Save_Point.push_back(x);
                } else if (c == RGB(0, 0, 0)) {
                    Save_Point x("DMLine", x1, y_1, x2, y2, 0, 0, 0);
                    Arr_Save_Point.push_back(x);
                } else {
                    Save_Point x("DMLine", x1, y_1, x2, y2, Rc, Gc, Bc);
                    Arr_Save_Point.push_back(x);
                }
                cout << "Draw mid point line Done !" << endl;
            } else if (m == 4) //parametric Line
            {
                cout << "Draw parametric Line" << endl;
                DrawLine_Parametric(hdc, x1, y_1, x2, y2, c);
                if (c == RGB(255, 0, 0)) {

                    Save_Point x("DPLine", x1, y_1, x2, y2, 255, 0, 0);
                    Arr_Save_Point.push_back(x);
                } else if (c == RGB(0, 255, 0)) {

                    Save_Point x("DPLine", x1, y_1, x2, y2, 0, 255, 0);
                    Arr_Save_Point.push_back(x);
                } else if (c == RGB(0, 0, 255)) {

                    Save_Point x("DPLine", x1, y_1, x2, y2, 0, 0, 255);
                    Arr_Save_Point.push_back(x);
                } else if (c == RGB(0, 0, 0)) {

                    Save_Point x("DPLine", x1, y_1, x2, y2, 0, 0, 0);
                    Arr_Save_Point.push_back(x);
                } else {

                    Save_Point x("DPLine", x1, y_1, x2, y2, Rc, Gc, Bc);
                    Arr_Save_Point.push_back(x);
                }
                cout << "Draw parametric line Done !" << endl;
            }
            break;
        case WM_LBUTTONDBLCLK:


            break;

        case WM_CREATE:
            addMenu(hwnd);
            break;
        case WM_COMMAND: {
            switch (wParam) {
                case (0):
                    cout << "Saving Process\n\n";
                    Save();
                    cout << "Save Done ! \n\n";
                    break;
                case (1):
                    cout << "Load process\n\n";
                    Load(hdc);
                    cout << "load Done !\n\n";
                    break;
                case (2):
                    m = 2;
                    cout << "Draw Line by DDA\n\n LeftClick then RightClick \n";
                    break;
                case (3):
                    m = 3;
                    cout << "Draw Line by Mid Point\n\n";
                    break;
                case (4):
                    m = 4;
                    cout << "Draw Line by Parametric\n\n";
                    break;
                case (5):
                    m = 5;
                    cout << " Draw Direct ellipse \n\n FirstClick for 1st R then 2nd for 2nd R then Center point\n";
                    break;
                case (6):
                    m = 6;
                    cout << "Draw polar ellipse \n\n FirstClick for 1st R then 2nd for 2nd R then Center point\n";
                    break;
                case (7):
                    m = 7;
                    cout << "Draw direct circle \n\n FirstClick for Center then R \n";
                    break;
                case (8):
                    m = 8;
                    cout << "Draw Polar circle \n\n  FirstClick for Center then R \n";
                    break;
                case (9):
                    m = 9;
                    cout << "Draw Mid Point circle \n\n  FirstClick for Center then R \n";
                    break;
                case (26):
                    m = 26;
                    cout << "Draw circle by Modified Mid point \n\n  FirstClick for Center then R \n";
                    break;
                case (10):
                    cout << "Draw by black color \n\n";
                    c = 0;
                    break;
                case (11):
                    cout << "Draw by red color \n\n";
                    c = RGB(255, 0, 0);
                    break;
                case (12):
                    cout << "Draw by blue color \n\n";
                    c = RGB(0, 0, 255);
                    break;
                case (13):
                    cout << "Draw by green color \n\n";
                    c = RGB(0, 255, 0);
                    break;
                case (14):
                    m = 14;
                    cout << "Fill circle by the Black color\n\n";
                    cF = 0;
                    break;
                case (15):
                    m = 15;
                    cout << "Fill circle by the Red color\n\n";
                    cF = RGB(255, 0, 0);
                    break;
                case (16):
                    m = 16;
                    cout << "Fill circle by  the blue color \n\n";
                    cF = RGB(0, 0, 255);
                    break;
                case (17):
                    m = 17;
                    cout << "Fill circle by the green color \n\n";
                    cF = RGB(0, 255, 0);
                    break;
                case (18):
                    m = 18;
                    cout << "Filling Circle with First Quarter \n\n";
                    quarter = 1;
                    break;
                case (19):
                    m = 19;
                    cout << "Filling circle with second Quarter \n\n";
                    quarter = 2;
                    break;
                case (20):
                    m = 20;
                    cout << "Filling circle with third Quarter \n\n";
                    quarter = 3;
                    break;
                case (21):
                    m = 21;
                    cout << "Filling circle with fourth Quarter \n\n";
                    quarter = 4;
                    break;
                case (22):
                    m = 22;
                    cout << "line after clipping Using Rectangle as Window\n\n";
                    break;
                case (23):
                    m = 23;
                    cout << "clipping point Using Rectangle as Window \n\n";
                    break;
                case (50):
                    m = 50;
                    cout << "clipping Polygon Using Rectangle as Window \n\n";
                    break;
                case (51):
                    m = 51;
                    cout << "clipping Line Using Circle as Window \n\n";
                    break;

                case (52):
                    m = 52;
                    cout << "clipping Point Using Circle as Window \n\n";
                    break;
                case (25):
                    cout << "Clean the window \n";
                    Clear();
                    InvalidateRect(hwnd, nullptr, TRUE);
                    break;
                case (27):
                    m = 27;
                    cout << "Filling Circle  with other circle in First Quarter \n\n";
                    quarter = 1;
                    break;
                case (28):
                    m = 28;
                    cout << "Filling Circle  with other circle in second Quarter \n\n";
                    quarter = 2;
                    break;
                case (29):
                    m = 29;
                    cout << "Filling Circle  with other circle in third Quarter \n\n";
                    quarter = 3;
                    break;
                case (30):
                    m = 30;
                    cout << "Filling Circle  with other circle in fourth Quarter \n\n";
                    quarter = 4;
                    break;
                case (31):
                    curs = IDC_HAND;
                    SetCursor(LoadCursor(nullptr, curs));
                    cout << "Setting Courser\n";
                    m = 31;
                    break;
                case (32):
                    curs = IDC_ARROW;
                    m = 32;
                    break;
                case (33):
                    curs = IDC_CROSS;
                    m = 33;
                    break;
                case (34):
                    curs = IDC_HELP;
                    m = 34;
                    break;
                case (35):
                    curs = IDC_UPARROW;
                    m = 35;
                    break;
                case (36):
                    curs = IDC_IBEAM;
                    m = 36;
                    break;
                case (37):
                    curs = IDC_WAIT;
                    m = 37;
                    break;
                case (57):
                    m = 57;
                    cout << "Draw square and Filling Black (Hermite) \n\n";
                    break;
                case (58):
                    m = 58;
                    cout << "Draw square and Filling Red (Hermite) \n\n";
                    break;
                case (59):
                    m = 59;
                    cout << "Draw square and Filling Blue (Hermite) \n\n";
                    break;
                case (60):
                    m = 60;
                    cout << "Draw square and Filling Green (Hermite) \n\n";
                    break;

                case (61):
                    m = 61;
                    cout << "Draw Rectangle and Filling Black (Bezier) \n\n";
                    break;
                case (62):
                    m = 62;
                    cout << "Draw Rectangle and Filling Red (Bezier) \n\n";
                    break;
                case (63):
                    m = 63;
                    cout << "Draw Rectangle and Filling Blue (Bezier) \n\n";
                    break;
                case (64):
                    m = 64;
                    cout << "Draw Rectangle and Filling Green (Bezier) \n\n";
                    break;

                case (65):
                    m = 65;
                    cout << "Convex filling Polygon Black \n\n";
                    break;
                case (66):
                    m = 66;
                    cout << "Convex filling Polygon Red \n\n";
                    break;
                case (67):
                    m = 67;
                    cout << "Convex filling Polygon Green \n\n";
                    break;
                case (68):
                    m = 68;
                    cout << "Convex filling Polygon Blue \n\n";
                    break;

                case (69):
                    m = 69;
                    cout << "Non Convex filling Polygon Black \n\n";
                    break;
                case (70):
                    m = 70;
                    cout << "Non Convex filling Polygon Red \n\n";
                    break;
                case (71):
                    m = 71;
                    cout << "Non Convex filling Polygon Green \n\n";
                    break;
                case (72):
                    m = 72;
                    cout << "Non Convex filling Polygon Blue \n\n";
                    break;

                case (74):
                    m = 74;
                    cout << "Flood Fill Recursive\n\n";
                    break;

                case (75):
                    m = 75;
                    cout << "Flood Fill Non Recursive\n\n";
                    break;

                case (76):
                    m = 76;
                    cout << " Cardinal Spline Curve\n\n";
                    break;
                default:
                    m = LOWORD(wParam);
                    break;

            }
        }
            break;
        case WM_SETCURSOR:
            if (LOWORD(lParam) == HTCLIENT) {
                SetCursor(LoadCursor(nullptr, curs));
                return TRUE;
            }
            break;
        case WM_DESTROY:
            PostQuitMessage(0);       /* send a WM_QUIT to the message queue */
            break;
        default:                      /* for messages that we don't deal with */
            return DefWindowProc(hwnd, message, wParam, lParam);
    }
    return 0;
}