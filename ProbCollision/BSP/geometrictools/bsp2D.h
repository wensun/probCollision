#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <deque>
#include <math.h>
#include <float.h>

class Edge
{
public:
    int i0, i1;

	Edge();
    Edge(int i0, int i1);
    bool operator< (const Edge& edge) const;
};


#include "Wm5Vector2.h"
using namespace Wm5;

// vertices
typedef std::map<Vector2d,int> VMap;
typedef std::map<Vector2d,int>::iterator VIterator;
typedef std::map<Vector2d,int>::const_iterator VCIterator;
typedef std::vector<Vector2d> VArray;

// edges
typedef std::map<Edge2,int> EMap;
typedef std::map<Edge2,int>::iterator EIterator;
typedef std::map<Edge2,int>::const_iterator ECIterator;
typedef std::vector<Edge2> EArray;

class BspPolygon;
class BspTree
{
public:
    // Construction and destruction.
    BspTree2 (BspPolygon2& polygon, const EArray& edges);
    ~BspTree2 ();

    BspTree2* GetCopy () const;

    // Polygon Boolean operation support.
    void Negate ();
    void GetPartition (const BspPolygon2& polygon, const Vector2d& v0,
        const Vector2d& v1, BspPolygon2& pos, BspPolygon2& neg,
        BspPolygon2& coSame, BspPolygon2& coDiff) const;

    // Point-in-polygon support (-1 outside, 0 on polygon, +1 inside).
    int PointLocation (const BspPolygon2& polygon, const Vector2d& vertex)
        const;

    void Indent (std::ofstream& outFile, int numSpaces);
    void Print (std::ofstream& outFile, int level, char type);

private:
    BspTree2 ()
    {
        // support for get copy
    }

    BspTree2 (const BspTree2&)
    {
        // not supported
    }

    BspTree2& operator= (const BspTree2&)
    {
        // not supported
        return *this;
    }

    enum
    {
        TRANSVERSE_POSITIVE,
        TRANSVERSE_NEGATIVE,
        ALL_POSITIVE,
        ALL_NEGATIVE,
        COINCIDENT
    };

    int Classify (const Vector2d& end0, const Vector2d& end1,
        const Vector2d& v0, const Vector2d& v1, Vector2d& intr) const;

    void GetPosPartition (const BspPolygon2& polygon, const Vector2d& v0,
        const Vector2d& v1, BspPolygon2& pos, BspPolygon2& neg,
        BspPolygon2& coSame,  BspPolygon2& coDiff) const;

    void GetNegPartition (const BspPolygon2& polygon, const Vector2d& v0,
        const Vector2d& v1, BspPolygon2& pos, BspPolygon2& neg,
        BspPolygon2& coSame, BspPolygon2& coDiff) const;

    class Interval
    {
    public:
        Interval (double t0, double t1, bool sameDir, bool touching)
        {
            T0 = t0;
            T1 = t1;
            SameDir = sameDir;
            Touching = touching;
        }

        double T0, T1;
        bool SameDir, Touching;
    };

    void GetCoPartition (const BspPolygon2& polygon, const Vector2d& v0,
        const Vector2d& v1, BspPolygon2& pos, BspPolygon2& neg,
        BspPolygon2& coSame, BspPolygon2& coDiff) const;

    // point-in-polygon support
    int Classify (const Vector2d& end0, const Vector2d& end1,
        const Vector2d& vertex) const;
    int CoPointLocation (const BspPolygon2& polygon, const Vector2d& vertex)
        const;

    EArray mCoincident;
    BspTree2* mPosChild;
    BspTree2* mNegChild;
};

#endif










// Geometric Tools, LLC
// Copyright (c) 1998-2011
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)

#ifndef BSPPOLYGON2_H
#define BSPPOLYGON2_H

#include "Types2.h"

class BspTree2;

class BspPolygon2
{
public:
    // Construction and destruction.
    BspPolygon2 ();
    BspPolygon2 (const BspPolygon2& polygon);
    ~BspPolygon2 ();

    // Support for deferred construction.
    int InsertVertex (const Vector2d& vertex);
    int InsertEdge (const Edge2& edge);
    void Finalize ();

    // Assignment.
    BspPolygon2& operator= (const BspPolygon2& polygon);

    // Member access.
    int GetNumVertices () const;
    bool GetVertex (int i, Vector2d& vertex) const;
    int GetNumEdges () const;
    bool GetEdge (int i, Edge2& edge) const;

    // negation
    BspPolygon2 operator~ () const;

    // intersection
    BspPolygon2 operator& (const BspPolygon2& polygon) const;

    // union
    BspPolygon2 operator| (const BspPolygon2& polygon) const;

    // difference
    BspPolygon2 operator- (const BspPolygon2& polygon) const;

    // exclusive or
    BspPolygon2 operator^ (const BspPolygon2& polygon) const;

    // point location (-1 inside, 0 on polygon, 1 outside)
    int PointLocation (const Vector2d& vertex) const;

    void Print (const char* filename) const;

protected:
    void SplitEdge (int v0, int v1, int vmid);
    void GetInsideEdgesFrom (const BspPolygon2& polygon, BspPolygon2& inside)
        const;

    // vertices
    VMap mVMap;
    VArray mVArray;

    // edges
    EMap mEMap;
    EArray mEArray;

    friend class BspTree2;
    BspTree2* mTree;
};


// Geometric Tools, LLC
// Copyright (c) 1998-2011
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)

#ifndef BOOLEAN2D_H
#define BOOLEAN2D_H

#include "Wm5WindowApplication2.h"
#include "BspPolygon2.h"
using namespace Wm5;

class Boolean2D : public WindowApplication2
{
    WM5_DECLARE_INITIALIZE;
    WM5_DECLARE_TERMINATE;

public:
    Boolean2D ();

    virtual bool OnInitialize ();
    virtual void OnTerminate ();
    virtual void OnDisplay ();
    virtual bool OnKeyDown (unsigned char key, int x, int y);

protected:
    BspPolygon2* ConstructInvertedEll ();
    BspPolygon2* ConstructPentagon ();
    BspPolygon2* ConstructSquare ();
    BspPolygon2* ConstructSShape ();
    BspPolygon2* ConstructPolyWithHoles ();

    void DoBoolean ();
    void DrawPolySolid (BspPolygon2& polygon, ColorRGB color);

    BspPolygon2 mIntersection, mUnion, mDiff01, mDiff10, mXor;
    BspPolygon2* mPoly0;
    BspPolygon2* mPoly1;
    BspPolygon2* mActive;
    int mChoice;
    int mSize;
};

WM5_REGISTER_INITIALIZE(Boolean2D);
WM5_REGISTER_TERMINATE(Boolean2D);

#endif


