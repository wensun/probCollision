// Geometric Tools, LLC
// Copyright (c) 1998-2011
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)

#include "Edge2.h"

//----------------------------------------------------------------------------
Edge2::Edge2 ()
{
    I0 = -1;
    I1 = -1;
}
//----------------------------------------------------------------------------
Edge2::Edge2 (int i0, int i1)
{
    I0 = i0;
    I1 = i1;
}
//----------------------------------------------------------------------------
bool Edge2::operator< (const Edge2& edge) const
{
    if (I1 < edge.I1)
    {
        return true;
    }

    if (I1 > edge.I1)
    {
        return false;
    }

    return I0 < edge.I0;
}
//----------------------------------------------------------------------------





// Geometric Tools, LLC
// Copyright (c) 1998-2011
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)

#include "BspTree2.h"
#include "BspPolygon2.h"
#include "Wm5Memory.h"

//----------------------------------------------------------------------------
BspTree2::BspTree2 (BspPolygon2& polygon, const EArray& edges)
{
    assertion(edges.size() > 0, "Invalid input.\n");

    // Construct splitting line from first edge.
    Vector2d end0 = polygon.mVArray[edges[0].I0];
    Vector2d end1 = polygon.mVArray[edges[0].I1];

    // Add edge to coincident list.
    mCoincident.push_back(edges[0]);

    // Split remaining edges.
    EArray posArray, negArray;
    int imax = (int)edges.size();
    for (int i = 1; i < imax; ++i)
    {
        int v0 = edges[i].I0;
        int v1 = edges[i].I1;
        Vector2d vertex0 = polygon.mVArray[v0];
        Vector2d vertex1 = polygon.mVArray[v1];

        Vector2d intr;
        int vmid;

        switch (Classify(end0, end1, vertex0, vertex1, intr))
        {
            case TRANSVERSE_POSITIVE:
                // modify edge <V0,V1> to <V0,I> and add new edge <I,V1>
                vmid = polygon.InsertVertex(intr);
                polygon.SplitEdge(v0, v1, vmid);
                posArray.push_back(Edge2(vmid, v1));
                negArray.push_back(Edge2(v0, vmid));
                break;
            case TRANSVERSE_NEGATIVE:
                // modify edge <V0,V1> to <V0,I> and add new edge <I,V1>
                vmid = polygon.InsertVertex(intr);
                polygon.SplitEdge(v0, v1, vmid);
                posArray.push_back(Edge2(v0, vmid));
                negArray.push_back(Edge2(vmid, v1));
                break;
            case ALL_POSITIVE:
                posArray.push_back(edges[i]);
                break;
            case ALL_NEGATIVE:
                negArray.push_back(edges[i]);
                break;
            default:  // COINCIDENT
                mCoincident.push_back(edges[i]);
                break;
        }
    }

    if (posArray.size() > 0)
    {
        mPosChild = new0 BspTree2(polygon, posArray);
    }
    else
    {
        mPosChild = 0;
    }

    if (negArray.size() > 0)
    {
        mNegChild = new0 BspTree2(polygon, negArray);
    }
    else
    {
        mNegChild = 0;
    }
}
//----------------------------------------------------------------------------
BspTree2::~BspTree2 ()
{
    delete0(mPosChild);
    delete0(mNegChild);
}
//----------------------------------------------------------------------------
BspTree2* BspTree2::GetCopy () const
{
    BspTree2* tree = new0 BspTree2();

    tree->mCoincident = mCoincident;

    if (mPosChild)
    {
        tree->mPosChild = mPosChild->GetCopy();
    }
    else
    {
        tree->mPosChild = 0;
    }

    if (mNegChild)
    {
        tree->mNegChild = mNegChild->GetCopy();
    }
    else
    {
        tree->mNegChild = 0;
    }

    return tree;
}
//----------------------------------------------------------------------------
void BspTree2::Negate ()
{
    // Reverse coincident edge directions.
    const int numEdges = (int)mCoincident.size();
    for (int i = 0; i < numEdges; ++i)
    {
        Edge2& edge = mCoincident[i];
        int save = edge.I0;
        edge.I0 = edge.I1;
        edge.I1 = save;
    }

    // Swap positive and negative subtrees.
    BspTree2* save = mPosChild;
    mPosChild = mNegChild;
    mNegChild = save;

    if (mPosChild)
    {
        mPosChild->Negate();
    }

    if (mNegChild)
    {
        mNegChild->Negate();
    }
}
//----------------------------------------------------------------------------
int BspTree2::Classify (const Vector2d& end0, const Vector2d& end1,
    const Vector2d& v0, const Vector2d& v1, Vector2d& intr) const
{
    // For numerical round-off error handling.
    const double epsilon0 = 0.00001;
    const double epsilon1 = 0.99999;

    Vector2d dir = end1 - end0;
    Vector2d nor = dir.Perp();
    Vector2d diff0 = v0 - end0;
    Vector2d diff1 = v1 - end0;

    double d0 = nor.Dot(diff0);
    double d1 = nor.Dot(diff1);

    if (d0*d1 < 0.0)
    {
        // Edge <V0,V1> transversely crosses line.  Compute point of
        // intersection I = V0 + t*(V1 - V0).
        double t = d0/(d0 - d1);
        if (t > epsilon0)
        {
            if (t < epsilon1)
            {
                intr = v0 + t*(v1 - v0);
                if (d1 > 0.0)
                {
                    return TRANSVERSE_POSITIVE;
                }
                else
                {
                    return TRANSVERSE_NEGATIVE;
                }
            }
            else
            {
                // T is effectively 1 (numerical round-off issue), so
                // set d1 = 0 and go on to other cases.
                d1 = 0.0;
            }
        }
        else
        {
            // T is effectively 0 (numerical round-off issue), so
            // set d0 = 0 and go on to other cases.
            d0 = 0.0;
        }
    }

    if (d0 > 0.0 || d1 > 0.0)
    {
        // edge on positive side of line
        return ALL_POSITIVE;
    }

    if (d0 < 0.0 || d1 < 0.0)
    {
        // edge on negative side of line
        return ALL_NEGATIVE;
    }

    return COINCIDENT;
}
//----------------------------------------------------------------------------
void BspTree2::GetPosPartition (const BspPolygon2& polygon,
    const Vector2d& v0, const Vector2d& v1, BspPolygon2& pos,
    BspPolygon2& neg, BspPolygon2& coSame, BspPolygon2& coDiff) const
{
    if (mPosChild)
    {
        mPosChild->GetPartition(polygon, v0, v1, pos, neg, coSame, coDiff);
    }
    else
    {
        int i0 = pos.InsertVertex(v0);
        int i1 = pos.InsertVertex(v1);
        pos.InsertEdge(Edge2(i0, i1));
    }
}
//----------------------------------------------------------------------------
void BspTree2::GetNegPartition (const BspPolygon2& polygon,
    const Vector2d& v0, const Vector2d& v1, BspPolygon2& pos,
    BspPolygon2& neg, BspPolygon2& coSame, BspPolygon2& coDiff) const
{
    if (mNegChild)
    {
        mNegChild->GetPartition(polygon, v0, v1, pos, neg, coSame, coDiff);
    }
    else
    {
        int i0 = neg.InsertVertex(v0);
        int i1 = neg.InsertVertex(v1);
        neg.InsertEdge(Edge2(i0, i1));
    }
}
//----------------------------------------------------------------------------
void BspTree2::GetCoPartition (const BspPolygon2& polygon, const Vector2d& v0,
    const Vector2d& v1, BspPolygon2& pos, BspPolygon2& neg,
    BspPolygon2& coSame, BspPolygon2& coDiff) const
{
    const double epsilon = 0.00001;

    // Segment the line containing V0 and V1 by the coincident intervals that
    // intersect <V0,V1>.
    Vector2d dir = v1 - v0;
    double tmax = dir.Dot(dir);

    Vector2d end0, end1;
    double t0, t1;
    bool sameDir;

    std::list<Interval> intervalList;
    std::list<Interval>::iterator iter;

    const int numEdges = (int)mCoincident.size();
    for (int i = 0; i < numEdges; ++i)
    {
        end0 = polygon.mVArray[mCoincident[i].I0];
        end1 = polygon.mVArray[mCoincident[i].I1];

        t0 = dir.Dot(end0 - v0);
        if (Mathd::FAbs(t0) <= epsilon)
        {
            t0 = 0.0;
        }
        else if (Mathd::FAbs(t0 - tmax) <= epsilon)
        {
            t0 = tmax;
        }

        t1 = dir.Dot(end1 - v0);
        if (Mathd::FAbs(t1) <= epsilon)
        {
            t1 = 0.0;
        }
        else if (Mathd::FAbs(t1 - tmax) <= epsilon)
        {
            t1 = tmax;
        }

        sameDir = (t1 > t0);
        if (!sameDir)
        {
            double save = t0;
            t0 = t1;
            t1 = save;
        }

        if (t1 > 0.0 && t0 < tmax)
        {
            if (intervalList.empty())
            {
                intervalList.push_front(Interval(t0, t1, sameDir, true));
            }
            else
            {
                iter = intervalList.begin();
                for (/**/; iter != intervalList.end(); ++iter)
                {
                    if (Mathd::FAbs(t1 - iter->T0) <= epsilon)
                    {
                        t1 = iter->T0;
                    }

                    if (t1 <= iter->T0)
                    {
                        // [t0,t1] is on the left of [I.t0,I.t1]
                        intervalList.insert(iter,
                            Interval(t0, t1, sameDir, true));
                        break;
                    }

                    // Theoretically, the intervals are disjoint or intersect
                    // only at an end point.  The assert makes sure that
                    // [t0,t1] is to the right of [I.t0,I.t1].
                    if (Mathd::FAbs(t0 - iter->T1) <= epsilon)
                    {
                        t0 = iter->T1;
                    }

                    assertion(t0 >= iter->T1, "Invalid ordering.\n");

                    std::list<Interval>::iterator last = intervalList.end();
                    --last;
                    if (iter == last)
                    {
                        intervalList.push_back(Interval(t0, t1, sameDir,
                            true));
                        break;
                    }
                }
            }
        }
    }

    if (intervalList.empty())
    {
        GetPosPartition(polygon, v0, v1, pos, neg, coSame, coDiff);
        GetNegPartition(polygon, v0, v1, pos, neg, coSame, coDiff);
        return;
    }

    // Insert outside intervals between the touching intervals.  It is
    // possible that two touching intervals are adjacent, so this is not just
    // a simple alternation of touching and outside intervals.
    Interval& front = intervalList.front();
    if (front.T0 > 0.0)
    {
        intervalList.push_front(Interval(0.0, front.T0, front.SameDir,
            false));
    }
    else
    {
        front.T0 = 0.0;
    }

    Interval& back = intervalList.back();
    if (back.T1 < tmax)
    {
        intervalList.push_back(Interval(back.T1, tmax, back.SameDir, false));
    }
    else
    {
        back.T1 = tmax;
    }

    std::list<Interval>::iterator iter0 = intervalList.begin();
    std::list<Interval>::iterator iter1 = intervalList.begin();
    for (++iter1; iter1 != intervalList.end(); ++iter0, ++iter1)
    {
        t0 = iter0->T1;
        t1 = iter1->T0;
        if (t1 - t0 > epsilon)
        {
            iter0 = intervalList.insert(iter1, Interval(t0, t1, true, false));
        }
    }

    // Process the segmentation.
    double invTMax = 1.0/tmax;
    t0 = intervalList.front().T0*invTMax;
    end1 = v0 + (intervalList.front().T0*invTMax)*dir;
    iter = intervalList.begin();
    for (/**/; iter != intervalList.end(); ++iter)
    {
        end0 = end1;
        t1 = iter->T1*invTMax;
        end1 = v0 + (iter->T1*invTMax)*dir;

        if (iter->Touching)
        {
            Edge2 edge;
            if (iter->SameDir)
            {
                edge.I0 = coSame.InsertVertex(end0);
                edge.I1 = coSame.InsertVertex(end1);
                if (edge.I0 != edge.I1)
                {
                    coSame.InsertEdge(edge);
                }
            }
            else
            {
                edge.I0 = coDiff.InsertVertex(end1);
                edge.I1 = coDiff.InsertVertex(end0);
                if (edge.I0 != edge.I1)
                {
                    coDiff.InsertEdge(edge);
                }
            }
        }
        else
        {
            GetPosPartition(polygon, end0, end1, pos, neg, coSame, coDiff);
            GetNegPartition(polygon, end0, end1, pos, neg, coSame, coDiff);
        }
    }
}
//----------------------------------------------------------------------------
void BspTree2::GetPartition (const BspPolygon2& polygon, const Vector2d& v0, const Vector2d& v1, BspPolygon2& pos, BspPolygon2& neg, BspPolygon2& coSame, BspPolygon2& coDiff) const
{
    // Construct splitting line from first coincident edge.
    Vector2d end0 = polygon.mVArray[mCoincident[0].I0];
    Vector2d end1 = polygon.mVArray[mCoincident[0].I1];

    Vector2d intr;

    switch (Classify(end0, end1, v0, v1, intr))
    {
    case TRANSVERSE_POSITIVE:
        GetPosPartition(polygon, intr, v1, pos, neg, coSame, coDiff);
        GetNegPartition(polygon, v0, intr, pos, neg, coSame, coDiff);
        break;
    case TRANSVERSE_NEGATIVE:
        GetPosPartition(polygon, v0, intr, pos, neg, coSame, coDiff);
        GetNegPartition(polygon, intr, v1, pos, neg, coSame, coDiff);
        break;
    case ALL_POSITIVE:
        GetPosPartition(polygon, v0, v1, pos, neg, coSame, coDiff);
        break;
    case ALL_NEGATIVE:
        GetNegPartition(polygon, v0, v1, pos, neg, coSame, coDiff);
        break;
    default:  // COINCIDENT
        GetCoPartition(polygon, v0, v1, pos, neg, coSame, coDiff);
        break;
    }
}
//----------------------------------------------------------------------------
int BspTree2::Classify (const Vector2d& end0, const Vector2d& end1, const Vector2d& vertex) const
{
    // For numerical round-off error handling.
    const double epsilon = 0.00001;

    Vector2d dir = end1 - end0;
    Vector2d nor = dir.Perp();
    Vector2d diff = vertex - end0;
    double c = nor.Dot(diff);

    if (c > epsilon)
    {
        return ALL_POSITIVE;
    }

    if (c < -epsilon)
    {
        return ALL_NEGATIVE;
    }

    return COINCIDENT;
}
//----------------------------------------------------------------------------
int BspTree2::CoPointLocation (const BspPolygon2& polygon, const Vector2d& vertex) const
{
    // For numerical round-off error handling.
    const double epsilon = 0.00001;

    const int numEdges = (int)mCoincident.size();
    for (int i = 0; i < numEdges; ++i)
    {
        Vector2d end0 = polygon.mVArray[mCoincident[i].I0];
        Vector2d end1 = polygon.mVArray[mCoincident[i].I1];
        Vector2d dir = end1 - end0;
        Vector2d diff = vertex - end0;
        double tmax = dir.Dot(dir);
        double t = dir.Dot(diff);

        if (-epsilon <= t && t <= tmax + epsilon)
        {
            return 0;
        }
    }

    // Does not matter which subtree you use.
    if (mPosChild)
    {
        return mPosChild->PointLocation(polygon, vertex);
    }

    if (mNegChild)
    {
        return mNegChild->PointLocation(polygon, vertex);
    }

    return 0;
}
//----------------------------------------------------------------------------
int BspTree2::PointLocation (const BspPolygon2& polygon, const Vector2d& vertex) const
{
    // Construct splitting line from first coincident edge.
    Vector2d end0 = polygon.mVArray[mCoincident[0].I0];
    Vector2d end1 = polygon.mVArray[mCoincident[0].I1];

    switch (Classify(end0, end1, vertex))
    {
    case ALL_POSITIVE:
        if (mPosChild)
        {
            return mPosChild->PointLocation(polygon, vertex);
        }
        else
        {
            return 1;
        }
    case ALL_NEGATIVE:
        if (mNegChild)
        {
            return mNegChild->PointLocation(polygon, vertex);
        }
        else
        {
            return -1;
        }
    default:  // COINCIDENT
        return CoPointLocation(polygon, vertex);
    }
}
//----------------------------------------------------------------------------
void BspTree2::Indent (std::ofstream& outFile, int numSpaces)
{
    for (int i = 0; i < numSpaces; ++i)
    {
        outFile << ' ';
    }
}
//----------------------------------------------------------------------------
void BspTree2::Print (std::ofstream& outFile, int level, char type)
{
    const int numEdges = (int)mCoincident.size();
    for (int i = 0; i < numEdges; ++i)
    {
        Indent(outFile, 4*level);
        outFile << type << " <" << mCoincident[i].I0 << ',' <<
            mCoincident[i].I1 << ">" << std::endl;
    }

    if (mPosChild)
    {
        mPosChild->Print(outFile, level + 1, 'p');
    }

    if (mNegChild)
    {
        mNegChild->Print(outFile, level + 1, 'n');
    }
}
//----------------------------------------------------------------------------



// Geometric Tools, LLC
// Copyright (c) 1998-2011
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)

#include "BspPolygon2.h"
#include "BspTree2.h"
#include "Wm5Memory.h"
using namespace Wm5;

//----------------------------------------------------------------------------
BspPolygon2::BspPolygon2 ()
{
    mTree = 0;
}
//----------------------------------------------------------------------------
BspPolygon2::BspPolygon2 (const BspPolygon2& polygon)
{
    mTree = 0;
    *this = polygon;
}
//----------------------------------------------------------------------------
BspPolygon2::~BspPolygon2 ()
{
    delete0(mTree);
}
//----------------------------------------------------------------------------
BspPolygon2& BspPolygon2::operator= (const BspPolygon2& polygon)
{
    mVMap = polygon.mVMap;
    mVArray = polygon.mVArray;
    mEMap = polygon.mEMap;
    mEArray = polygon.mEArray;
    delete0(mTree);
    mTree = (polygon.mTree ? polygon.mTree->GetCopy() : 0);
    return *this;
}
//----------------------------------------------------------------------------
int BspPolygon2::InsertVertex (const Vector2d& vertex)
{
    VIterator iter = mVMap.find(vertex);
    if (iter != mVMap.end())
    {
        // Vertex already in map, just return its unique index.
        return iter->second;
    }

    // Vertex not in map, insert it and assign it a unique index.
    int i = (int)mVArray.size();
    mVMap.insert(std::make_pair(vertex, i));
    mVArray.push_back(vertex);
    return i;
}
//----------------------------------------------------------------------------
int BspPolygon2::InsertEdge (const Edge2& edge)
{
    assertion(edge.I0 != edge.I1, "Degenerate edges not allowed.\n");

    EIterator iter = mEMap.find(edge);
    if (iter != mEMap.end())
    {
        // Edge already in map, just return its unique index.
        return iter->second;
    }

    // Edge not in map, insert it and assign it a unique index.
    int i = (int)mEArray.size();
    mEMap.insert(std::make_pair(edge, i));
    mEArray.push_back(edge);
    return i;
}
//----------------------------------------------------------------------------
void BspPolygon2::SplitEdge (int v0, int v1, int vmid)
{
    // Find the edge in the map to get the edge-array index.
    EIterator iter = mEMap.find(Edge2(v0, v1));
    assertion(iter != mEMap.end(), "Edge does not exist in the map.\n");
    int eIndex = iter->second;

    // Delete edge <V0,V1>.
    mEMap.erase(iter);

    // Insert edge <V0,VM>.
    mEArray[eIndex].I1 = vmid;
    mEMap.insert(std::make_pair(mEArray[eIndex], eIndex));

    // Insert edge <VM,V1>.
    InsertEdge(Edge2(vmid, v1));
}
//----------------------------------------------------------------------------
void BspPolygon2::Finalize ()
{
    delete0(mTree);
    mTree = new0 BspTree2(*this, mEArray);
}
//----------------------------------------------------------------------------
int BspPolygon2::GetNumVertices () const
{
    return (int)mVMap.size();
}
//----------------------------------------------------------------------------
bool BspPolygon2::GetVertex (int i, Vector2d& vertex) const
{
    if (0 <= i && i < (int)mVArray.size())
    {
        vertex = mVArray[i];
        return true;
    }
    return false;
}
//----------------------------------------------------------------------------
int BspPolygon2::GetNumEdges () const
{
    return (int)mEMap.size();
}
//----------------------------------------------------------------------------
bool BspPolygon2::GetEdge (int i, Edge2& edge) const
{
    if (0 <= i && i < (int)mEArray.size())
    {
        edge = mEArray[i];
        return true;
    }
    return false;
}
//----------------------------------------------------------------------------
void BspPolygon2::GetInsideEdgesFrom (const BspPolygon2& polygon,
    BspPolygon2& inside) const
{
    assertion(mTree != 0, "Tree must exist.\n");

    BspPolygon2 ignore;
    const int numEdges = polygon.GetNumEdges();
    for (int i = 0; i < numEdges; ++i)
    {
        int v0 = polygon.mEArray[i].I0;
        int v1 = polygon.mEArray[i].I1;
        Vector2d vertex0 = polygon.mVArray[v0];
        Vector2d vertex1 = polygon.mVArray[v1];
        mTree->GetPartition(*this, vertex0, vertex1, ignore, inside, inside,
            ignore);
    }
}
//----------------------------------------------------------------------------
BspPolygon2 BspPolygon2::operator~ () const
{
    assertion(mTree != 0, "Tree must exist.\n");

    // negation
    BspPolygon2 neg;
    neg.mVMap = mVMap;
    neg.mVArray = mVArray;
    ECIterator iter = mEMap.begin();
    ECIterator end = mEMap.end();
    for (/**/; iter != end; ++iter)
    {
        neg.InsertEdge(Edge2(iter->first.I1, iter->first.I0));
    }

    neg.mTree = mTree->GetCopy();
    neg.mTree->Negate();
    return neg;
}
//----------------------------------------------------------------------------
BspPolygon2 BspPolygon2::operator& (const BspPolygon2& polygon) const
{
    assertion(mTree != 0, "Tree must exist.\n");

    // intersection
    BspPolygon2 intersect;
    GetInsideEdgesFrom(polygon, intersect);
    polygon.GetInsideEdgesFrom(*this, intersect);
    intersect.Finalize();
    return intersect;
}
//----------------------------------------------------------------------------
BspPolygon2 BspPolygon2::operator| (const BspPolygon2& polygon) const
{
    // union
    const BspPolygon2& thisPolygon = *this;
    return ~(~thisPolygon & ~polygon);
}
//----------------------------------------------------------------------------
BspPolygon2 BspPolygon2::operator- (const BspPolygon2& polygon) const
{
    // difference
    const BspPolygon2& thisPolygon = *this;
    return thisPolygon & ~polygon;
}
//----------------------------------------------------------------------------
BspPolygon2 BspPolygon2::operator^ (const BspPolygon2& polygon) const
{
    // exclusive or
    const BspPolygon2& thisPolygon = *this;
    return (thisPolygon - polygon) | (polygon - thisPolygon);
}
//----------------------------------------------------------------------------
int BspPolygon2::PointLocation (const Vector2d& vertex) const
{
    assertion(mTree != 0, "Tree must exist.\n");
    return mTree->PointLocation(*this, vertex);
}
//----------------------------------------------------------------------------
void BspPolygon2::Print (const char* filename) const
{
    std::ofstream outFile(filename);

    const int numVertices = (int)mVArray.size();
    outFile << "vquantity = " << numVertices << std::endl;
    int i;
    for (i = 0; i < numVertices; ++i)
    {
        outFile << i << "  (" << mVArray[i].X() << ',' << mVArray[i].Y()
            << ')' << std::endl;
    }
    outFile << std::endl;

    const int numEdges = (int)mEArray.size();
    outFile << "equantity = " << numEdges << std::endl;
    for (i = 0; i < numEdges; ++i)
    {
        outFile << "  <" << mEArray[i].I0 << ',' << mEArray[i].I1
            << '>' << std::endl;
    }
    outFile << std::endl;

    outFile << "bsp tree" << std::endl;
    if (mTree)
    {
        mTree->Print(outFile, 0, 'r');
    }
    outFile << std::endl;
}
//----------------------------------------------------------------------------




// Geometric Tools, LLC
// Copyright (c) 1998-2011
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)

#include "Boolean2D.h"

WM5_WINDOW_APPLICATION(Boolean2D);

//----------------------------------------------------------------------------
Boolean2D::Boolean2D ()
    :
    WindowApplication2("SampleMathematics/Boolean2D", 0, 0, 256, 256,
        Float4(1.0f, 1.0f, 1.0f, 1.0f))
{
    mActive = 0;
    mPoly0 = 0;
    mPoly1 = 0;
    mSize = GetWidth();
}
//----------------------------------------------------------------------------
bool Boolean2D::OnInitialize ()
{
    if (!WindowApplication2::OnInitialize())
    {
        return false;
    }

    mChoice = 0;
    mPoly0 = ConstructInvertedEll();
    mPoly1 = ConstructPentagon();
    DoBoolean();

    OnDisplay();
    return true;
}
//----------------------------------------------------------------------------
void Boolean2D::OnTerminate ()
{
    delete0(mPoly0);
    delete0(mPoly1);

    WindowApplication2::OnTerminate();
}
//----------------------------------------------------------------------------
void Boolean2D::OnDisplay ()
{
    ClearScreen();

    DrawPolySolid(*mPoly0, ColorRGB(255, 0, 0));
    DrawPolySolid(*mPoly1, ColorRGB(0, 255, 0));
    if (mActive)
    {
        DrawPolySolid(*mActive, ColorRGB(0, 0, 255));
    }

    WindowApplication2::OnDisplay();
}
//----------------------------------------------------------------------------
bool Boolean2D::OnKeyDown (unsigned char key, int x, int y)
{
    if (WindowApplication2::OnKeyDown(key, x, y))
    {
        return true;
    }

    switch (key)
    {
    case 'n':
    case 'N':
        delete0(mPoly0);
        delete0(mPoly1);
        mActive = 0;

        mChoice = (mChoice + 1) % 3;
        switch (mChoice)
        {
        case 0:
            mPoly0 = ConstructInvertedEll();
            mPoly1 = ConstructPentagon();
            break;
        case 1:
            mPoly0 = ConstructSquare();
            mPoly1 = ConstructSShape();
            break;
        case 2:
            mPoly0 = ConstructPolyWithHoles();
            mPoly1 = ConstructPentagon();
            break;
        }
        DoBoolean();
        break;

    case 'p':
    case 'P':
        mActive = 0;
        break;
    case 'u':
    case 'U':
        mActive = &mUnion;
        break;
    case 'i':
    case 'I':
        mActive = &mIntersection;
        break;
    case 'd':
    case 'D':
        mActive = &mDiff01;
        break;
    case 'e':
    case 'E':
        mActive = &mDiff10;
        break;
    case 'x':
    case 'X':
        mActive = &mXor;
        break;
    }

    OnDisplay();
    return true;
}
//----------------------------------------------------------------------------
BspPolygon2* Boolean2D::ConstructInvertedEll ()
{
    double w = (double)GetWidth();
    double d1d8 = 0.125*w;
    double d2d8 = 0.250*w;
    double d3d8 = 0.375*w;
    double d5d8 = 0.625*w;
    double d6d8 = 0.750*w;
    double d7d8 = 0.875*w;

    const int numVertices = 10;
    Vector2d vertices[numVertices] =
    {
        Vector2d(d1d8, d1d8),
        Vector2d(d3d8, d1d8),
        Vector2d(d3d8, d3d8),
        Vector2d(d2d8, d3d8),
        Vector2d(d2d8, d6d8),
        Vector2d(d5d8, d6d8),
        Vector2d(d5d8, d5d8),
        Vector2d(d7d8, d5d8),
        Vector2d(d7d8, d7d8),
        Vector2d(d1d8, d7d8)
    };

    BspPolygon2* poly = new0 BspPolygon2();
    for (int i0 = numVertices - 1, i1 = 0; i1 < numVertices; i0 = i1++)
    {
        poly->InsertVertex(vertices[i1]);
        poly->InsertEdge(Edge2(i0, i1));
    }
    poly->Finalize();
    return poly;
}
//----------------------------------------------------------------------------
BspPolygon2* Boolean2D::ConstructPentagon ()
{
    const int numVertices = 5;

    double primitiveAngle = Mathd::TWO_PI/numVertices;
    double radius = 0.35*GetWidth();
    double cx = 0.5*GetWidth(), cy = 0.5*GetWidth();

    Vector2d vertices[numVertices];
    for (int i = 0; i < numVertices; ++i)
    {
        double angle = i*primitiveAngle;
        vertices[i].X() = cx + radius*Mathd::Cos(angle);
        vertices[i].Y() = cy + radius*Mathd::Sin(angle);
    }

    BspPolygon2* poly = new0 BspPolygon2();
    for (int i0 = numVertices - 1, i1 = 0; i1 < numVertices; i0 = i1++)
    {
        poly->InsertVertex(vertices[i1]);
        poly->InsertEdge(Edge2(i0, i1));
    }
    poly->Finalize();
    return poly;
}
//----------------------------------------------------------------------------
BspPolygon2* Boolean2D::ConstructSquare ()
{
    double w = (double)GetWidth();
    double d2d8 = 0.250*w;
    double d6d8 = 0.750*w;

    const int numVertices = 4;
    Vector2d vertices[numVertices] =
    {
        Vector2d(d2d8, d2d8),
        Vector2d(d6d8, d2d8),
        Vector2d(d6d8, d6d8),
        Vector2d(d2d8, d6d8)
    };

    BspPolygon2* poly = new0 BspPolygon2();
    for (int i0 = numVertices - 1, i1 = 0; i1 < numVertices; i0 = i1++)
    {
        poly->InsertVertex(vertices[i1]);
        poly->InsertEdge(Edge2(i0, i1));
    }
    poly->Finalize();
    return poly;
}
//----------------------------------------------------------------------------
BspPolygon2* Boolean2D::ConstructSShape ()
{
    double w = (double)GetWidth();
    double d10d32 = 10.0*w/32.0;
    double d12d32 = 12.0*w/32.0;
    double d13d32 = 13.0*w/32.0;
    double d16d32 = 16.0*w/32.0;
    double d19d32 = 19.0*w/32.0;
    double d20d32 = 20.0*w/32.0;
    double d22d32 = 22.0*w/32.0;
    double d24d32 = 24.0*w/32.0;
    double d26d32 = 26.0*w/32.0;
    double d28d32 = 28.0*w/32.0;

    const int numVertices = 12;
    Vector2d vertices[numVertices] =
    {
        Vector2d(d24d32, d10d32),
        Vector2d(d28d32, d10d32),
        Vector2d(d28d32, d16d32),
        Vector2d(d22d32, d16d32),
        Vector2d(d22d32, d19d32),
        Vector2d(d24d32, d19d32),
        Vector2d(d24d32, d22d32),
        Vector2d(d20d32, d22d32),
        Vector2d(d20d32, d13d32),
        Vector2d(d26d32, d13d32),
        Vector2d(d26d32, d12d32),
        Vector2d(d24d32, d12d32)
    };

    BspPolygon2* poly = new0 BspPolygon2();
    for (int i0 = numVertices - 1, i1 = 0; i1 < numVertices; i0 = i1++)
    {
        poly->InsertVertex(vertices[i1]);
        poly->InsertEdge(Edge2(i0, i1));
    }
    poly->Finalize();
    return poly;
}
//----------------------------------------------------------------------------
BspPolygon2* Boolean2D::ConstructPolyWithHoles ()
{
    double w = (double)GetWidth();
    double d2d16 = 2.0*w/16.0;
    double d3d16 = 3.0*w/16.0;
    double d4d16 = 4.0*w/16.0;
    double d6d16 = 6.0*w/16.0;
    double d14d16 = 14.0*w/16.0;

    const int numVertices = 6;
    Vector2d vertices[numVertices] =
    {
        // outer boundary
        Vector2d(d2d16, d2d16),
        Vector2d(d14d16, d2d16),
        Vector2d(d2d16, d14d16),

        // inner boundary
        Vector2d(d4d16, d3d16),
        Vector2d(d6d16, d6d16),
        Vector2d(d6d16, d3d16)
    };

    BspPolygon2* poly = new0 BspPolygon2();
    for (int i = 0; i < numVertices; ++i)
    {
        poly->InsertVertex(vertices[i]);
    }

    poly->InsertEdge(Edge2(0, 1));
    poly->InsertEdge(Edge2(1, 2));
    poly->InsertEdge(Edge2(2, 0));
    poly->InsertEdge(Edge2(3, 4));
    poly->InsertEdge(Edge2(4, 5));
    poly->InsertEdge(Edge2(5, 3));

    poly->Finalize();
    return poly;
}
//----------------------------------------------------------------------------
void Boolean2D::DrawPolySolid (BspPolygon2& polygon, ColorRGB color)
{
    Vector2d  v0, v1;
    Edge2 edge;
    int i, x0, y0, x1, y1;

    // Draw the edges.
    for (i = 0; i < polygon.GetNumEdges(); ++i)
    {
        polygon.GetEdge(i, edge);
        polygon.GetVertex(edge.I0, v0);
        polygon.GetVertex(edge.I1, v1);
        
        x0 = (int)(v0.X() + 0.5f);
        y0 = GetWidth() - 1 - (int)(v0.Y() + 0.5f);
        x1 = (int)(v1.X() + 0.5f);
        y1 = GetWidth() - 1 - (int)(v1.Y() + 0.5f);

        DrawLine(x0, y0, x1, y1, color);
    }

    // Draw the vertices.
    ColorRGB black(0, 0, 0);
    for (i = 0; i < polygon.GetNumVertices(); ++i)
    {
        polygon.GetVertex(i, v0);
        x0 = (int)(v0.X() + 0.5f);
        y0 = GetWidth() - 1 - (int)(v0.Y() + 0.5f);
        SetThickPixel(x0, y0, 1, black);
    }
}
//----------------------------------------------------------------------------
void Boolean2D::DoBoolean ()
{
    BspPolygon2& P = *mPoly0;
    BspPolygon2& Q = *mPoly1;

    mIntersection = P & Q;
    mUnion        = P | Q;
    mDiff01       = P - Q;
    mDiff10       = Q - P;
    mXor          = P ^ Q;
}
//----------------------------------------------------------------------------