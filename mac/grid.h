#ifndef GRID_H
#define GRID_H
#include "types.h"
#include "traits.h"
#include "lerp.h"
#include <vector>

#include <iostream>
template <typename Derived>
class GridBase{
    enum {dim = mtao::internal::traits<Derived>::embed_dim};
    typedef typename mtao::internal::traits<Derived> derived_traits;
    typedef typename derived_traits::Scalar Scalar;
    typedef typename mtao::dim_types<dim> _dt;
    public:
    typedef typename _dt::template scalar_types<Scalar>::Vec Vec;
    typedef typename _dt::Veci Veci;
    typedef typename Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1> > MapVec;
    typedef typename Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1> > ConstMapVec;
    protected:
    GridBase(const Veci & dim,const Vec & origin,  const Vec & dx): m_N(dim),m_origin(origin),  m_dx(dx), m_data(m_N.prod()), m_lerp(*this) {}
    public:

    typedef Scalar value_type;
    int size() const {return stlVector().size();}
    

    //Derived handlers
    Derived& derived(){return *static_cast<Derived*>(this);}
    const Derived & derived()const {return *static_cast<const Derived*>(this);}

    //Internal data accessors
    const Veci& N() const {return m_N;}
    int N(int idx) const {return m_N(idx);}
    const Vec dx() const {return m_dx;}
    Scalar dx(int idx) const {return m_dx(idx);}
    const Vec& origin() const {return m_origin;}
    Scalar origin(int idx) const {return m_origin(idx);}

    //There are several nice ways to access the data...
    std::vector<Scalar> & stlVector(){return m_data;}
    const std::vector<Scalar> & stlVector()const{return m_data;}
    Scalar * data() {return m_data.data();}
    const Scalar * data() const {return m_data.data();}
    MapVec mapVector(){return MapVec(data());}
    ConstMapVec mapVector()const{return ConstMapVec(data());}

    //Some transformations
    Vec worldToIndex(const Vec & v) {return (v-origin()).cwiseQuotient(dx());}
    Vec indexToWorld(const Vec & v) {return v.cwiseProduct(dx()) + origin();}
    Vec indexToWorld(const Veci & v) {return v.template cast<Scalar>().cwiseProduct(dx()) + origin();}
    Vec indexToWorld(int x,int y) {return indexToWorld(Veci(x,y));}
    Vec indexToWorld(int x,int y, int z) {return indexToWorld(Veci(x,y,z));}

    //linear interpolation for a index-space vector
    Scalar ilerp(const Vec & iv) {
        Veci idx = iv.template cast<int>();
        Vec bary = iv - idx.template cast<Scalar>();
        return m_lerp(bary,idx);
    }
    //linear interpolation for a world-space vector
    Scalar lerp(const Vec & v) {
        return ilerp(worldToIndex(v));
    }
    //Maps from 3d to 2d elements
    int idx2ind(int a) const{
        return a;
    }
    int idx2ind(int a, int b) const{
        return m_N(1) * a + b;
    }
    int idx2ind(int a, int b, int c)const {
        return m_N(2) * (m_N(1) * a + b) + c;
    }
    //Accessors
    Scalar& operator()(int a) {return m_data[a];}
    Scalar& operator()(int a, int b) {return m_data[idx2ind(a,b)];}
    Scalar& operator()(int a, int b, int c) {return m_data[idx2ind(a,b,c)];}
    Scalar operator()(int a)const {return m_data[a];}
    Scalar operator()(int a, int b) const{return m_data[idx2ind(a,b)];}
    Scalar operator()(int a, int b, int c)const {return m_data[idx2ind(a,b,c)];}
    //Iterators
    typename std::vector<Scalar>::iterator begin() {return m_data.begin();}
    typename std::vector<Scalar>::iterator end() {return m_data.end();}
    typename std::vector<Scalar>::const_iterator cbegin() {return m_data.cbegin();}
    typename std::vector<Scalar>::const_iterator cend() {return m_data.cend();}

    private:
    Veci m_N;
    Vec m_origin;
    Vec m_dx;
    std::vector<Scalar> m_data;
    const typename mtao::lerp::template LinearInterpolator<Vec::RowsAtCompileTime, GridBase<Derived> > m_lerp;

        
};

#endif
