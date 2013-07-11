#ifndef MAC_H
#define MAC_H
#include "grid.h"
#include "traits.h"


template <typename Scalar, int EmbedDim,int FormDim,int WhichForm>
class MACGrid;
//Introduce some more tradition names than the exterior calculus form declarations
template <typename Scalar, int D>
struct grid_types{};
template <typename Scalar>
struct grid_types<Scalar,2>{
    typedef MACGrid<Scalar,2,0,0> NGrid;
    typedef MACGrid<Scalar,2,1,0> UGrid;
    typedef MACGrid<Scalar,2,1,1> VGrid;
    typedef MACGrid<Scalar,2,2,0> CGrid;
};
template <typename Scalar>
struct grid_types<Scalar,3> {
    typedef MACGrid<Scalar,3,0,0> NGrid;
    typedef MACGrid<Scalar,3,1,0> UGrid;
    typedef MACGrid<Scalar,3,1,1> VGrid;
    typedef MACGrid<Scalar,3,1,2> WGrid;
    typedef MACGrid<Scalar,3,2,0> DUGrid;
    typedef MACGrid<Scalar,3,2,1> DVGrid;
    typedef MACGrid<Scalar,3,2,2> DWGrid;
    typedef MACGrid<Scalar,3,3,0> CGrid;
};
//TODO: Make sure grids live in index space

namespace mtao{namespace internal{
template <typename Scalar_, int EmbedDim,int FormDim,int WhichForm>
struct traits<MACGrid<Scalar_,EmbedDim,FormDim,WhichForm> > {
    enum {embed_dim = EmbedDim, form_dim = FormDim, which_form = WhichForm};
    typedef Scalar_ Scalar;
};


//All of the ugly case work for staggered grids...
template <typename GridType>
struct mac_offsets {
};

//NGrid
template <typename Scalar>
struct mac_offsets<MACGrid<Scalar,2,0,0> > {
     constexpr static Scalar offsets[2] = {0,0};
     constexpr static int extra_cells[2] = {1,1};
};
//UGrid
template <typename Scalar>
struct mac_offsets<MACGrid<Scalar,2,1,0> > {
     constexpr static Scalar offsets[2] = {0,0.5};
     constexpr static int extra_cells[2] = {1,0};
};
//VGrid
template <typename Scalar>
struct mac_offsets<MACGrid<Scalar,2,1,1> > {
     constexpr static Scalar offsets[2] = {0.5,0.0};
     constexpr static int extra_cells[2] = {0,1};
};
//CGrid
template <typename Scalar>
struct mac_offsets<MACGrid<Scalar,2,2,0> > {
     constexpr static Scalar offsets[2] = {.5,.5};
     constexpr static int extra_cells[2] = {0,0};
};

//NGrid
template <typename Scalar>
struct mac_offsets<MACGrid<Scalar,3,0,0> > {
     constexpr static Scalar offsets[3] = {0,0,0};
     constexpr static int extra_cells[3] = {1,1,1};
};
//UGrid
template <typename Scalar>
struct mac_offsets<MACGrid<Scalar,3,1,0> > {
     constexpr static Scalar offsets[3] = {0,0.5,0.5};
     constexpr static int extra_cells[3] = {1,0,0};
};
//VGrid
template <typename Scalar>
struct mac_offsets<MACGrid<Scalar,3,1,1> > {
     constexpr static Scalar offsets[3] = {0.5,0.0,0.5};
     constexpr static int extra_cells[3] = {0,1,0};
};
//WGrid
template <typename Scalar>
struct mac_offsets<MACGrid<Scalar,3,1,2> > {
     constexpr static Scalar offsets[3] = {0.5,0.5,0.0};
     constexpr static int extra_cells[3] = {0,0,1};
};
//CGrid
template <typename Scalar>
struct mac_offsets<MACGrid<Scalar,3,3,0> > {
     constexpr static Scalar offsets[3] = {.5,.5,.5};
     constexpr static int extra_cells[3] = {0,0,0};
};
//DUGrid
template <typename Scalar>
struct mac_offsets<MACGrid<Scalar,3,2,0> > {
     constexpr static Scalar offsets[3] = {0.5,0.0,0.0};
     constexpr static int extra_cells[3] = {0,1,1};
};
//DVGrid
template <typename Scalar>
struct mac_offsets<MACGrid<Scalar,3,2,1> > {
     constexpr static Scalar offsets[3] = {0.0,0.5,0.0};
     constexpr static int extra_cells[3] = {1,0,1};
};
//DWGrid
template <typename Scalar>
struct mac_offsets<MACGrid<Scalar,3,2,2> > {
     constexpr static Scalar offsets[3] = {0.0,0.0,0.5};
     constexpr static int extra_cells[3] = {1,1,0};
};
}}
template <typename Scalar, int EmbedDim>
class MACGridFactory;
template <typename Scalar_, int EmbedDim,int FormDim,int WhichForm>
class MACGrid: public GridBase<MACGrid<Scalar_,EmbedDim,FormDim,WhichForm> > {
    typedef MACGrid<Scalar_,EmbedDim,FormDim,WhichForm> MyType;
    public:
    //Make template parameters externally available
    typedef Scalar_ Scalar;
    enum {embed_dim = EmbedDim, form_dim = FormDim, which_form = WhichForm};
    //Some helpful typedefs, not necessarily for external use
    typedef GridBase<MACGrid<Scalar,EmbedDim,FormDim,WhichForm> > Base;
    typedef typename mtao::dim_types<EmbedDim>::template scalar_types<Scalar>::Vec Vec;
    typedef typename mtao::dim_types<EmbedDim>::Veci Veci;

    protected:
    friend class MACGridFactory<Scalar,EmbedDim>;
    MACGrid(const Veci & dim,const Vec & origin,  const Vec & dx): Base(dim,origin,dx) {}
    //using Base::Base;
};
#endif
