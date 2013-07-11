#ifndef LERP_H
#define LERP_H
namespace mtao{ namespace lerp {
    //Some easy interpolation functions...
    template <typename T>
        T lerp(T v0, T v1, T a) {
            return (1-a) * v0 + a * v1;
        }
    template <typename T>
        T bilerp(T v00, T v01, T v10, T v11, T a, T b) {
            return (1-b) * lerp(v00,v01,a)
                + b * lerp(v10,v11,b);
        }
    template <typename T>
        T trilerp(
                T v000, T v001, T v010, T v011,
                T v100, T v101, T v110, T v111
                , T a, T b, T c) {
            return (1-c) * bilerp(v000,v001,v010,v011,a,b)
                + c * bilerp(v100,v101,v110,v111,a,b);
        }
    //Have to jump some hoops to allow for different dimensioned grids
    template <int Dim, typename GridType>
        struct LinearInterpolator {
        };
    template <typename GridType>
        struct LinearInterpolator<1,GridType> {
            LinearInterpolator(const GridType & g): g(g) {}
            const GridType & g;
            typename GridType::value_type operator()(const typename GridType::Vec & v, const typename GridType::Veci & idx) const {
                const typename GridType::value_type& x=v(0);
                const int& a=idx(0);
                return lerp(g(a),g(a+1),x);
            }
        };
    template <typename GridType>
        struct LinearInterpolator<2,GridType> {
            LinearInterpolator(const GridType & g): g(g) {}
            const GridType & g;
            typename GridType::value_type operator()(const typename GridType::Vec & v, const typename GridType::Veci & idx) const {
                const typename GridType::value_type& x=v(0);
                const typename GridType::value_type& y=v(1);
                const int& a=idx(0);
                const int& b=idx(1);
                return bilerp(g(a,b),g(a+1,b),g(a,b+1),g(a+1,b+1),x,y);
            }
        };
    template <typename GridType>
        struct LinearInterpolator<3,GridType> {
            LinearInterpolator(const GridType & g): g(g) {}
            const GridType & g;
            typename GridType::value_type operator()(const typename GridType::Vec & v, const typename GridType::Veci & idx) const {
                const typename GridType::value_type& x=v(0);
                const typename GridType::value_type& y=v(1);
                const typename GridType::value_type& z=v(2);
                const int& a=idx(0);
                const int& b=idx(1);
                const int& c=idx(2);
                return trilerp(
                        g(a,b,c),g(a+1,b,c),g(a,b+1,c),g(a+1,b+1,c),
                        g(a,b,c+1),g(a+1,b,c+1),g(a,b+1,c+1),g(a+1,b+1,c+1)
                        ,x,y,z);
            }
        };

}
}


#endif
