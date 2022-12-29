namespace SEMBA::cudg3d::dg {

template <class T, std::size_t D>
class Field {
public:
    using Axis = std::size_t;
    using NodeIndex = std::size_t;

    Field() = default;
    Field(std::size_t n);

    T* get(const Axis& d, const NodeIndex& i);
    const T* get(const Axis& d, const NodeIndex& i) const;
    
    void setSize(std::size_t newSize);

private:
    std::array<std::vector<T>, D> val_;
};

using FieldR3 = Field<Math::Real, 3>;

template<class T, std::size_t D>
void Field<T, D>::setSize(std::size_t newSize)
{
    for (auto& v : val_) {
        v.resize(newSize);
    }
}

template<class T, std::size_t D>
Field<T, D>::Field(std::size_t n) {
    for (const auto& v : val_) {
        v = std::vector<T>(n, 0.0);
    }
}

template<class T, std::size_t D>
inline T* Field<T, D>::get(const Axis& d, const std::size_t& i) 
{
    return &val_[d][i];
}

template<class T, std::size_t D>
inline const T* Field<T, D>::get(const Axis& d, const std::size_t& i) const
{
    return &val_[d][i];
}

}