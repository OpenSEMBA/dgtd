#include "math/vector/Cartesian.h"

namespace SEMBA::dgtd::dg {

template <class T, std::size_t D, std::size_t N>
class Field {
public:

    T* operator()(const std::size_t i);
    const T* operator()(const std::size_t i) const;

    Vector::Cartesian<T, D> getCVec(const std::size_t i) const;

private:
    std::array<T, D* N> val_{ 0 };
};

using FieldR3 = Field<Math::Real, 3>

template<class T, std::size_t D>
Field<T, D>::Field(std::size_t size) {
    size_ = size;
    val_ = new (T) (size_ * D);
}

template<class T, std::size_t D>
inline T* Field<T, D>::operator()(const std::size_t i) {
    assert(i < D);
    return &val_[size_ * i];
}

template<class T, std::size_t D>
inline const T* Field<T, D>::operator()(const std::size_t i) const {
    assert(i < D);
    return &val_[size_ * i];
}

template<class T, std::size_t D>
inline Vector::Cartesian<T, D> Field<T, D>::getCVec(const std::size_t i) const {
    Vector::Cartesian<T, D> res;
    for (std::size_t d = 0; d < D; d++) {
        res(d) = (*this)(d)[i];
    }
    return res;
}

template<class T, std::size_t D>
inline T* Field<T, D>::set(const std::size_t i) const {
    assert(i < D);
    return &val_[size_ * i];
}

template<class T, std::size_t D>
inline T Field<T, D>::operator[](const std::size_t i) const {
    return val_[i];
}

template<class T, std::size_t D>
inline void Field<T, D>::setSize(const std::size_t siz) {
    size_ = siz;
    val_ = new T[D * siz];
}

template<class T, std::size_t D>
inline void Field<T, D>::set(const std::size_t i,
    const Vector::Cartesian<T, D>& vec) {
    for (std::size_t j = 0; j < D; j++) {
        val_[j * size_ + i] = vec(j);
    }
}

template<class T, std::size_t D>
inline void Field<T, D>::set(const std::size_t i, const T& num) {
    for (std::size_t j = 0; j < D; j++) {
        val_[j * size_ + i] = num;
    }
}

template<class T, std::size_t D>
inline void Field<T, D>::setAll(const T& num) {
    for (std::size_t i = 0; i < size_ * D; i++) {
        val_[i] = (T)num;
    }
}

template<class T, std::size_t D>
inline void Field<T, D>::prod(const std::size_t init, const std::size_t end,
    const T param) {
    for (std::size_t d = 0; d < D; d++) {
        for (std::size_t i = init; i < end; i++) {
            val_[d * size_ + i] *= param;
        }
    }
}

template<class T, std::size_t D>
inline void Field<T, D>::copy(const std::size_t init, const std::size_t end,
    const Field<T, D>& field) {
    for (std::size_t d = 0; d < D; d++) {
        for (std::size_t i = init; i < end; i++) {
            val_[d * size_ + i] = field.val_[d * size_ + i];
        }
    }
}

template<class T, std::size_t D>
inline void Field<T, D>::addProd(const std::size_t init, const std::size_t end,
    const Field<T, D>& field, const T param) {
    for (std::size_t d = 0; d < D; d++) {
        for (std::size_t i = init; i < end; i++) {
            val_[d * size_ + i] += field.val_[d * size_ + i] * param;
        }
    }
}

template<class T, std::size_t D>
inline void Field<T, D>::addProd_omp(const std::size_t init,
    const std::size_t end,
    const Field<T, D>& field, const T param) {
    std::size_t i;
    for (std::size_t d = 0; d < D; d++) {
#pragma omp parallel for private(i)
        for (i = init; i < end; i++) {
            val_[d * size_ + i] += field.val_[d * size_ + i] * param;
        }
    }
}

template<class T, std::size_t D>
inline std::size_t Field<T, D>::getDOFs() const {
    return (D * size_);
}

template<class T, std::size_t D>
inline std::size_t Field<T, D>::size() const {
    return size_;
}

template<class T, std::size_t D>
inline void Field<T, D>::swap(Field<T, D>& param,
    const std::size_t first,
    const std::size_t last) {
    for (std::size_t i = 0; i < D; i++) {
        for (std::size_t k = first; k < last; k++) {
            T aux = val_[i * size_ + k];
            val_[i * size_ + k] = param(i)[k];
            param.set(i)[k] = aux;
        }
    }
}

}