#pragma once

#include <mfem.hpp>

namespace maxwell {
namespace mfemExtension{

using namespace mfem;
class BilinearFormTF : public BilinearForm
{
public:

void Assemble(int skip_zeros = 1);

	BilinearFormTF(FiniteElementSpace* f);

private:
	/// Copy construction is not supported; body is undefined.
	BilinearFormTF(const BilinearFormTF&);

	/// Copy assignment is not supported; body is undefined.
	BilinearFormTF& operator=(const BilinearFormTF&);

};

}
}