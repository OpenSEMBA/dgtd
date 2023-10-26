#include "Geometry.h"

namespace maxwell {

using namespace mfem;

bool elementsHaveSameOrientation(const mfem::Element* e1, const mfem::Element* e2)
{
	if (e1->GetType() != e2->GetType()) {
		throw std::runtime_error(
			"Elements with different orientations can not be compared");
	}

	Array<int> e1Vertices;
	e1->GetVertices(e1Vertices);

	Array<int> e2Vertices;
	e2->GetVertices(e2Vertices);

	switch (e1->GetType()) {
	case Element::SEGMENT:
		return e1Vertices == e2Vertices;
		break;
	case Element::TRIANGLE:
	{
		auto it{ std::find(e2Vertices.begin(), e2Vertices.end(), e1Vertices[0]) };
		if (it == e2Vertices.end()) {
			return false;
		}
		std::rotate(e2Vertices.begin(), it, e2Vertices.end());
		return e1Vertices == e2Vertices;
		break;
	}
	default:
		throw std::runtime_error(
			"Unsupported orientation comparison."
		);
	}



}

}