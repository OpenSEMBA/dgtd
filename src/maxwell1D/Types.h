
namespace maxwell1D {

	enum class FieldType {
		Electric,
		Magnetic
	};

	enum class FluxType {
		Centered,
		Upwind
	};

	enum class BdrCond {
		PEC,
		PMC,
		SMA
	};


}