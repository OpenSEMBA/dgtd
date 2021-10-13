#include "gtest/gtest.h"

#include "SquareResonantCavity2D.h"

#include <fstream>

class TestDriver: public ::testing::Test {
};

TEST_F(TestDriver, SquareResonantCavity2D) {

	// Posibilidades: las guardo aquí por ahora. Pero voy a tomar la OPCION 3.
	// --> Opcion 1: Comprobación iterativa: Caculo en todos los tiempos y que salte FALSE si en uno no se cumple?
	// --> Opcion 2: Guarda todo: Guardo el campo en un vector con su tiempo y luego lo recorro?
	// --> Opcion 3: Solo un tiempo: Guardo el campo en un vector a tiempo final, total, se arrastra el error.
	// 
	// Medir la Condición Inicial del test. 

	// Calcular la evolución temporal mediante MFEM usando la Condición Inicial del test.

	// Medir la evolución temporal mediante la función analítica.


	// Calcular el error mediante L2. 
	// Solución del test ejemplo: EXPECT_NEAR(error, 0.0, 1e-8);
}