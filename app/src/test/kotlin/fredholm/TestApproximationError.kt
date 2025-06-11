package fredholm

import kotlin.test.Test

class TestApproximationError {
    @Test
    fun `approximation error is correct`() {
        val f1 = { x: Double -> x * x }
        val f2 = { x: Double -> if (x >= 3) x * x + 0.2 else x * x - 0.2 }
        val err = approxError(f1, f2, Grid(0.0, 5.0, 10))
        assert(err in 0.19..0.21)
    }
}
