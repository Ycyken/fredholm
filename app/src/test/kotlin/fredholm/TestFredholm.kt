package fredholm

import kotlin.math.cos
import kotlin.math.sin
import kotlin.test.Test

class TestFredholm {
    @Test
    fun `fredholm with K = sintcosx, f = sint is correct`() {
        val a = 0.0
        val b = kotlin.math.PI / 2.0
        val gridSizes = listOf(16, 31, 46)
        val apprFuncs = listOf(ApprFunctional.GREVILLE, ApprFunctional.MU)

        for (gridSize in gridSizes) {
            for (apprFunc in apprFuncs) {
                val grid = Grid(a, b, gridSize)
                val kernel = { t: Double -> { x: Double -> sin(t) * cos(x) } }
                val f = { x: Double -> sin(x) }
                val fredholm = Fredholm(kernel, f)

                val uAppr = fredholm.solve(grid, apprFunc)
                val original = { x: Double -> 2 * sin(x) }
                val err = approxError(original, uAppr, Grid(a, b, gridSize * 10))
                assert(err < 0.1)
            }
        }
    }
}
