package fredholm

import kotlin.test.Test
import kotlin.test.assertEquals

class TestGrid {
    @Test
    fun grid() {
        val a = -1.0
        val b = 5.0
        val size = 20
        val grid = Grid(a, b, size)
        assert(grid[-2] < grid[0])
        assert(grid[-1] < grid[0])
        assert(grid[grid.n + 1] > grid[grid.n])
        assert(grid[grid.n + 2] > grid[grid.n])

        val step = (b - a) / (size - 1)
        for (i in 0..<size) {
            println(i)
            assertEquals(grid[i], a + step * i)
        }
    }
}
