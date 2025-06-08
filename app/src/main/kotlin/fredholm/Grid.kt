package fredholm

class Grid(a: Double, b: Double, size: Int) {
    private val grid: DoubleArray
    val n: Int

    init {
        require(size >= 2)
        val step = (b - a) / (size - 1)
        grid =
            DoubleArray(size) { i ->
                a + i * step
            }
        n = grid.size - 1
    }

    operator fun get(i: Int): Double {
        require(i in -2..n + 2)

        // There are two fictional knots on the left and two on the right, with a step size of eps from borders.
        // The same points cannot be used because NaN is returned when evaluating 0/0 in B-splines.
        val eps = 1e-5
        if (i == -1 || i == -2) {
            return grid[0] + i * eps
        }
        if (i == n + 1 || i == n + 2) {
            return grid[n] + i * eps
        }
        return grid[i]
    }

    fun toSequence(): Sequence<Double> {
        val arr = doubleArrayOf(grid[0], grid[0]) + grid + doubleArrayOf(grid[n], grid[n])
        return arr.asSequence()
    }
}
