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

        if (i == -1 || i == -2) {
            return grid[0] + i * 1e-3
        }
        if (i == n + 1 || i == n + 2) {
            return grid[n] + i * 1e-3
        }
        return grid[i]
    }

    fun toSequence(): Sequence<Double> {
        val arr = doubleArrayOf(grid[0], grid[0]) + grid + doubleArrayOf(grid[n], grid[n])
        return arr.asSequence()
    }
}
