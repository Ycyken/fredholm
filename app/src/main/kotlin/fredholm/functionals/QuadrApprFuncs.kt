package fredholm.functionals

import fredholm.Grid

class QuadrApprFuncs(private val grid: Grid) {
    private val theta = 0.5

    private fun yGrid(j: Int): Double {
        require(j in -2..<grid.n)

        if (j == -2) {
            return grid[0]
        }
        if (j in -1..grid.n - 2) {
            return grid[j + 1] + theta * (grid[j + 2] - grid[j + 1])
        }
        return grid[grid.n]
    }

    fun muApproximation(
        f: (Double) -> Double,
        j: Int,
    ): Double {
        require(j in -2..<grid.n)

        if (j == -2) {
            return f(yGrid(-2))
        }
        if (j in -1..grid.n - 2) {
            val t = -1.0 / 8.0 * (f(yGrid(j - 1)) - 10.0 * f(yGrid(j)) + f(yGrid(j + 1)))
            return t
        }
        return f(yGrid(grid.n - 1))
    }

    fun grevilleApproximation(
        f: (Double) -> Double,
        j: Int,
    ): Double {
        val d = 2
        val abscissa = (j + 1..j + d).map { i -> grid[i] }.sum() / d
        return f(abscissa)
    }
}
