package functionals

import Grid
import kotlin.math.pow

class QuadraticBSplineBasis(private val grid: Grid) {
    fun splines(): List<(Double) -> Double> {
        return (-2..<grid.n).map { j -> spline(j) }
    }

    private fun spline(j: Int): (Double) -> Double {
        val x_j = grid[j]
        val x_j1 = grid[j + 1]
        val x_j2 = grid[j + 2]
        val x_j3 = grid[j + 3]
        return { x ->
            if (x >= x_j && x < x_j1) {
                (x - x_j).pow(2) / (x_j1 - x_j) / (x_j2 - x_j)
            } else if (x >= x_j1 && x < x_j2) {
                1 / (x_j1 - x_j) * ((x - x_j).pow(2) / (x_j2 - x_j) - (x - x_j1).pow(2) * (x_j3 - x_j) / (x_j2 - x_j1) / (x_j3 - x_j1))
            } else if (x >= x_j2 && x < x_j3) {
                (x - x_j3).pow(2) / (x_j3 - x_j1) / (x_j3 - x_j2)
            } else {
                0.0
            }
        }
    }
}
