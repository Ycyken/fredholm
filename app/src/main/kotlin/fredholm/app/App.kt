package fredholm.app

import fredholm.ApprFunctional
import fredholm.Fredholm
import fredholm.Grid
import fredholm.Kernel
import fredholm.approxError
import kotlin.math.cos
import kotlin.math.sin

fun main() {
    val interval = Pair(0.0, kotlin.math.PI / 2.0)
    val gridSizes = listOf(16, 31, 46)
    val apprFuncs = listOf(ApprFunctional.GREVILLE, ApprFunctional.MU)

    printExperimentTable(gridSizes, apprFuncs, interval)
}

fun printExperimentTable(
    gridSizes: List<Int>,
    apprFuncs: List<ApprFunctional>,
    interval: Pair<Double, Double>,
) {
    val a = interval.first
    val b = interval.second

    print("".padEnd(12))
    gridSizes.forEach { print("| %8s ".format("n=${it - 1}")) }
    println("\n" + "-".repeat(12 + gridSizes.size * 12))

    apprFuncs.forEach { method ->
        print(method.name.padEnd(12))
        gridSizes.forEach { gridSize ->
            val grid = Grid(a, b, gridSize)
            val kernel = Kernel({ x -> sin(x) }, { x -> cos(x) })
            val f = { x: Double -> sin(x) }

            // u(t) - [integral from a to b (K(t,x) * u(x) dx)] = sint
            val fredholm = Fredholm(kernel, f)

            val uAppr = fredholm.solve(grid, method)
            val original = { x: Double -> 2 * sin(x) }
            val err = approxError(original, uAppr, Grid(a, b, gridSize * 10))
            print("| %8.5f ".format(err))
        }
        println()
    }
}
