package fredholm.app

import fredholm.*
import kotlin.math.cos
import kotlin.math.sin

fun main() {
    val a = 0.0
    val b = kotlin.math.PI / 2.0
    val gridSize = 46
    val grid = Grid(a, b, gridSize)

    val kernel = Kernel({ x -> sin(x) }, { x -> cos(x) })
    val f = { x: Double -> sin(x) }

    // u(t) - [integral from a to b (K(t,x) * u(x) dx)] = sint
    val fredholm = Fredholm(kernel, f)

    val uAppr = fredholm.solve(grid, ApprFunctional.MU)

    val original = { x: Double -> 2 * sin(x) }
    val err = approxError(original, uAppr, Grid(a, b, gridSize * 10))
    println(err)
}


