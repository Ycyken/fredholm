package app

import Grid
import functionals.QuadrApprFuncs
import functionals.QuadraticBSplineBasis
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator
import org.jetbrains.kotlinx.multik.api.identity
import org.jetbrains.kotlinx.multik.api.linalg.solve
import org.jetbrains.kotlinx.multik.api.mk
import org.jetbrains.kotlinx.multik.api.ndarray
import org.jetbrains.kotlinx.multik.ndarray.operations.minus
import org.jetbrains.kotlinx.multik.ndarray.operations.toList
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

    val splineBasis = QuadraticBSplineBasis(grid)
    val ws = splineBasis.splines()

    val integrator =
        IterativeLegendreGaussIntegrator(
            gridSize,
            1e-5,
            1e-5,
        )
    // Functions w~_j(t) = Kappa w_j(t) = [integral from a to b (K(t,x) * w_j(x) dx)], j=-2..n-1
    val wsTilda =
        ws.map { w ->
            val intergralKappaWj =
                integrator.integrate(10000, { x -> fredholm.ker.kx(x) * w(x) }, a, b);
            { t: Double -> fredholm.ker.kt(t) * intergralKappaWj }
        }

    val apprFuncs = QuadrApprFuncs(grid)
    // Vector of mu_j_Y(f), j=-2..n-1
    val mu =
        (-2..<grid.n).map { j -> apprFuncs.muApproximation(fredholm.f, j) }.let { mk.ndarray(it) }

    // Matrix M, where M_(i,j) = mu_Y_j(w~_i)
    val M =
        List(grid.n + 2) { j ->
            List(grid.n + 2) { i ->
                apprFuncs.muApproximation(
                    wsTilda[i],
                    j - 2,
                )
            }
        }.let { mk.ndarray(it) }
    val I_minus_M = (mk.identity<Double>(grid.n + 2) - M)
    val c = mk.linalg.solve(I_minus_M, mu)

    val uAppr =
        { x: Double -> fredholm.f(x) + c.toList().zip(wsTilda).sumOf { (c, g) -> c * g(x) } }

    val sinTwice = { x: Double -> 2 * sin(x) }
    println("Approximation error = ${kotlin.math.abs(sinTwice(kotlin.math.PI / 4.0) - uAppr(kotlin.math.PI / 4.0))}")
}

// For simplicity, we assume that the kernel function is separable
data class Kernel(val kt: (Double) -> Double, val kx: (Double) -> Double)

data class Fredholm(val ker: Kernel, val f: (Double) -> Double)
