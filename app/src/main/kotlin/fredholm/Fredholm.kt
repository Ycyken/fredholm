package fredholm

import fredholm.functionals.QuadrApprFuncs
import fredholm.functionals.QuadraticBSplineBasis
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator
import org.jetbrains.kotlinx.multik.api.identity
import org.jetbrains.kotlinx.multik.api.linalg.solve
import org.jetbrains.kotlinx.multik.api.mk
import org.jetbrains.kotlinx.multik.api.ndarray
import org.jetbrains.kotlinx.multik.ndarray.operations.minus
import org.jetbrains.kotlinx.multik.ndarray.operations.toList

/*
 * Type of approximation functionals to use.
 */
enum class ApprFunctional {
    GREVILLE,
    MU,
}

/**
 * Solves the Fredholm integral equation of the second kind:
 * u(t) = f(t) + ∫ K(t, x) * u(x) dx
 *
 * The kernel function must take arguments in the order: (t, x)
 */
class Fredholm(val ker: (Double) -> (Double) -> Double, val f: (Double) -> Double) {
    fun solve(
        grid: Grid,
        apprFunc: ApprFunctional,
    ): (Double) -> Double {
        val a = grid[0]
        val b = grid[grid.n]

        val splineBasis = QuadraticBSplineBasis(grid)
        val ws = splineBasis.splines()

        val integrator =
            IterativeLegendreGaussIntegrator(
                grid.n,
                1e-5,
                1e-5,
            )
        // Functions w~_j(t) = Kappa w_j(t) = [integral from a to b (K(t,x) * w_j(x) dx)], j=-2..n-1
        val wsTilda =
            ws.map { w ->
                { t: Double ->
                    integrator.integrate(10000, { x -> ker(t)(x) * w(x) }, a, b)
                }
            }

        val apprFuncs = QuadrApprFuncs(grid)
        // Vector of mu_j_Y(f), j=-2..n-1
        val mu =
            (-2..<grid.n).map { j ->
                when (apprFunc) {
                    ApprFunctional.MU -> apprFuncs.muApproximation(f, j)
                    ApprFunctional.GREVILLE -> apprFuncs.grevilleApproximation(f, j)
                }
            }.let { mk.ndarray(it) }

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
            { x: Double -> f(x) + c.toList().zip(wsTilda).sumOf { (c, g) -> c * g(x) } }
        return uAppr
    }
}
