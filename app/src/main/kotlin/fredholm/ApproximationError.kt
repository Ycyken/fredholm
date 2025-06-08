package fredholm

import kotlin.math.abs

fun approxError(origin: (Double) -> Double, appr: (Double) -> Double, grid: Grid): Double {
    return grid.toSequence().map { x -> abs(origin(x) - appr(x)) }.max()
}