package fredholm

// For simplicity, we assume that the kernel function is separable
class Kernel(val kt: (Double) -> Double, val kx: (Double) -> Double)

class Fredholm(val ker: Kernel, val f: (Double) -> Double) {
}