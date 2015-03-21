package sampler

import breeze.numerics._
import breeze.stats.distributions._

import scala.util.control.Breaks._

class DistributionsV1 {
  var cache: String = ""

  def gamrnd(a: Double, b: Double) = Gamma(shape = a, scale = 1 / b).draw()

  def tnorm(mu: Double, variance: Double, a: Double, b: Double): Double = {
    /* Suitable for general purpose truncated normal
     * Based on inverse cdf
     * ! Experiments show this method is not numerical stable.
     */
    val sigma = sqrt(variance)
    val eng = Gaussian(mu, sigma)
    val phia = eng.cdf(a)
    val phib = eng.cdf(b)

    val tmp = phia + Uniform(0.0, 1.0).sample() * (phib - phia)

    cache += "a: " + a.toString + "b: " + b.toString + "phia: " + phia.toString + "phib: " + phib.toString + "tmp: " + tmp.toString + "\n"
    try {
      eng.icdf(tmp)

    } catch {
      case _ => {
        throw new Exception("Fuck off!")
      }
    }
  }

  def tuniform(a: Double, b: Double): Double = {
    Uniform(a, b).draw()
  }

  def halfNorm(mu: Double, variance: Double): Double = {
    if (mu < 0) throw new Exception("Error using half norm!")
    var ans = 0.0


    breakable {
      Gaussian(mu, sqrt(variance)).samples.foreach {
        u => if (u >= 0) {
          ans = u;
          break
        }
      }
    }
    /////debug////////
    if (ans == 0.0) throw new Exception("ans == 0.0")
    ///////////
    ans
  }

  def tailNorm(mu: Double, variance: Double): Double = {
    /* Especially suitable for tail norm in extreme conditions.
     * The algorithm is reject sampling based on exponential distritbuions.
     */
    val lambda: Double = -mu / variance
    var r = 0.0

    // Count the efficiency
    //var rej = 0

    breakable {
      Uniform(0.0, 1.0).samples.foreach {
        u => {
          r = -1 / lambda * log(1 - u)
          if (Uniform(0.0, 1.0).sample <= exp(-1.0 / 2.0 / variance * (r - mu) * (r - mu) + lambda * r + mu * mu / 2.0 / variance))
            break

//          rej += 1
        }
      }
    }
    //println(s"Reject number in tailNorm: $rej")
    r
  }

  def unitUniform(m: Int, n: Int, len: Double = 0.9): Array[Array[Double]] = {
    val ans = Array.ofDim[Double](m + 1, n + 1)
    val sumj = Array.ofDim[Double](n + 1)
    for (i <- 1 to m)
      for (j <- 1 to n) {
        ans(i)(j) = Uniform(0.0, 1.0).draw()
        sumj(j) += ans(i)(j) * ans(i)(j)
      }
    for (i <- 1 to m)
      for (j <- 1 to n) {
        ans(i)(j) /= (sqrt(sumj(j)) / len)
      }
    ans
  }
}

object Distributions extends DistributionsV1 {
  override def tnorm(mu: Double, variance: Double, a: Double, b: Double): Double = {
    /*
     Uses direct reject sampling with uniform proposal distribution
     Specially adjusted for our task!
     */
    val mode = if (mu <= a) (a - mu) * (a - mu)
    else if (mu >= b) (b - mu) * (b - mu)
    else 0
    var r = 0.0

    // Count efficiency
    //var rej = 0

    breakable {
      Uniform(0.0, 1.0).samples.foreach {
        u => {
          r = a + (b - a) * u
          if (Uniform(0.0, 1.0).draw() <= exp(1.0 / 2.0 / variance * (mode - (r - mu) * (r - mu))))
            break
      //    rej += 1
        }
      }
    }
    //println(s"Reject in tnorm: $rej")
    r
  }
}
